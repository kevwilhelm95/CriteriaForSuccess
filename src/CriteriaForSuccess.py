'''
@author: Kevin Wilhelm
Currently missing:
    - STRING Enrichment
    - FUMA Enrichment
    - Pathway Odds Ratios
    - Structures

Notes:
- Run from CriteriaForSuccess/src directory

-Files need in CriteriaForSuccess/refs
    - HMD_HumanPhenotype_v2021.rpt (MGI High-Level Enrichment)
'''


import pandas as pd
import os
import tempfile
import argparse
from datetime import date
import multiprocessing as mp
import networkx as nx

from GetInputs import *
from GoldStandardOverlap import *
from nDiffusion.src.run_Diffusion_Class_v2 import *
from MGIEnrichment import *
from OddsRatios import *
from Pharmacology import *

# Get starting time for time calculation
import time
startTime = time.time()

def parse_args():
    '''
    Parses CriteriaForSuccess Inputs
    '''
    parser = argparse.ArgumentParser(description = "Criteria For Success Arguments")
    parser.add_argument('--ExperimentName', nargs='?', default = 'Criteria for Success', help = 'Name of disease and/or cohort')
    parser.add_argument('--InputPath', nargs='?', default = './', help = 'path to Big Pipeline Results directory')
    parser.add_argument('--InputList', nargs='?', default = './', help = 'Path to .txt file of gene list')
    parser.add_argument('--Analysis', nargs='?', choices=('BigPipeline', 'InputList')) ### Can hard code this in based on if values above are empty or not
    parser.add_argument('--PickExperiments', nargs='?', default = 'All', help = "No-space, comma-separated list of experiments to run (i.e. GS Overlap,OR)")
    parser.add_argument('--PickNetwork', nargs = '?', choices = ('STRINGv10', 'STRINGv11', 'MeTEOR', 'toy'), help = 'Network to use for nDiffusion')
    parser.add_argument('--GSPath', nargs='?', default = './', help = 'Path to CSV of Gold Standard Lists')
    parser.add_argument('--CaseControlPath', nargs='?', default = './', help = 'Path to CSV with Sample IDs and 1/0 (Case/Control) - No header')
    parser.add_argument('--ExactTestPath', nargs='?', default = './', help = 'Path to .txt output from ExactTest.sh')
    parser.add_argument('--OutPutPath', nargs='?', default = './', help = 'Path to output directory')
    parser.add_argument('--AC_Threshold', nargs='?', type = int, default =5, help = 'Select the Allele Count Threshold to include in Consensus2')
    parser.add_argument('--cores', nargs='?', type = int, default = 1, help = 'Number of cores used to run the program')

    return parser.parse_args()

def ParseExperiments(experiments_str):
    """Function to parse the experiments to run

    Args:
        experiments_str (str): comma-separated, no space sequence of experiments to run - All will be filled in

    Returns:
        list: parsed list of experiments needing to be run
    """
    experiments_lst = experiments_str.split(",")
    if 'All' in experiments_lst:
        experiments_lst = ['GS Overlap', 'nDiffusion', 'MGI', 'OR', 'Pharmacology']
    print("Experiments to run: ", experiments_lst)
    return experiments_lst

def ParseNetwork(network_str):
    """Function to pick and load the network of choice

    Args:
        network_str (str): String detailing which network to use - loaded using argument "--PickNetwork"
    
    Returns:
        type????: network file containing node1, node2, and weight of the network to use
    """
    network_opts = {'STRINGv10': os.getcwd() + '/nDiffusion/data/networks/STRING_v10.txt',
                   'STRINGv11': os.getcwd() + '/nDiffusion/data/networks/STRING_v11.txt',
                   'MeTEOR': os.getcwd() + '/nDiffusion/data/networks/MeTEOR.txt',
                   'toy': os.getcwd() + '/nDiffusion/data/networks/toy_network.txt'}
    networkPath = network_opts[network_str]
    G_main = nx.read_edgelist(open(networkPath, 'r'), data=(('weight', float),))
    return G_main

def ParseInputFiles(arguments, experiments_lst):
    """Load files based on experiments needing to be run

    Args:
        arguments (???): Arguments object
        experiments_lst (list): Parsed experiments to run list

    Returns:
        pd.DataFrame: Dataframes of supplementary data needed for experiments
    """
    fileDict = {'Gold Standards' : None,
                'PPI Network' : None,
                'MGI' : None,
                'CaseControl' : None,
                'ExactTest' : None}

    if any(check in ['GS Overlap', 'nDiffusion'] for check in experiments_lst):
        fileDict['Gold Standards'] = GetInputs(arguments.GSPath, None, None, None).GoldStandards()
    if any(check in ['nDiffusion'] for check in experiments_lst):
        ppi_network = ParseNetwork(arguments.PickNetwork)
    if any(check in ['MGI'] for check in experiments_lst):
        mgi = GetInputs(None, None, None, None).MGI()
    if any(check in ['OR'] for check in experiments_lst):
        caseControl = GetInputs(args.CaseControlPath, None,None, None).CaseControl()
        exactTest = GetInputs(args.ExactTestPath, None, None, None).ExactTest()
    return fileDict

# Function create intermediate files for and run nDiffusion
def RunnDiffusion(df, df_name, G_main, goldStandards, interst_list, nDiffOutPutPath):
    # Parse network input from args
    network_opts = {'STRINGv10':os.getcwd() + '/nDiffusion/data/networks/STRING_v10.txt',
                        'STRINGv11':os.getcwd() + '/nDiffusion/data/networks/STRINGv11.txt',
                        'MeTEOR':os.getcwd() + '/nDiffusion/data/networks/MeTEOR.txt',
                        'toy':os.getcwd() + '/nDiffusion/data/networks/toy_network.txt'}
    networkPath = network_opts[args.nDiffusionGraph]

    # Get the genes for multimodal networks
    graph_gene = []
    if args.nDiffusionGraph == 'MeTeOR':
        for line in open(networkPath).readlines():
            line = line.strip('\n').split('\t')
            for i in line[:2]:
                if ':' not in i:
                    graph_gene.append(i)

    #Prepare input files
    processes = []

    #pool = mp.Pool(int(args.cores))
    for gL in interst_list:
        for gS in goldStandards.columns:
            gL_hold = df[gL].dropna()
            gS_hold = goldStandards[gS].dropna()

            # Set new output path
            nDiffOutPutPath = args.OutPutPath + df_name + \
                    "/nDiffusion/" + gL + ' v. ' + gS + "/"
            os.makedirs(nDiffOutPutPath, exist_ok=True)

            print('registering process %s' % gS)
            processes.append(Process(target = nDiffusion, args = ( (G_main, graph_gene, gL_hold, gL, gS_hold, gS, nDiffOutPutPath) )))

    for process in processes:
        process.start()

    for process in processes:
        process.join()

# Function holding calls to external scripts
def RunCriteriaForSuccess(df, df_name, interst_list, num_genes, experiments, input_file_dict, arguments):
    # Make note of experimentName, outPath, cores
    # --- Calculate Overlap with Known Disease Gold Standards --- #
    if "GS Overlap" in experiments:
        print("... Calculating Gold Standard Overlap for " + df_name + "...\n")
        # Create new directory for experiment
        gsOutPutPath = arguments.OutPutPath + df_name + "/Gold Standard Overlap/"
        os.makedirs(gsOutPutPath, exist_ok= True)
        GoldStandardOverlap(input_file_dict['Gold Standards'], df, df_name, num_genes, interst_list, arguments.ExperimentName, gsOutPutPath)

    # --- nDiffusion --- #
    if "nDiffusion" in experiments:
        print("... Running nDiffusion for " + df_name + "...\n")
        # Create output path
        nDiffOutPutPath = arguments.OutPutPath + df_name + "/nDiffusion/"
        os.makedirs(nDiffOutPutPath, exist_ok = True)
        RunnDiffusion(df, df_name, input_file_dict['PPI Network'], input_file_dict['Gold Standards'], interst_list, nDiffOutPutPath)

    # --- MGI Enrichment --- #
    if "MGI" in experiments:
        print('... Analyzing MGI Enrichment for ' + df_name + '... \n')
        # Create output path
        MGIOutPutPath = arguments.OutPutPath + df_name + "/MGI Enrichment/"
        os.makedirs(MGIOutPutPath, exist_ok = True)
        # Make function call
        MGIEnrichment(input_file_dict['MGI'], df, df_name, arguments.ExperimentName, MGIOutPutPath, interst_list, arguments.cores)

    # --- Odds Ratios --- #
    if "OR" in experiments:
        print('... Calculating Odds Ratios for ' + df_name + '... \n')
        # Create output path
        OROutPutPath = arguments.OutPutPath + df_name + "/Odds Ratios/"
        os.makedirs(OROutPutPath, exist_ok = True)
        # Make function call
        GetOddsRatios(df, df_name, interst_list, input_file_dict['CaseControl'], input_file_dict['ExactTest'], arguments.ExperimentName, OROutPutPath, arguments.cores)

    # --- Pharmacology Analysis --- #
    if "Pharmacology" in experiments:
        print('... Pulling Drug-Gene Interaction data for ' + df_name + '... \n')
        # Create output path
        PharmaOutPutPath = arguments.OutPutPath + df_name + "/Pharmacology/"
        os.makedirs(PharmaOutPutPath, exist_ok = True)
        # Make function call
        GetGeneDrugInteractions(df, df_name, interst_list, arguments.ExperimentName, PharmaOutPutPath, arguments.cores)

def main(args):
    # Determine which experiments to run
    ExpToRun = ParseExperiments(args.PickExperiments)

    # Load required files based on experiments chosen
    inputFileDict = ParseInputFiles(args, ExpToRun)
    print(type(inputFileDict['Gold Standards']))

    # Load BigPipeline Output and Create Consensus Lists
    print("... Loading, Cleaning, and Preparing Big Pipeline Input...\n")
    if args.Analysis == 'BigPipeline':
        # Load and parse the input files
        consensus_fdr01, consensus_fdr001, num_genes = GetInputs(
            args.InputPath, args.AC_Threshold, args.OutPutPath, args.ExperimentName).BigPipeline(args.Analysis)
        # Define parameters for input data
        df_dict = {'FDR_0.1': consensus_fdr01, 'FDR_0.01': consensus_fdr001}
        interst_list = ['EPI', 'EAML', 'EAW', 'Reactome', 'STRING', 'Consensus3', 'Consensus2']

        # Run criteria for success
        for df_name, df in df_dict.items():
            RunCriteriaForSuccess(df, df_name, interst_list, num_genes, ExpToRun, inputFileDict, args)

    # 1B) Load input gene list and run Criteria for Success
    if args.Analysis == 'InputList':
        # Get input list
        gL_input = GetInputs(
            args.InputList, args.AC_Threshold, args.OutPutPath, args.ExperimentName).BigPipeline(args.Analysis)
        # Set number of discoverable genes - Currently from MS
        num_genes = 19872

        # Run criteria for success
        RunCriteriaForSuccess(gL_input, 'Input-List', gL_input.columns, num_genes, ExpToRun, inputFileDict, args)
    


if __name__ == '__main__':
    args = parse_args()
    main(args)
    
    # Analyze time
    executionTime = (time.time() - startTime)
    print('Execution time in minutes: ' + str(executionTime/60))



