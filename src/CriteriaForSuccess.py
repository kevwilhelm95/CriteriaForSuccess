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
from helper_functions import *

# Get starting time for time calculation
import time
startTime = time.time()

def parse_args():
    '''
    Parses CriteriaForSuccess Inputs
    '''
    parser = argparse.ArgumentParser(description = "Criteria For Success Arguments")
    parser.add_argument('--ExperimentName', nargs='?', default = 'Criteria for Success', help = 'Name of disease and/or cohort')
    parser.add_argument('--InputPath', nargs='?', default = None, help = 'path to Big Pipeline Results directory')
    parser.add_argument('--InputList', nargs='?', default = None, help = 'Path to .txt file of gene list')
    parser.add_argument('--VCF', nargs='?', help = 'Path to cohort VCF')
    parser.add_argument('--PickExperiments', nargs='?', default = 'All', help = "No-space, comma-separated list of experiments to run (i.e. GS Overlap,OR)")
    parser.add_argument('--PickNetwork', nargs = '?', choices = ('STRINGv10', 'STRINGv11', 'MeTEOR', 'toy'), help = 'Network to use for nDiffusion')
    parser.add_argument('--GSPath', nargs='?', default = './', help = 'Path to CSV of Gold Standard Lists')
    parser.add_argument('--CaseControlPath', nargs='?', default = './', help = 'Path to CSV with Sample IDs and 1/0 (Case/Control) - No header')
    parser.add_argument('--ExactTestPath', nargs='?', default = './', help = 'Path to .txt output from ExactTest.sh')
    parser.add_argument('--OutPutPath', nargs='?', default = './', help = 'Path to output directory')
    parser.add_argument('--AC_Threshold', nargs='?', type = int, default =5, help = 'Select the Allele Count Threshold to include in Consensus2')
    parser.add_argument('--ref', choices = ('GRCh37', 'GRCh38', 'hg19', 'hg38'), help = 'Gene Location file. Choices = GRCh37, GRCh38, hg19, hg38')
    parser.add_argument('--cores', nargs='?', type = int, default = 1, help = 'Number of cores used to run the program')

    return parser.parse_args()

# Function create intermediate files for and run nDiffusion
def RunnDiffusion(df, df_name, G_main, graph_gene, goldStandards, interst_list, nDiffOutPutPath):
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
        gsOutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/Gold Standard Overlap/')
        GoldStandardOverlap(input_file_dict['Gold Standards'], df, df_name, num_genes, interst_list, arguments.ExperimentName, gsOutPutPath)

    # --- nDiffusion --- #
    if "nDiffusion" in experiments:
        print("... Running nDiffusion for " + df_name + "...\n")
        # Create output path
        nDiffOutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/nDiffusion/')
        RunnDiffusion(df, df_name, input_file_dict['PPI Network'], input_file_dict['PPI Network-GraphGene'], input_file_dict['Gold Standards'], interst_list, nDiffOutPutPath)

    # --- MGI Enrichment --- #
    if "MGI" in experiments:
        print('... Analyzing MGI Enrichment for ' + df_name + '... \n')
        # Create output path
        MGIOutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/MGI Enrichment/')
        # Make function call
        MGIEnrichment(input_file_dict['MGI'], df, df_name, arguments.ExperimentName, MGIOutPutPath, interst_list, arguments.cores)

    # --- Odds Ratios --- #
    if "OR" in experiments:
        print('... Calculating Odds Ratios for ' + df_name + '... \n')
        # Create output path
        OROutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/Odds Ratios/')
        # Make function call
        GetOddsRatios(df, df_name, interst_list, input_file_dict['CaseControl'], arguments.CaseControl, arguments.ref, arguments.VCF, arguments.ExperimentName, OROutPutPath, arguments.cores)

    # --- Pharmacology Analysis --- #
    if "Pharmacology" in experiments:
        print('... Pulling Drug-Gene Interaction data for ' + df_name + '... \n')
        # Create output path
        PharmaOutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/Pharmacology/')
        # Make function call
        GetGeneDrugInteractions(df, df_name, interst_list, arguments.ExperimentName, PharmaOutPutPath, arguments.cores)

def main(args):
    # Determine which experiments to run
    ExpToRun = ParseExperiments(args.PickExperiments)

    # Load required files based on experiments chosen
    inputFileDict = ParseInputFiles(args, ExpToRun)

    # Determine which analysis to run based on declaration of InputPath or InputList
    analysis = ParseInputPaths(args)

    # Load BigPipeline Output and Create Consensus Lists
    if analysis == 'BigPipeline':
        print("... Preparing BigPipeline Input...\n")
        # Load and parse the input files
        consensus_fdr01, consensus_fdr001, num_genes = GetInputs(
            args.InputPath, args.AC_Threshold, args.OutPutPath, args.ExperimentName).BigPipeline(analysis)
        # Define parameters for input data
        df_dict = {'FDR_0.1': consensus_fdr01, 'FDR_0.01': consensus_fdr001}
        interst_list = ['EPI', 'EAML', 'EAW', 'Reactome', 'STRING', 'Consensus3', 'Consensus2']

        # Run criteria for success
        for df_name, df in df_dict.items():
            RunCriteriaForSuccess(df, df_name, interst_list, num_genes, ExpToRun, inputFileDict, args)

    # 1B) Load input gene list and run Criteria for Success
    if analysis == 'InputList':
        print("... Preparing InputList...\n")
        # Get input list
        gL_input = GetInputs(
            args.InputList, args.AC_Threshold, args.OutPutPath, args.ExperimentName).BigPipeline(analysis)
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



