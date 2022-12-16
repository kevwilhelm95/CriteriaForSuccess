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
print('--- STARTING PYCFS---', flush = True)
print('... Loading Packages and Functions ...', flush = True)

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
from PubMed_Enrichment import *
from helper_functions import *
from VariantsBySample import *
from Intermethod_STRING_connectivity_enrichment import *

# Get starting time for time calculation
import time
startTime = time.time()

# Unbuffer pooling commands
#unbuffered = os.fdopen(sys.stdout.fileno(), 'w', 0)
#sys.stdout = unbuffered

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
    parser.add_argument('--PubMedKeywords', nargs = '?', help = 'comma-separated list of key words to query co-mentions for (e.g. "Type 2 Diabetes,Insulin,Obesity")')
    parser.add_argument('--InterConnectivity_Evidences', nargs = '?', help = 'comma-separated list of evidences to score network (e.g. "fusion,coexpression,database")')
    parser.add_argument('--InterConnectivity_Confidence', nargs = '?', choices=('all', 'medium', 'high', 'highest'), help = 'edge combined score confidence level')
    parser.add_argument('--GSPath', nargs='?', default = './', help = 'Path to CSV of Gold Standard Lists')
    parser.add_argument('--CaseControlPath', nargs='?', default = './', help = 'Path to CSV with Sample IDs and 1/0 (Case/Control) - No header')
    parser.add_argument('--OutPutPath', nargs='?', default = './', help = 'Path to output directory')
    parser.add_argument('--AC_Threshold', nargs='?', type = int, default =5, help = 'Select the Allele Count Threshold to include in Consensus2')
    parser.add_argument('--ref', choices = ('GRCh37', 'GRCh38', 'hg19', 'hg38'), help = 'Gene Location file. Choices = GRCh37, GRCh38, hg19, hg38')
    parser.add_argument('--cores', nargs='?', type = int, default = 1, help = 'Number of cores used to run the program')

    return parser.parse_args()


# Skip buffering of Terminal outputs
def Unbuffer():
    class Unbuffered(object):
        def __init__(self, stream):
            self.stream = stream
        def write(self, data):
            self.stream.write(data)
            self.stream.flush()
        def writelines(self,datas):
            self.stream.writelines(datas)
            self.stream.flush()
        def __getattr__(self, attr):
            return getattr(self.stream, attr)
    sys.stdout=Unbuffered(sys.stdout)

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
        GetOddsRatios(df, df_name, interst_list, input_file_dict['CaseControl'], arguments.CaseControlPath, arguments.ref, arguments.VCF, arguments.ExperimentName, OROutPutPath, arguments.cores)

    # --- PubMed Enrichment --- #
    if "PubMed Enrichment" in experiments:
        print(f'... Querying PubMed for Keywords Co-mentioned with {df_name} genes...')
        # Create output path
        pubmedOutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/PubMed Enrichment/')
        PubMed_Enrichment(df, df_name, interst_list, arguments.ref, arguments.PubMedKeywords, pubmedOutPutPath)

    # --- Variants By Sample --- #
    if "Variants By Sample" in experiments:
        print(f'... Parsing VCF to Get Variants By Sample ...')
        # Create output path
        variantsOutPutPath = CreateDir(arguments.OutPutPath, f'IntermediateFiles')
        # Call
        VariantsBySample(df, df_name, arguments.VCF, input_file_dict['CaseControl'], arguments.cores, variantsOutPutPath)

    # --- InterMethodConnectivity --- #
    if "InterMethod Connectivity" in experiments:
        print(f'... Analyzing interconnectivity between methods ...')
        # Create output path 
        connectivityOutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/InterMethod_Connectivity/')
        # Call
        InterMethod_Connectivity(input_file_dict['Gold Standards'], arguments.InterConnectivity_Evidences, arguments.InterConnectivity_Confidence, df, df_name, interst_list, connectivityOutPutPath)

    # --- Phenotype Associations --- #
    if "HOLD_PHENOTYPE_ASSOCIATION" in experiments:
        if ("Variants By Sample" not in experiments) & (arguments.HOLD_VBS_FILE != None):
            x = "Run"

    # --- Pharmacology Analysis --- #
    if "Pharmacology" in experiments:
        print('... Pulling Drug-Gene Interaction data for ' + df_name + '... \n')
        # Create output path
        PharmaOutPutPath = CreateDir(arguments.OutPutPath, f'{df_name}/Pharmacology/')
        # Make function call
        GetGeneDrugInteractions(df, df_name, interst_list, arguments.ExperimentName, PharmaOutPutPath, arguments.cores)

def main(args):
    Unbuffer()
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



