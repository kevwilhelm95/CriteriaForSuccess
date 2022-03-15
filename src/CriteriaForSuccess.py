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

from GetInputs import *
from GoldStandardOverlap import *
from nDiffusion.src.run_Diffusion_Class import *
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
    parser.add_argument('--Analysis', nargs='?', choices=('BigPipeline', 'InputList'))
    parser.add_argument('--GSPath', nargs='?', default = './', help = 'Path to CSV of Gold Standard Lists')
    parser.add_argument('--CaseControlPath', nargs='?', default = './', help = 'Path to CSV with Sample IDs and 1/0 (Case/Control) - No header')
    parser.add_argument('--ExactTestPath', nargs='?', default = './', help = 'Path to .txt output from ExactTest.sh')
    parser.add_argument('--OutPutPath', nargs='?', default = './', help = 'Path to output directory')
    parser.add_argument('--nDiffusionGraph', nargs = '?', choices = ('STRINGv10', 'STRINGv11', 'MeTEOR', 'toy'), help = 'Network to use for nDiffusion')
    parser.add_argument('--AC_Threshold', nargs='?', type = int, default =5, help = 'Select the Allele Count Threshold to include in Consensus2')
    parser.add_argument('--cores', nargs='?', type = int, default = 1, help = 'Number of cores used to run the program')

    return parser.parse_args()


# Function create intermediate files for and run nDiffusion
def RunnDiffusion(df, df_name, goldStandards, interst_list, nDiffOutPutPath):
    # Parse network input from args
    network_opts = {'STRINGv10':os.getcwd() + '/nDiffusion/data/networks/STRING_v10.txt',
                        'STRINGv11':os.getcwd() + '/nDiffusion/data/networks/STRING_v11.txt',
                        'MeTEOR':os.getcwd() + '/nDiffusion/data/networks/MeTEOR.txt',
                        'toy':os.getcwd() + '/nDiffusion/data/networks/toy_network.txt'}
    networkPath = network_opts[args.nDiffusionGraph]

    #Prepare input files
    pool = mp.Pool(int(args.cores))
    for gL in interst_list:
        for gS in goldStandards.columns:
            gL_hold = df[gL].dropna()
            gS_hold = goldStandards[gS].dropna()

            # Create input text files for nDiffusion input
            temp_gL = tempfile.NamedTemporaryFile(
                    prefix=gL + "_", suffix='.txt', delete=False)
            temp_gL.write(bytes("\n".join(list(gL_hold)), 'utf-8'))
            temp_gL.seek(0)
            temp_gL.close()

            temp_gS = tempfile.NamedTemporaryFile(
                    prefix=gS + "_", suffix='.txt', delete=False)
            temp_gS.write(bytes("\n".join(list(gS_hold)), 'utf-8'))
            temp_gS.seek(0)
            temp_gS.close()

            # Set new output path
            nDiffOutPutPath = args.OutPutPath + df_name + \
                    "/nDiffusion/" + gL + ' v. ' + gS + "/"
            os.makedirs(nDiffOutPutPath, exist_ok=True)

            # Run async pooling for nDiffusion
            pool.apply(nDiffusion, args=(
                    networkPath, temp_gL.name, gL, temp_gS.name, gS, nDiffOutPutPath))

    pool.close()
    pool.join()
    pool.terminate()

# Function holding calls to external scripts
def RunCriteriaForSuccess(df, df_name, interst_list, num_genes, goldStandards, mgi, caseControl, exactTest, experimentName, outPath, cores):

    # --- 2) Calculate Overlap with Known Disease Gold Standards --- #
    print("... Calculating Gold Standard Overlap for " + df_name + "...\n")
    # Create new directory for experiment
    gsOutPutPath = outPath + df_name + "/Gold Standard Overlap/"
    os.makedirs(gsOutPutPath, exist_ok= True)
    GoldStandardOverlap(goldStandards, df, df_name, num_genes, interst_list, experimentName, gsOutPutPath)


    # --- 3) nDiffusion --- #
    print("... Running nDiffusion for " + df_name + "...\n")
    # Create output path
    nDiffOutPutPath = args.OutPutPath + df_name + "/nDiffusion/"
    os.makedirs(nDiffOutPutPath, exist_ok = True)
    RunnDiffusion(df, df_name, goldStandards, interst_list, nDiffOutPutPath)


    # --- 4) MGI Enrichment --- #
    print('... Analyzing MGI Enrichment for ' + df_name + '... \n')
    # Create output path
    MGIOutPutPath = args.OutPutPath + df_name + "/MGI Enrichment/"
    os.makedirs(MGIOutPutPath, exist_ok = True)
    # Make function call
    MGIEnrichment(mgi, df, df_name, args.ExperimentName, MGIOutPutPath, interst_list, cores)


    # --- 5) Odds Ratios --- #
    ### 5) Odds Ratios
    print('... Calculating Odds Ratios for ' + df_name + '... \n')
    # Create output path
    OROutPutPath = args.OutPutPath + df_name + "/Odds Ratios/"
    os.makedirs(OROutPutPath, exist_ok = True)
    # Make function call
    GetOddsRatios(df, df_name, interst_list, caseControl, exactTest, experimentName, OROutPutPath, cores)


    # --- 6) Pharmacology Analysis --- #
    print('... Pulling Drug-Gene Interaction data for ' + df_name + '... \n')
    # Create output path
    PharmaOutPutPath = args.OutPutPath + df_name + "/Pharmacology/"
    os.makedirs(PharmaOutPutPath, exist_ok = True)
    # Make function call
    GetGeneDrugInteractions(df, df_name, interst_list, experimentName, PharmaOutPutPath, cores)

def main(args):
    ## General Workflow 
    # Load required external files
    print("... Loading Secondary Files ... \n")
    goldStandards = GetInputs(args.GSPath, None, None, None).GoldStandards()
    mgi = GetInputs(None, None, None, None).MGI()
    caseControl = GetInputs(args.CaseControlPath, None,None, None).CaseControl()
    exactTest = GetInputs(args.ExactTestPath, None, None, None).ExactTest()

    # 1A) Load BigPipeline Output and Create Consensus Lists
    print("... Loading, Cleaning, and Preparing Big Pipeline Input...\n")
    if args.Analysis == 'BigPipeline':
        consensus_fdr01, consensus_fdr001, num_genes = GetInputs(
            args.InputPath, args.AC_Threshold, args.OutPutPath, args.ExperimentName).BigPipeline(args.Analysis)
        df_dict = {'FDR 0.1': consensus_fdr01, 'FDR 0.01': consensus_fdr001}
        interst_list = ['EPI', 'EAML', 'EAW', 'Reactome', 'STRING', 'Consensus3', 'Consensus2', 'UniqLinked']

        # Run criteria for success
        for df_name, df in df_dict.items():
            RunCriteriaForSuccess(df, df_name, interst_list, num_genes, goldStandards, mgi, caseControl, exactTest, args.ExperimentName, args.OutPutPath, args.cores)

    # 1B) Load input gene list and run Criteria for Success
    if args.Analysis == 'InputList':
        # Get input list
        gL_input = GetInputs(
            args.InputList, args.AC_Threshold, args.OutPutPath, args.ExperimentName).BigPipeline(args.Analysis)
        # Set number of discoverable genes - Currently from MS
        num_genes = 19872

        # Run criteria for success
        RunCriteriaForSuccess(gL_input, 'Input-List', gL_input.columns, num_genes, goldStandards, mgi, caseControl, exactTest, args.ExperimentName, args.OutPutPath, args.cores)
    


if __name__ == '__main__':
    args = parse_args()
    main(args)
    
    # Analyze time
    executionTime = (time.time() - startTime)
    print('Execution time in minutes: ' + str(executionTime/60))


