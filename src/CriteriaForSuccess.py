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
import sys
import multiprocessing as mp
import itertools
from functools import partial
from joblib import Parallel, delayed

from GetInputs import *
from GoldStandardOverlap import *
from nDiffusion.src.run_Diffusion_Class import *
from MGIEnrichment import *
from OddsRatios import *
from Pharmacology import *

import time
startTime = time.time()

def parse_args():
    '''
    Parses CriteriaForSuccess Inputs
    '''
    parser = argparse.ArgumentParser(description = "Criteria For Success Arguments")
    parser.add_argument('--ExperimentName', nargs='?', default = 'Criteria for Success', help = 'Name of disease and/or cohort')
    parser.add_argument('--InputPath', nargs='?', default = './', help = 'path to Big Pipeline Results directory')
    parser.add_argument('--GSPath', nargs='?', default = './', help = 'Path to CSV of Gold Standard Lists')
    parser.add_argument('--CaseControlPath', nargs='?', default = './', help = 'Path to CSV with Sample IDs and 1/0 (Case/Control) - No header')
    parser.add_argument('--ExactTestPath', nargs='?', default = './', help = 'Path to .txt output from ExactTest.sh')
    parser.add_argument('--OutPutPath', nargs='?', default = './', help = 'Path to output directory')
    parser.add_argument('--nDiffusionGraph', nargs = '?', choices = ('STRINGv10', 'STRINGv11', 'MeTEOR', 'toy'), help = 'Network to use for nDiffusion')
    parser.add_argument('--cores', nargs='?', type = int, default = 1, help = 'Number of cores used to run the program')

    return parser.parse_args()

def poolcontext(processes):#*args, **kwargs):
    pool = mp.Pool(processes)
    yield pool
    pool.terminate()

def main(args):
    ## General Workflow 

    # 1) Load BigPipeline Output and Create Consensus Lists
    print("... Loading, Cleaning, and Preparing Big Pipeline Input...\n")
    consensus_fdr01, consensus_fdr001, num_genes = GetInputs(args.InputPath, args.OutPutPath, args.ExperimentName).BigPipeline()

    # Load needed files for the loop
    goldStandards = GetInputs(args.GSPath, None, None).GoldStandards()
    mgi = GetInputs(None, None, None).MGI()
    caseControl = GetInputs(args.CaseControlPath, None, None).CaseControl()
    exactTest = GetInputs(args.ExactTestPath, None, None).ExactTest()
    
    # Loop through each FDR threshold here
    df_dict = {'FDR 0.1': consensus_fdr01, 'FDR 0.01': consensus_fdr001}
    interst_list = ['EPI', 'EAML', 'EAW', 'Reactome', 'STRING', 'Consensus3', 'Consensus2']
    for df_name, df in df_dict.items():

        ### 2) Calculate overlap with known disease gold standards
        print("... Calculating Gold Standard Overlap for " + df_name + "...\n")
        gsOutPutPath = args.OutPutPath + df_name + "/Gold Standard Overlap/"
        os.makedirs(gsOutPutPath, exist_ok= True)
        GoldStandardOverlap(goldStandards, df, df_name, num_genes, interst_list, args.ExperimentName, gsOutPutPath)

        ### 3) nDiffusion
        print("... Running nDiffusion for " + df_name + "...\n")
        # Create output path
        nDiffOutPutPath = args.OutPutPath + df_name + "/nDiffusion/"
        os.makedirs(nDiffOutPutPath, exist_ok = True)

        # Parse network input from args
        network_opts = {'STRINGv10':os.getcwd() + '/nDiffusion/data/networks/STRING_v10.txt',
                        'STRINGv11':os.getcwd() + '/nDiffusion/data/networks/STRING_v11.txt',
                        'MeTEOR':os.getcwd() + '/nDiffusion/data/networks/MeTEOR.txt',
                        'toy':os.getcwd() + '/nDiffusion/data/networks/toy_network.txt'}
        networkPath = network_opts[args.nDiffusionGraph]

        # Prepare input files
        nDiff_pool = []
        pool = mp.Pool(int(args.cores))
        for gL in interst_list:
            for gS in goldStandards.columns:
                gL_hold = df[gL].dropna()
                gS_hold = goldStandards[gS].dropna()

                # Create input text files for nDiffusion input
                #temp_dir = tempfile.TemporaryDirectory()
                temp_gL = tempfile.NamedTemporaryFile(prefix = gL + "_", suffix='.txt', delete=False)
                temp_gL.write(bytes("\n".join(list(gL_hold)), 'utf-8'))
                temp_gL.seek(0)
                temp_gL.close()

                temp_gS = tempfile.NamedTemporaryFile(prefix = gS + "_", suffix = '.txt', delete = False)
                temp_gS.write(bytes("\n".join(list(gS_hold)), 'utf-8'))
                temp_gS.seek(0)
                temp_gS.close()

                # Set new output path
                nDiffOutPutPath = args.OutPutPath + df_name + "/nDiffusion/" + gL + ' v. ' + gS + "/"
                os.makedirs(nDiffOutPutPath, exist_ok=True)

                # Run async pooling for nDiffusion
                pool.apply_async(nDiffusion, args = (networkPath, temp_gL.name, gL, temp_gS.name, gS, nDiffOutPutPath))

        pool.close()
        pool.join()

        ### 4) MGI Enrichment
        print('... Analyzing MGI Enrichment for ' + df_name + '... \n')
        # Create output path
        MGIOutPutPath = args.OutPutPath + df_name + "/MGI Enrichment/"
        os.makedirs(MGIOutPutPath, exist_ok = True)
        # Make function call
        MGIEnrichment(mgi, df, df_name, args.ExperimentName, MGIOutPutPath, interst_list, args.cores)

        ### 5) Odds Ratios
        print('... Calculating Odds Ratios for ' + df_name + '... \n')
        # Create output path
        OROutPutPath = args.OutPutPath + df_name + "/Odds Ratios/"
        os.makedirs(OROutPutPath, exist_ok = True)

        # Make function call
        GetOddsRatios(df, df_name, interst_list, caseControl, exactTest, args.ExperimentName, OROutPutPath, args.cores)

        ### 6) Pharmacology Analysis
        print('... Pulling Drug-Gene Interaction data for ' + df_name + '... \n')
        # Create output path
        PharmaOutPutPath = args.OutPutPath + df_name + "/Pharmacology/"
        os.makedirs(PharmaOutPutPath, exist_ok = True)

        # Make function call
        GetGeneDrugInteractions(df, df_name, interst_list, args.ExperimentName, PharmaOutPutPath, args.cores)


if __name__ == '__main__':
    args = parse_args()
    main(args)
    
    # Analyze time
    executionTime = (time.time() - startTime)
    print('Execution time in minutes: ' + str(executionTime/60))


