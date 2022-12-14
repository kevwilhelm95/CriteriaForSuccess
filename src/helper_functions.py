import os
import pandas as pd
from GetInputs import *


# Test which experiments to run
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

# Load the PPI Network file
def ParseNetwork(network_str):
    """Function to pick and load the network of choice

    Args:
        network_str (str): String detailing which network to use - loaded using argument "--PickNetwork"
    
    Returns:
        type????: network file containing node1, node2, and weight of the network to use
    """
    network_opts = {'STRINGv10': os.path.dirname(os.getcwd()) + '/refs/STRINGv10.txt',
                   'STRINGv11': os.path.dirname(os.getcwd()) + '/refs/STRINGv11.txt',
                   'MeTEOR': os.path.dirname(os.getcwd()) + '/refs/MeTEOR.txt',
                   'toy': os.path.dirname(os.getcwd()) + '/refs/toy_network.txt'}
    networkPath = network_opts[network_str]
    G_main = nx.read_edgelist(open(networkPath, 'r'), data=(('weight', float),))
    
    graph_gene = []
    if network_str == 'MeTEOR':
        for line in open(networkPath).readlines():
            line = line.strip('\n').split('\t')
            for i in line[:2]:
                if ':' not in i:
                    graph_gene.append(i)

    return G_main, graph_gene

# Test which files we need based on the experiments to be run
def ParseInputFiles(arguments, experiments_lst):
    """Load files based on experiments needing to be run

    Args:
        arguments (???): Arguments object
        experiments_lst (list): Parsed experiments to run list

    Returns:
        dict: dictionary of pd.DataFrames containing secondary data
    """
    fileDict = {}

    if any(check in ['GS Overlap', 'nDiffusion'] for check in experiments_lst):
        fileDict['Gold Standards'] = GetInputs(arguments.GSPath, None, None, None).GoldStandards()
    if any(check in ['nDiffusion'] for check in experiments_lst):
        fileDict['PPI Network'], fileDict['PPI Network-GraphGene'] = ParseNetwork(arguments.PickNetwork)
    if any(check in ['MGI'] for check in experiments_lst):
        fileDict['MGI'] = GetInputs(None, None, None, None).MGI()
    if any(check in ['OR'] for check in experiments_lst):
        fileDict['CaseControl'] = GetInputs(arguments.CaseControlPath, None,None, None).CaseControl()
    return fileDict

# Test which analysis to run based on declaration of InputPath or InputList
def ParseInputPaths(arguments):
    if arguments.InputPath == None:
        analysis = "InputList"
    elif arguments.InputList == None:
        analysis = "BigPipeline"
    return analysis

def ParseGeneLocationFile(ref):
    ref_opts = {'GRCh37': os.path.dirname(os.getcwd()) + '/refs/ENSEMBL-lite_GRCh37.v75.txt',
                'GRCh38': os.path.dirname(os.getcwd()) + '/refs/ENSEMBL-lite_GRCh38.v94.txt',
                'hg19': os.path.dirname(os.getcwd()) + '/refs/refGene-lite_hg19.May2013.txt',
                'hg38' : os.path.dirname(os.getcwd()) + '/refs/refGene-lite_hg38.June2017.txt'}
    ref_path = ref_opts[ref]
    ref_df = pd.read_csv(ref_path, sep='\t')
    return ref_df

# Create new directories for experiments
def CreateDir(output_path, dir_name):
    hold_outpath = f"{output_path}/{dir_name}/"
    os.makedirs(hold_outpath, exist_ok = True)
    return hold_outpath

# Create intermediate file containing only sample IDs for ExactTest
def CreateSampleOnlyFile(CaseControl_file, output_path):
    samp_df = CaseControl_file.iloc[:,0]
    samp_df.to_csv(output_path + "CaseControl_SampleOnly.txt", sep = '\t', index=False, header=False)

def CreateSampleCaseTabFile(CaseControl_file, output_path):
    CaseControl_file.to_csv(output_path + "CaseControl_SampleCase.txt", sep='\t', index = False, header = False)

# Create intermediate .fam file for ExactTest
def CreateSampleFamFile(CaseControl_file, output_path):
    samp_fam_df = pd.DataFrame({'famid':CaseControl_file.iloc[:,0],
                                'sampid':CaseControl_file.iloc[:,0],
                                'fatherid':0,
                                'motherid':0,
                                'sex':-9,
                                'disease':CaseControl_file.iloc[:,1]})
    samp_fam_df.to_csv(output_path + "CaseControl_fam.fam", sep='\t', index=False, header=False)

def CreateGeneRegionFile(series, ref, output_path):
    ref_df = ParseGeneLocationFile(ref)
    ref_sub_df = ref_df[ref_df.gene.isin(series)]
    ref_sub_df = ref_sub_df[['chrom', 'start', 'end']]
    ref_sub_df.to_csv(output_path + "AllUniqueGenesLocationFile.txt", sep = '\t', index = False, header = False)