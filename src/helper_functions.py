import os
import pandas as pd

def CreateDir(output_path, dir_name):
    hold_outpath = f"{output_path}/{dir_name}/"
    os.makedirs(hold_outpath, exist_ok = True)
    return hold_outpath

def CreateSampleOnlyFile(CaseControl_file, output_path):
    samp_df = CaseControl_file.iloc[:,0]
    samp_df.to_csv(output_path + "IntermediateFiles/CaseControl_SampleOnly.txt", sep = '\t', index=False, header=False)

def CreateSampleFamFile(CaseControl_file, output_path):
    samp_fam_df = pd.DataFrame({'famid':CaseControl_file.iloc[:,0],
                                'sampid':CaseControl_file.iloc[:,0]},
                                'fatherid':0,
                                'motherid':0,
                                'sex':-9,
                                'disease':CaseControl_file.iloc[:,1])
    samp_fam_df.to_csv(output_path + "IntermediateFiles/CaseControl_fam.fam", sep='\t', index=False, header=False)