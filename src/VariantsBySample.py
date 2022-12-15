import pandas as pd
import numpy as np
from pysam import VariantFile
from joblib import Parallel, delayed
from tqdm import tqdm
import time


class VariantsBySample():
    def __init__(self, df, df_name, VCF_path, CaseControl, cores, outpath):
        self.df = df
        self.df_name = df_name
        self.VCF_path = VCF_path
        self.CaseControl = CaseControl
        self.cores = cores
        self.outpath = outpath
        self.main()

    def GetGenes(self):
        if self.df_name == 'Input-List':
            genes = self.df['Genes']
        else:
            genes = self.df['AllUnique']
        return genes

    # Split samples by cases and controls
    def parse_samples(samples):
        cases = samples[samples.iloc[:,0] == 1].index.astype(str).tolist()
        conts = samples[samples.iloc[:,0] == 0].index.astype(str).tolist()
        return cases,conts

    # Read through vcf and fetch columns of interest
    def parse_VEP(vcf_fn, genes, sample, min_af = None, max_af = None, af_field = 'AF', EA_parser = 'canonical'):
        row=[]   
        sample=[sample]
        vcf = VariantFile(vcf_fn)
        vcf.subset_samples(sample)
            
        for rec in vcf.fetch():
            zyg = convert_zygo(rec.samples[sample[0]]['GT'])
            if (zyg!=0) and (rec.info['SYMBOL'][0] in list(genes)):
                all_ea = rec.info.get('EA', (None,))
                all_ensp = rec.info.get('Ensembl_proteinid', (rec.info['ENSP'][0],))
                ea = fetch_EA_VEP(all_ea, rec.info['ENSP'][0], all_ensp, rec.info['Consequence'][0], EA_parser='canonical')
                if not np.isnan(ea):
                    row.append([rec.chrom, rec.pos, rec.ref, rec.alts[0], ea, rec.info['SYMBOL'][0], sample[0],zyg, rec.info['AF'][0] ])
        cols = ['chr','pos','ref','alt','EA','gene','sample','zyg','AF']
        col_type = {'chr': str, 'pos': str, 'ref': str, 'alt': str}
        df = pd.DataFrame(row, columns = cols)
        df = df.astype(col_type)
        return df
    

    def Parse_Variants(self, args_vcf, args_samples, genes, args_cores):
        cases, conts = self.parse_samples(args_samples)
        conts_dfs = Parallel(n_jobs=int(args_cores))(delayed(self.parse_VEP)(args_vcf, genes, sample, min_af=None, max_af=None, af_field='AF', EA_parser='canonical') for sample in tqdm(conts))
        conts_var = pd.concat(conts_dfs, axis=0)
        conts_var.reset_index(drop = True, inplace = True)
                            
        cases_dfs = Parallel(n_jobs=int(args_cores))(delayed(parse_VEP)(args_vcf, genes, sample, min_af=None, max_af=None, af_field='AF', EA_parser='canonical') for sample in tqdm(cases))
        cases_var = pd.concat(cases_dfs, axis=0)
        cases_var.reset_index(drop = True, inplace = True)

        return cases, conts, cases_var, conts_var   


    def main(self):
        # Get genes
        genes = self.GetGenes()

        # Parse variants for outputting
        cases, controls, case_variants, control_variants = self.Parse_Variants(self.VCF_path, self.CaseControl, genes, self.cores)

        case_variants.to_csv(self.outpath + "Cases_Intermediate.csv")
        control_variants.to_csv(self.outpath + "Controls_Intermediate.csv")
