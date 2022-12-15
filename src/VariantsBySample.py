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
    def parse_samples(self, samples):
        cases = samples[samples.iloc[:,0] == 1].index.astype(str).tolist()
        conts = samples[samples.iloc[:,0] == 0].index.astype(str).tolist()
        return cases,conts

    # Convert a genotype tuple to a zygosity integer
    def convert_zygo(self, genotype):
        """
        Convert a genotype tuple to a zygosity integer
        Args:
            genotype (tuple): The genotype of a variant for a sample
        Returns:
            int: The zygosity of the variant (0/1/2)
        """
        if genotype in [(1, 0), (0, 1)]:
            zygo = 1
        elif genotype == (1, 1):
            zygo = 2
        else:
            zygo = 0
        return zygo


    # Validate the EA score pulled from VCF
    def validate_EA(self, ea):
        """
        Checks for valid EA score
        Args:
            ea (str/float/None): EA score as string
        Returns:
            float: EA score between 0-100 if valid, otherwise returns NaN
        """
        try:
            ea = float(ea)
        except ValueError:
            if type(ea) == str and (ea == 'fs-indel' or 'STOP' in ea):
                ea = 100
            else:
                ea = np.nan
        except TypeError:
            ea = np.nan
        return ea

    # Fetch EA score of canonical ENSP id
    def fetch_EA_VEP(self, EA, canon_ensp, all_ensp, csq, EA_parser = 'canonical'):
        """
        Fetch the EA score index based on the canonical ENSP id index
        Args:
            EA (_type_): _description_
            canon_ensp (_type_): _description_
            all_ensp (_type_): _description_
            csq (_type_): _description_
            EA_parser (str, optional): _description_. Defaults to 'canonical'.
        """
        if 'stop_gained' in csq or 'frameshift_variant' in csq or 'stop_lost' in csq or 'splice_donor_variant' in csq or 'splice_acceptor_variant' in csq:
            return 100
        if EA_parser == 'canonical':
            try:
                canon_idx = all_ensp.index(canon_ensp)
            except ValueError:
                return np.nan
            else:
                return self.validate_EA(EA[canon_idx])
        else:
            newEA = []
            for score in EA:
                newEA.append(self.validate_EA(score))
            if np.isnan(newEA).all():
                return np.nan
            elif EA_parser == 'mean':
                return np.nanmean(newEA)
            elif EA_parser == 'max':
                return np.nanmax(newEA)
            else:
                return newEA

    # Read through vcf and fetch columns of interest
    def parse_VEP(self, vcf_fn, genes, sample, min_af = None, max_af = None, af_field = 'AF', EA_parser = 'canonical'):
        row=[]   
        sample=[sample]
        vcf = VariantFile(vcf_fn)
        vcf.subset_samples(sample)
            
        for rec in vcf.fetch():
            zyg = self.convert_zygo(rec.samples[sample[0]]['GT'])
            if (zyg!=0) and (rec.info['SYMBOL'][0] in list(genes)):
                all_ea = rec.info.get('EA', (None,))
                all_ensp = rec.info.get('Ensembl_proteinid', (rec.info['ENSP'][0],))
                ea = self.fetch_EA_VEP(all_ea, rec.info['ENSP'][0], all_ensp, rec.info['Consequence'][0], EA_parser='canonical')
                if not np.isnan(ea):
                    row.append([rec.chrom, rec.pos, rec.ref, rec.alts[0], ea, rec.info['SYMBOL'][0], sample[0],zyg, rec.info['AF'][0] ])
        cols = ['chr','pos','ref','alt','EA','gene','sample','zyg','AF']
        col_type = {'chr': str, 'pos': str, 'ref': str, 'alt': str}
        df = pd.DataFrame(row, columns = cols)
        df = df.astype(col_type)
        return df
    

    def Parse_Variants(self, args_vcf, args_samples, genes, args_cores):
        cases, conts = self.parse_samples(args_samples)
        print(cases)
        print(conts)
        conts_dfs = Parallel(n_jobs=int(args_cores))(delayed(self.parse_VEP)(args_vcf, genes, sample, min_af=None, max_af=None, af_field='AF', EA_parser='canonical') for sample in tqdm(conts))
        conts_var = pd.concat(conts_dfs, axis=0)
        conts_var.reset_index(drop = True, inplace = True)
                            
        cases_dfs = Parallel(n_jobs=int(args_cores))(delayed(self.parse_VEP)(args_vcf, genes, sample, min_af=None, max_af=None, af_field='AF', EA_parser='canonical') for sample in tqdm(cases))
        cases_var = pd.concat(cases_dfs, axis=0)
        cases_var.reset_index(drop = True, inplace = True)

        return cases, conts, cases_var, conts_var   


    def main(self):
        # Get genes
        genes = self.GetGenes()

        # Parse variants for outputting
        cases, controls, case_variants, control_variants = self.Parse_Variants(self.VCF_path, self.CaseControl, genes, self.cores)
        all_variants = pd.concat([case_variants, control_variants], axis = 0, ignore_index = True)

        case_variants.to_csv(self.outpath + "Cases_VariantsBySample.csv", sep = ',', index = False)
        control_variants.to_csv(self.outpath + "Controls_VariantsBySample.csv", sep = ',', index = False)
        all_variants.to_csv(self.outpath + "CaseControl_VariantsBySample.csv", sep = ',', index = False)
