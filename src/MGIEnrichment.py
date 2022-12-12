import pandas as pd
from scipy.stats import hypergeom, fisher_exact
import numpy as np
import os
from statsmodels.stats.multitest import fdrcorrection as fdr
from helper_functions import CreateDir

class MGIEnrichment():
    def __init__(self, mgi, df, df_name, experiment_name, output_path, interst_list, cores):
        self.MGI = mgi
        self.consDf = df
        self.consDfName = df_name
        self.expName = experiment_name
        self.oPath = output_path
        self.interstList = interst_list
        self.Cores = cores
        
        self.main()



    def main(self):
        # Clean MGI file and create pheno_dict for mapping - self.mgi_genes, self.interst, self.pheno_dict
        def CleanMGI(self):
            self.mgi_genes = list(self.MGI[0])

            # Get unique high-level phenotype codes for MGI
            pheno = self.MGI[4].dropna()
            unique_pheno = []
            for i in pheno:
                hold_pheno = list(i.split(', '))
                for x in hold_pheno:
                    if x not in unique_pheno:
                        unique_pheno.append(x)
            self.interst = unique_pheno

            # Define the unique phenotypes and their names
            self.pheno_dict = {"MP:0005387": "Immune System",
                               "MP:0003631": "Nervous System",
                               "MP:0005369": "Muscle",
                               "MP:0010768": "Mortality/Aging",
                               "MP:0005386": "Behavior/Neurological",
                               "MP:0005375": "Adipose Tissue",
                               "MP:0005385": "Cardiovascular",
                               "MP:0005381": "Digestive/Alimentary",
                               "MP:0005379": "Endocrine/Exocrine",
                               "MP:0005378": "Growth/Size/Body",
                               "MP:0005376": "Homeostasis/Metabolism",
                               "MP:0005384": "Cellular",
                               "MP:0005382": "Craniofacial",
                               "MP:0005380": "Embryo",
                               "MP:0005377": "Hearing/Vestibular System",
                               "MP:0005397": "Hematopoietic System",
                               "MP:0010771": "Integument",
                               "MP:0005371": "Limbs/Digits/Tail",
                               "MP:0005370": "Liver/Biliary System",
                               "MP:0002006": "Neoplasm",
                               "MP:0001186": "Pigmentation",
                               "MP:0005367": "Renal/Urinary System",
                               "MP:0005389": "Reproductive System",
                               "MP:0005388": "Respiratory System",
                               "MP:0005390": "Skeleton",
                               "MP:0005394": "Taste/Olfaction",
                               "MP:0005391": "Vision/Eye"}
        
        # Calculate p-value for enrichment against MGI
        def CalculateP(self, phen_interst, mgi, all_repeated):
            phen = 0
            for i, r in mgi.iterrows():
                phenlist = str(mgi.loc[i, 4])
                phenlist = list(phenlist.split(', '))
                # and phen[1] in phenlist: ##### phenotype conditionals
                if phen_interst in phenlist:
                    phen += 1

            query = all_repeated
            query = list(set(query).intersection(set(self.mgi_genes)))
            query = mgi[mgi[0].isin(query)].copy()
            num_genes_in_mgi = query.shape[0]

            query_phen = 0
            cand_list = []
            for i, r in query.iterrows():
                phenquery = str(query.loc[i, 4])
                phenquery = list(phenquery.split(', '))
                # and phen[1] in phenquery: ##### phenotype conditionals
                if phen_interst in phenquery:
                    query_phen += 1
                    cand_hold = str(query.loc[i, 0])
                    cand_list.append(cand_hold)

            pval = hypergeom.sf(query_phen-1,
                        mgi.shape[0], query.shape[0], phen)
            # Fishers Exact Test
            #table = |TP  FP|
            #        |FN  TN|
            tp = query_phen
            fp = query.shape[0] - query_phen
            fn = phen - query_phen
            tn = mgi.shape[0] - phen + query_phen - query.shape[0]

            table = np.array([[tp, fp], [fn, tn]])

            orr, f_pval = fisher_exact(table, alternative='greater')
            
            return_df = pd.DataFrame( {"Phenotype": [idx],
                                    "Phenotype Name": [self.pheno_dict[idx]],
                                    "Total Genes": [len(all_repeated)],
                                    "MGI Genes with Phenotype": [phen],
                                    "Candidate Genes in MGIdb": [num_genes_in_mgi],
                                    "Enrichment - p value": [f_pval],
                                    "# of Candidates with Phenotypic Model": [query_phen],
                                    "Ratio of genes with phenotypic model": [query_phen / len(all_repeated)],
                                    "Genes with MGI Phenotype": [str(cand_list)]} )

            return return_df


        CleanMGI(self)
        # Loop through each method gene list and calculate the enrichment for a phenotype in MGI
        for gL in self.interstList:
            # Create new output location for each method's enrichment
            newMGI_outpath = CreateDir(self.oPath, f'{gL}/')
            
            # Define output df
            out_df = pd.DataFrame(columns=['Phenotype', 'Phenotype Name', 'Total Genes', 'MGI genes with Phenotype',
                                           'Candidate Genes in MGIdb', '# of Candidates with Phenotypic Model',
                                           'Enrichment - p value', 'Ratio of genes with phenotypic model',
                                           'Genes with MGI Phenotype'])
            for idx in self.interst:
                # Subset gL
                gL_hold = self.consDf[gL].dropna()

                p_df = CalculateP(self, idx, self.MGI, gL_hold)
                out_df = pd.concat([out_df, p_df], axis=0, ignore_index=True)
            
            # Calculate fdr for out_df
            out_df = out_df.sort_values("Enrichment - p value", ascending = True).reset_index(drop=True)
            out_df['FDR'] = fdr(list(out_df['Enrichment - p value']), is_sorted = True)[1]

            # Export file to .csv
            out_df.to_csv(newMGI_outpath + gL + "_MGI Enrichment.csv", index=False)
