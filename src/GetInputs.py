'''
author: @kevwilhelm95
Built from version 2 -
Purpose - Load and clean input files needed for Criteria For Success
'''
import pandas as pd
from scipy.stats import hypergeom as hg
import numpy as np
import os
import upsetplot as up
import matplotlib.pyplot as plt
import requests


class GetInputs():
    def __init__(self, input_path, ac_thresh, output_path, exp_name):
        '''
        Takes as input the argument path to a directory containing input data
        '''
        self.path = input_path
        self.acThresh = ac_thresh
        self.opath = output_path
        self.expName = exp_name

    def BigPipeline(self, analysis):
        ## Define functions needed for BigPipeline
        # Loads and merges gene lists from machine learning methods
        def LoadLists(self):
            # EAML
            eaml = pd.read_csv(self.path +
                               "EAML_output/meanMCC-results.nonzero-stats.rankings", sep = ',')
            self.num_genes = eaml.shape[0]
            eaml_fdr01 = eaml[(eaml['qvalue'] <= 0.1)
            eaml_fdr001 = eaml[(eaml['qvalue'] <= 0.01)]

            # EPIMUTESTR
            epi = pd.read_csv(self.path + "EPIMUTESTR_output/EPI_output.tsv",
                            sep='\t', header=None, names=['Genes', 'pvalue', 'qvalue'])
            epi_fdr01 = epi[epi["qvalue"] <= 0.1]
            epi_fdr001 = epi[epi["qvalue"] <= 0.01]

            # EA-Wavelet
            eaw = pd.read_csv(self.path + "EAWavelet_output/wavelet_output.csv")
            eaw_fdr01 = eaw[eaw['fdr'] <= 0.1]
            eaw_fdr001 = eaw[eaw['fdr'] <= 0.01]

            # Reactome and STRING Genes
            reactome = pd.read_csv(self.path + "Pathway_output/AC" + str(self.acThresh) + "/Reactome_Cases/Step8_SigCoreGenesAbove5thPercentile.txt", sep='\t', header=None, skiprows=2, names=['Genes'])
            string = pd.read_csv(self.path + "Pathway_output/AC" + str(self.acThresh) + "/STRING_Cases/Step8_SigCoreGenesAbove5thPercentile.txt", sep="\t", header=None, skiprows=2, names=['Genes'])

            # Create dataframes
            self.fdr01_df = pd.DataFrame({'EAML': eaml_fdr01.gene.dropna(),
                                        'EPI': epi_fdr01.Genes.dropna(), 
                                        'EAW': eaw_fdr01.gene.dropna(), 
                                        'Reactome': reactome.Genes.dropna(), 
                                        'STRING': string.Genes.dropna()})
            self.fdr001_df = pd.DataFrame({'EAML': eaml_fdr001.gene, 
                                            'EPI': epi_fdr001.Genes, 
                                            'EAW': eaw_fdr001.gene, 
                                            'Reactome': reactome.Genes.dropna(), 
                                            'STRING': string.Genes.dropna()})


        # Creates three sheets - Meta sheet + Consensus lists, Stats for pairwise method overlap, Count for # of methods that find a certain gene
        def consensus(self):
            # Find overlapping genes
            def consensus5_2(self, df):
                ### Getting consensus lists
                list_5 = list(df.EPI.dropna()) + list(df.EAML.dropna()) + list(df.EAW.dropna()) + list(df.Reactome.dropna()) + list(df.STRING.dropna())
                list_5_df = pd.DataFrame({'List5': list_5})
                list_5_group = pd.DataFrame(list_5_df.groupby("List5", as_index=False)[
                                            'List5'].agg({'Count': 'count'}))
                list_5_group = list_5_group.rename({"List5": "Count"})
                # All5
                self.all5 = list_5_group[list_5_group.Count == 5]
                self.allunique = np.unique(list_5)
                # Consensus4
                self.consensus4 = list_5_group[list_5_group.Count >= 4]
                # Consensus3
                self.consensus3 = list_5_group[list_5_group.Count >= 3]
                #Consensus2
                self.consensus2 = list_5_group[list_5_group.Count >= 2]


            # Get pairwise overlap between gene lists
            def pairwise_method_overlap(self, df):
                # From EPIMUTESTR
                epi_eaw = df['EPI'][df['EPI'].isin(
                    df['EAW'])].dropna().reset_index(drop=True)
                epi_eaw_p = hg(M=self.num_genes, n=len(df['EPI'].dropna()), N=len(
                    df['EAW'].dropna())).sf(len(epi_eaw)-1)

                epi_eaml = df['EPI'][df['EPI'].isin(
                    df['EAML'])].dropna().reset_index(drop=True)
                epi_eaml_p = hg(M=self.num_genes, n=len(df['EPI'].dropna()), N=len(
                    df['EAW'].dropna())).sf(len(epi_eaml)-1)

                epi_reactome = df['EPI'][df['EPI'].isin(
                    df['Reactome'])].dropna().reset_index(drop=True)
                epi_reactome_p = hg(M=self.num_genes, n=len(df['EPI'].dropna()), N=len(
                    df['Reactome'].dropna())).sf(len(epi_reactome)-1)

                epi_string = df['EPI'][df['EPI'].isin(
                    df['STRING'])].dropna().reset_index(drop=True)
                epi_string_p = hg(M=self.num_genes, n=len(df['EPI'].dropna()), N=len(
                    df['STRING'].dropna())).sf(len(epi_string)-1)

                # From EAW
                eaw_eaml = df['EAW'][df['EAW'].isin(
                    df['EAML'])].dropna().reset_index(drop=True)
                eaw_eaml_p = hg(M=self.num_genes, n=len(df['EAW'].dropna()), N=len(
                    df['EAML'].dropna())).sf(len(eaw_eaml)-1)

                eaw_reactome = df['EAW'][df['EAW'].isin(
                    df['Reactome'])].dropna().reset_index(drop=True)
                eaw_reactome_p = hg(M=self.num_genes, n=len(df['EAW'].dropna()), N=len(
                    df['Reactome'].dropna())).sf(len(eaw_reactome)-1)

                eaw_string = df['EAW'][df['EAW'].isin(
                    df['STRING'])].dropna().reset_index(drop=True)
                eaw_string_p = hg(M=self.num_genes, n=len(df['EAW'].dropna()), N=len(
                    df['STRING'].dropna())).sf(len(eaw_string)-1)

                # From EAML
                eaml_reactome = df['EAML'][df['EAML'].isin(
                    df['Reactome'])].dropna().reset_index(drop=True)
                eaml_reactome_p = hg(M=self.num_genes, n=len(df['EAML'].dropna()), N=len(
                    df['Reactome'].dropna())).sf(len(eaml_reactome)-1)

                eaml_string = df['EAML'][df['EAML'].isin(
                    df['STRING'])].dropna().reset_index(drop=True)
                eaml_string_p = hg(M=self.num_genes, n=len(df['EAML'].dropna()), N=len(
                    df['STRING'].dropna())).sf(len(eaml_string)-1)

                # From Reactome
                reactome_string = df['Reactome'][df['Reactome'].isin(
                    df['STRING'])].dropna().reset_index(drop=True)
                reactome_string_p = hg(M=self.num_genes, n=len(df['Reactome'].dropna()), N=len(
                    df['STRING'].dropna())).sf(len(reactome_string)-1)

                # Compile into a dataframe with each experimental list, pairwise overlap
                consensus_out_dict = {"EPI": df.EPI.dropna(),
                                    "EAW": df.EAW.dropna(),
                                    "EAML": df.EAML.dropna(),
                                    "Reactome": df.Reactome.dropna(),
                                    "STRING": df.STRING.dropna(),
                                    "EPI-EAW": epi_eaw,
                                    "EPI_EAML": epi_eaml,
                                    "EPI_Reactome": epi_reactome,
                                    "EPI_STRING": epi_string,
                                    "EAW_EAML": eaw_eaml,
                                    "EAW_Reactome": eaw_reactome,
                                    "EAW_STRING": eaw_string,
                                    "EAML_Reactome": eaml_reactome,
                                    "EAML_STRING": eaml_string,
                                    "Reactome_STRING": reactome_string,
                                    "All5": self.all5['List5'].dropna().reset_index(drop=True),
                                    "Consensus4": self.consensus4['List5'].dropna().reset_index(drop=True),
                                    "Consensus3": self.consensus3['List5'].dropna().reset_index(drop=True),
                                    "Consensus2": self.consensus2['List5'].dropna().reset_index(drop=True),
                                    "AllUnique": pd.Series(self.allunique)}
                self.consensus_df = pd.DataFrame(consensus_out_dict)

                overlap_out_dict = {"List1": ['EPI', 'EPI', 'EPI', 'EPI', 'EAW', 'EAW', 'EAW',
                                              'EAML', 'EAML', 'Reactome'],
                                    "List2": ['EAW', 'EAML', 'Reactome', 'STRING', 'EAML', 'Reactome',
                                              'STRING', 'Reactome', 'STRING', 'STRING'],
                                    "Overlap": [len(epi_eaw), len(epi_eaml), len(epi_reactome), len(epi_string), len(eaw_eaml), len(eaw_reactome), len(eaw_string), len(eaml_reactome), len(eaml_string), len(reactome_string)],
                                    "Pval": [epi_eaw_p, epi_eaml_p, epi_reactome_p, epi_string_p, eaw_eaml_p, eaw_reactome_p, eaw_string_p, eaml_reactome_p, eaml_string_p, reactome_string_p]}
                self.overlap_df = pd.DataFrame(overlap_out_dict)


            # Count how many methods find a particular gene
            def geneMethodCount(self):
                # Create holder
                method_count_df = pd.DataFrame(columns=['Genes', '# of Methods', 'Methods'])

                # Loop through Consensus2 and check if gene in each gene list
                for gene in self.consensus_df.Consensus2.dropna():
                    count = 0
                    method_hold = []
                    
                    for method in ['EPI', 'EAML', 'EAW', 'Reactome', 'STRING']:
                        if gene in list(self.consensus_df[method].dropna()):
                            count += 1
                            method_hold.append(method)

                    #method_count_df = method_count_df.append({'Genes': gene, '# of Methods': count, 'Methods': str(method_hold)}, ignore_index=True)
                    hold_df = pd.DataFrame.from_dict({'Genes':[gene], '# of Methods':[count], 'Methods':[str(method_hold)]})
                    method_count_df = pd.concat([method_count_df, hold_df], axis = 0, ignore_index=True)

                self.method_count_df = method_count_df.sort_values(
                    '# of Methods', ascending=False)


            # UpsetPlot 
            def ConsensusUpSet(self, df_name):
                contents = {'EPI': self.consensus_df.EPI.dropna(),
                            'EAW': self.consensus_df.EAW.dropna(),
                            'EAML': self.consensus_df.EAML.dropna(),
                            'Reactome': self.consensus_df.Reactome.dropna(),
                            'STRING': self.consensus_df.STRING.dropna()}
                from_contents = up.from_contents(contents)
                up.plot(from_contents, sort_by='degree', show_counts = True)
                plt.savefig(self.opath + df_name + "/UpSetPlot.png")
                plt.close()

            # Loop to create consensus lists for each threshold
            df_dict = {'FDR_0.1': self.fdr01_df, 'FDR_0.01': self.fdr001_df}
            for df_name, df in df_dict.items():
                # Create directories for FDR threshold dfs
                os.makedirs(self.opath+df_name, exist_ok=True)
                os.makedirs(self.opath + df_name + '/STRING', exist_ok=True)

                # Get consensus lists - self.all5, self.consensus4, self.consensus3, self.consensus2
                consensus5_2(self, df)

                # Get overlap stats - self.consensus_df, self.overlap_df
                pairwise_method_overlap(self, df)

                # Get the number of methods that identify each gene in Consensus2
                geneMethodCount(self)

                # Create and write upset plot
                ConsensusUpSet(self, df_name)

                # Write Excel file for consensus_df, overlap_df, method_count_df
                with pd.ExcelWriter(self.opath + df_name + "/" + self.expName + "_Method+ConsensusFiles.xlsx") as writer:
                    self.consensus_df.to_excel(
                        writer, sheet_name="OutputLists+Consensus", index=False)
                    self.overlap_df.to_excel(
                        writer, sheet_name="Overlap Summary", index=False)
                    self.method_count_df.to_excel(
                        writer, sheet_name="Gene-Method Counts", index=False)

                # Change names of consensus_df for continuing in script
                if df_name == 'FDR_0.1':
                    self.consensus_fdr01 = self.consensus_df
                elif df_name == 'FDR_0.01':
                    self.consensus_fdr001 = self.consensus_df

        def LoadInputList(self):
            genes = pd.read_csv(self.path, sep='\t', header = None, names = ['Genes'])

            return genes

        # Function call
        if analysis == 'BigPipeline':
            LoadLists(self)
            consensus(self)
            
            return self.consensus_fdr01, self.consensus_fdr001, self.num_genes


        if analysis == 'InputList':
            genes = LoadInputList(self)

            return genes


        # Return consensus_dfs and the background # of genes to CriteriaForSuccess.py
        return self.consensus_fdr01, self.consensus_fdr001, self.num_genes

    def GoldStandards(self):
        ''' 
        Required Inputs:
        Self with self.path defined as the path to the .csv of the GoldStandard Lists
        '''

        self.goldstandards = pd.read_csv(self.path, sep=',')
        return self.goldstandards

    def MGI(self):
        '''
        Required Inputs:
        None
        - Path to Enrichment file is stationary in CriteriaForSuccess/refs
        '''
        self.mgi = pd.read_csv("../refs/HMD_HumanPhenotype_v2021.rpt", delimiter = "\t", header = None)
        return self.mgi

    def CaseControl(self):
        '''
        Required Inputs:
        Self with self.path defined as path to .csv of CaseControl file
        - 2 column file (1-Sample IDs, 2-Case/Control(1/0))
        '''
        self.caseControl = pd.read_csv(self.path, sep=',', header=None)
        return self.caseControl

