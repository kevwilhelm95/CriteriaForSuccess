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
                            "EAML_output/meanMCC-results.nonzero-stats.rankings.csv")
            self.num_genes = eaml.shape[0]
            eaml_fdr01 = eaml[(eaml['qvalue'] <= 0.1) & (eaml['zscore'] >= 0)]
            eaml_fdr001 = eaml[(eaml['qvalue'] <= 0.01) & (eaml['zscore'] >= 0)]

            # EPIMUTESTR
            epi = pd.read_csv(self.path + "EPIMUTESTR_output/EPI_output.tsv",
                            sep='\t', header=None, names=['Genes', 'pvalue', 'qvalue'])
            epi_fdr01 = epi[epi["qvalue"] <= 0.1]
            epi_fdr001 = epi[epi["qvalue"] <= 0.01]

            # EA-Wavelet
            eaw = pd.read_csv(self.path + "EAWavelet_output/wavelet_output.csv")
            eaw_fdr01 = eaw[eaw['fdr'] <= 0.1]
            eaw_fdr001 = eaw[eaw['fdr'] <= 0.01]

            # Reactome and Genes
            reactome = pd.read_csv(self.path + "Pathway_output/AC" + str(self.acThresh) + "/Reactome_Cases/Step8_SigCoreGenesAbove5thPercentile.txt", sep='\t', header=None, skiprows=2, names=['Genes'])
            string = pd.read_csv(self.path + "Pathway_output/AC" + str(self.acThresh) + "/STRING_Cases/Step8_SigCoreGenesAbove5thPercentile.txt", sep="\t", header=None, skiprows=2, names=['Genes'])

            # Create dataframes
            self.fdr01_df = pd.DataFrame(
                {'EAML': eaml_fdr01.gene.dropna(), 'EPI': epi_fdr01.Genes.dropna(), 'EAW': eaw_fdr01.gene.dropna(), 'Reactome': reactome.Genes.dropna(), 'STRING': string.Genes.dropna()})
            self.fdr001_df = pd.DataFrame({'EAML': eaml_fdr001.gene, 'EPI': epi_fdr001.Genes, 'EAW': eaw_fdr001.gene, 'Reactome': reactome.Genes.dropna(), 'STRING': string.Genes.dropna()})


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
                                    "Consensus2": self.consensus2['List5'].dropna().reset_index(drop=True)}
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


            # Unique linked genes between methods (Tier 2 genes)
            def UniqueLinked(self):
                # Generalizes the API Call
                def API_STRING(method, params):
                    string_api_url = "https://version-11-5.string-db.org/api"
                    output_format = 'json'

                    ## Construct URL
                    request_url = "/".join([string_api_url, output_format, method])

                    ## Call STRING
                    results = requests.post(request_url, data=params)

                    ## Read and parse the results
                    matchedTerms_STRING = pd.json_normalize(results.json())
                    #print(matchedTerms_STRING)
                    return(matchedTerms_STRING)

                # Recalculates the combined score only using database, experimental, and coexpression evidence
                def ReScore(df):
                    # https://string-db.org/cgi/help.pl?&subpage=faq%23how-are-the-scores-computed
                    # nscore = gene neighborhood score
                    # fscore = gene fusion score
                    # pscore = phylogenetic profile score
                    # ascore = coexpression score
                    # escore = experimental score
                    # dscore = database score
                    # tscore = textmining score

                    # Set prior probability of random gene having an interaction
                    p = 0.041

                    df['ascore'] = [p if x < p else x for x in df.ascore]
                    df['escore'] = [p if x < p else x for x in df.escore]
                    df['dscore'] = [p if x < p else x for x in df.dscore]

                    # Remove prior probability
                    a_noprior = (df.ascore - p) / (1 - p)
                    e_noprior = (df.escore - p) / (1 - p)
                    d_noprior = (df.dscore - p) / (1 - p)

                    # Combine scores of the channels of interest
                    s_tot_nop = 1 - ((1 - a_noprior) *
                                    (1 - e_noprior) * (1 - d_noprior))

                    # Add prior back
                    s_tot = s_tot_nop + p * (1 - s_tot_nop)

                    df['NewCombinedScore'] = s_tot

                    return df

                # Get gene lists
                consensus2 = list(self.consensus_df.Consensus2.dropna())
                genelist = np.unique(list(self.consensus_df.EAML.dropna()) + list(self.consensus_df.EPI.dropna()) + list(self.consensus_df.EAW.dropna()) + list(self.consensus_df.Reactome.dropna()) + list(self.consensus_df.STRING.dropna()))
                genelist = [x for x in genelist if x not in consensus2] # The unique genes identified that are not in Consensus list

                # Get STRING protein ids for gene list
                params = {
                    "identifiers" : "\r".join(genelist),
                    "species" : 9606,
                    "limit" : 1,
                    "echo_query" : 1
                }
                string_ids = API_STRING('get_string_ids', params)
                
                # Get STRING interaction network
                params = {
                    'identifiers' : '\r'.join(list(string_ids.preferredName)),
                    'species' : 9606,
                    'network_type' : 'functional',
                    'add_nodes' : 0
                }
                interactions = API_STRING('network', params)
                interactions.drop_duplicates(inplace = True)
                
                # Re-score interactions using only Databases, Experiments, Coexpression
                interactions_filt = ReScore(interactions)
                interactions_filt = interactions_filt[interactions_filt.NewCombinedScore >= 0.9].reset_index(drop=True)

                # Map our gene names to STRING preferred Names
                interactions_filt = pd.merge(interactions_filt, string_ids[['preferredName', 'queryItem']], left_on = 'preferredName_A', right_on = 'preferredName', suffixes = (None, 'A'))
                interactions_filt = pd.merge(interactions_filt, string_ids[['preferredName', 'queryItem']], left_on = 'preferredName_B', right_on = 'preferredName', suffixes = (None, 'B'))
                interactions_filt.rename(columns = {'queryItem': 'analysisName_A', 'queryItemB': 'analysisName_B'}, inplace = True)
                interactions_filt.drop(columns = ['preferredName', 'preferredNameB'], inplace = True)

                # Remove interactions where interaction is with a gene from same method
                out_df = pd.DataFrame(
                    columns=['Node1', 'Method1', 'Node2', 'Method2'])
                uniq_linked = []
                intra_linked_index = []

                for row in range(interactions_filt.shape[0]):
                    # Get gene names from STRING API Interaction network
                    node1 = interactions_filt.analysisName_A[row]
                    node2 = interactions_filt.analysisName_B[row]

                    # Find the method that each gene is identified by
                    try:
                        method1 = self.consensus_df[self.consensus_df == node1].stack().index.tolist()[0][1]
                        method2 = self.consensus_df[self.consensus_df == node2].stack().index.tolist()[0][1]
                    except:
                        continue

                    if method1 == method2:
                        intra_linked_index.append(row)

                    else:
                        uniq_linked.append(node1)
                        uniq_linked.append(node2)
                        out_df = out_df.append({'Node1': node1,
                                                'Method1': method1,
                                                'Node2': node2,
                                                'Method2': method2
                                                }, ignore_index=True)

                interactions_filt.drop(index = intra_linked_index, inplace = True)

                # Merge genes + method in which it is found into interactions dataframe
                columns = ['Node', 'Method']
                part1 = out_df[['Node1', 'Method1']]
                part2 = out_df[['Node2', 'Method2']]

                part1.columns, part2.columns = columns, columns

                gene_method = pd.concat([part1, part2], ignore_index = True)
                gene_method.drop_duplicates(inplace = True)

                interactions_filt = pd.merge(
                    interactions_filt, gene_method, left_on='analysisName_A', right_on='Node', suffixes=(None, '_A'))
                interactions_filt = pd.merge(interactions_filt, gene_method,
                                            left_on='analysisName_B', right_on='Node', suffixes=(None, '_B'))
                interactions_filt.drop_duplicates(inplace=True)
                interactions_filt.rename(
                    columns={'Method': 'Method_A', 'Method_B': 'Method_B'}, inplace=True)
                interactions_filt.drop(columns=['Node', 'Node_B'], inplace=True)

                # Calculate how many genes are linked to multiple methods
                gene_connection = pd.DataFrame(columns = ['Gene', '# of Methods Connected to'])
                for a in list(set(uniq_linked)):
                    hold = out_df[(out_df['Node1'] == a) | (out_df['Node2'] == a)]
                    methods = set(list(hold.Method1) + list(hold.Method2))

                    gene_connection = gene_connection.append({'Gene': a,
                                                              '# of Methods Connected to': len(methods) - 1}, ignore_index=True)

                # Merge number of methods for each gene
                interactions_filt = pd.merge(
                    interactions_filt, gene_connection, left_on='analysisName_A', right_on='Gene', suffixes=(None, '_A'))
                interactions_filt = pd.merge(interactions_filt, gene_connection,
                                            left_on='analysisName_B', right_on='Gene', suffixes=(None, '_B'))
                interactions_filt.drop_duplicates(inplace=True)
                interactions_filt.rename(columns={"# of Methods Connected to": "Connections_A",
                                                "# of Methods Connected to_B": "Connections_B"}, inplace=True)
                interactions_filt.drop(columns=['Gene', 'Gene_B'], inplace=True)


                # Prepare files for output
                self.consensus_df['UniqLinked'] = gene_connection.Gene
                self.gene_connection = gene_connection.sort_values(
                    by="# of Methods Connected to", ascending = False)
                self.interaction_df = interactions_filt



            # Loop to create consensus lists for each threshold
            df_dict = {'FDR 0.1': self.fdr01_df, 'FDR 0.01': self.fdr001_df}
            for df_name, df in df_dict.items():
                # Create directories for FDR threshold dfs
                os.makedirs(self.opath+df_name, exist_ok=True)
                os.makedirs(self.opath + df_name + '/STRING', exist_ok=True)

                # Get consensus lists - self.all5, self.consensus4, self.consensus3, self.consensus2
                consensus5_2(self, df)

                # Get overlap stats - self.consensus_df, self.overlap_df
                pairwise_method_overlap(self, df)

                # Get the unique linked genes - self.consensus_df, self.gene_connection, self.interaction_df
                UniqueLinked(self)

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
                    self.gene_connection.to_excel(writer, sheet_name = 'UniqueLinkedCounts', index=False)

                self.interaction_df.to_csv(self.opath + df_name + "/STRING/" + "UniqueLinked_STRING-InteractionFile.csv", index = False)

                # Change names of consensus_df for continuing in script
                if df_name == 'FDR 0.1':
                    self.consensus_fdr01 = self.consensus_df
                elif df_name == 'FDR 0.01':
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

    def ExactTest(self):
        '''
        Required Inputs:
        Self with self.path defined as path to .txt of ExactTest.sh output
        '''
        self.exactTest = pd.read_csv(self.path, sep='\t', header=None,
                                     names=['chrom', 'pos', 'ref', 'alt', 'gene', 'ENSP', 'Consequence',
                                            'HGVSp', 'EA', 'AN_0', 'AN_1', 'Cases', 'Controls', 'AC_1', 'AC_Het_1', 'AC_Hom_1', 'AC_0', 'AC_Het_0', 'AC_Hom_0', 'CC_ALL', 'CC_DOM', 'CC_REC', 'AF', 'AC'], low_memory=False)
        # Clean the file now
        self.exactTest = self.exactTest[self.exactTest.Consequence.str.contains("frameshift_variant|missense_variant|stop_gained|stop_lost", case = False)].reset_index() # Need to include splice sites
        self.exactTest['EA-Clean'] = [x.split(',')[0] for x in self.exactTest['EA']]
        self.exactTest['EA-Clean'] = [-1 if x == '.' else x for x in self.exactTest['EA-Clean']]
        self.exactTest['EA-Clean'] = self.exactTest['EA-Clean'].astype(float)
        
        return self.exactTest

