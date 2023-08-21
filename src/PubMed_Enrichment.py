from Bio import Entrez
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from urllib.error import HTTPError
from http.client import IncompleteRead
import time
import concurrent.futures
from itertools import repeat
from helper_functions import CreateDir, ParseGeneLocationFile

plt.rcParams.update({
    'font.family': 'Avenir',
    'font.size': 14
})

# Define Class
class PubMed_Enrichment():
    def __init__(self, df, df_name, interest_list, ref, keywords, output_path):
        self.email = "kwilhelm95@gmail.com"
        self.api_key = '3a82b96dc21a79d573de046812f2e1187508'
        self.df = df
        self.df_name = df_name
        self.interest_list = interest_list
        try: self.key_word_list = keywords.split(",")
        except AttributeError:
            print("Please define '--PubMedKeywords'")
        self.ref = ref
        self.output_path = output_path
        self.main()

    def LoadBackgroundGenes(self):
        ref_df = ParseGeneLocationFile(self.ref)
        background_genes = ref_df.gene
        return background_genes

    # API Call to PubMed
    def search(self, gene, disease, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        new_query = f'"{gene}" AND ("{disease}") AND ("gene" or "protein")'
        try: handle = Entrez.esearch(db = 'pubmed',
                                sort = 'relevance',
                                retmax = '100000',
                                retmode = 'xml',
                                term = new_query)
        except IndexError:
            print(f"{gene} - Not Found")
        except HTTPError:
            print('....Network Error-Waiting 10s')
            time.sleep(10)
            handle = Entrez.esearch(db = 'pubmed',
                                sort = 'relevance',
                                retmax = '100000',
                                retmode = 'xml',
                                term = new_query)
        except IncompleteRead:
            print('....Network Error-Waiting 10s')
            time.sleep(10)
            try: handle = Entrez.esearch(db = 'pubmed',
                                sort = 'relevance',
                                retmax = '100000',
                                retmode = 'xml',
                                term = new_query)
            except IncompleteRead:
                return "IncompleteReadError"

        results = Entrez.read(handle)
        handle.close()
        return results

    def GetEnrichment(self, query, disease_query, background_genes, max_workers, outpath):
        ## Perform PubMed query on our query genes
        print(f'Pulling Publications for {len(query)} Query Genes + {disease_query}')
        df = pd.DataFrame(columns=['Gene and paper', 'Gene and disease'], index=query)
        out_df = pd.DataFrame(columns=['Count', 'PMID for Gene + ' + str(disease_query)], index=query)

        # Thread api calls to Pubmed
        with concurrent.futures.ThreadPoolExecutor(max_workers = max_workers) as executor:
            # Using map, you do not need list iteration
            results = executor.map(self.search, query, repeat(disease_query), repeat(self.email), repeat(self.api_key))

            for result in results:
                try: 
                    gene = result['TranslationStack'][0]['Term'].split('"')[1]
                    n_paper_dis = int(result['Count'])
                except:
                    try: 
                        gene = result['WarningList']['QuotedPhraseNotFound'][0].split('"')[1]
                        n_paper_dis = np.nan
                    except: 
                        gene = result['QueryTranslation'].split('"')[1]
                        n_paper_dis = int(result['Count'])
                df.loc[gene, 'Gene and disease'] = n_paper_dis

                # Append values to out_df
                out_df.loc[gene, 'PMID for Gene + ' + disease_query] = "; ".join(result['IdList'])
                out_df.loc[gene, 'Count'] = n_paper_dis

        # Write summary of query genes to table
        out_df.sort_values(by = 'Count', ascending = False).to_csv(outpath + "/PMIDresults_query+" + disease_query + ".csv")

        # Perform pubmed query on random genes
        trials = 10
        randfs = []
        print('Pulling Publications for {} random gene sets of {} genes'.format(trials, len(query)))
        for i in range(trials):
            if i % 10 == 0:
                print(" Random Trial : ", i)
            randgenes = np.random.choice(background_genes, len(query), replace = False)
            tempdf = pd.DataFrame(columns = ['Gene and paper', 'Gene and disease'])

            with concurrent.futures.ThreadPoolExecutor(max_workers = max_workers) as executor:
                results = executor.map(self.search, randgenes, repeat(disease_query), repeat(self.email), repeat(self.api_key))

                for result in results:
                    try: gene = result['TranslationStack'][0]['Term'].split('"')[1]
                    except: 
                        gene = result['QueryTranslation'].split('"')[1]
                    n_paper_dis = result['Count']
                    tempdf.loc[gene, 'Gene and disease'] = int(n_paper_dis)

            randfs.append(tempdf)

        # Calculate Z-score for observation relative to background
        thrshlds = [[-1, 0], [0, 5], [5, 15], [15, 50], [50, 10000]]
        for paper_thrshld in thrshlds:
            observation = df[(df['Gene and disease'] > paper_thrshld[0]) & (df['Gene and disease'] <= paper_thrshld[1])].shape[0]
            background = [tmp[(tmp['Gene and disease'] > paper_thrshld[0]) & (tmp['Gene and disease'] <= paper_thrshld[1])].shape[0] for tmp in randfs]
            savefile = f"Enrichment_query+{disease_query}_>{paper_thrshld[0]}-<={paper_thrshld[1]}.png"
            back_mean = np.mean(background)
            back_std = np.std(background)
            z = (observation - back_mean)/back_std

            # Print results
            try:
                if -1 in paper_thrshld:
                    xlabel = '# of Genes with {} Co-Mentions with "{}"'.format(paper_thrshld[0]+1, disease_query)
                    observation = "{}/{} genes had 0 co-mentions with {} -- Z = {}".format(observation, len(query), disease_query, z)
                else: 
                    xlabel = '# of Genes with {}-{} Co-Mentions with "{}"'.format(paper_thrshld[0]+1, paper_thrshld[1], disease_query)
                    observation = "{}/{} genes had >{} & <={} co-mentions with {} -- Z = {}".format(observation, len(query), paper_thrshld[0], paper_thrshld[1], disease_query, z)

                with open(outpath + f"{disease_query}_Results.txt", 'rw') as f:
                    f.write(observation + "\n")

                # Plot Observation and Random Tests
                fig, ax = plt.subplots(
                    figsize=(6, 3.5), facecolor='white', frameon=False)

                y, x, _ = plt.hist(background, color='dimgrey', edgecolor='white')
                plt.axvline(x=observation, ymax=0.5, linestyle='dotted', color='red')
                plt.text(x=observation*0.99, y=(y.max()/1.95), s='{}/{} (Z = {})'.format(observation, len(query), round(z, 2)), color='red', ha='right')
                plt.xlabel(xlabel, fontsize=15)
                plt.ylabel('Count', fontsize=15)
                plt.savefig(outpath + savefile)
                plt.close()
            except:
                print("PubMed Enrichment - Error Plotting")
                continue

        
    def main(self):
        # Load background genes
        background_genes = self.LoadBackgroundGenes()

        # Loop through gene lists
        for gL in self.interest_list:
            gL_hold = self.df[gL].dropna()
            # Loop through keyword queries
            for keyword in self.key_word_list:
                print(f"--- Querying {gL} genes for co-mentions with {keyword}---")
                # Create directory
                hold_OutPutPath = CreateDir(self.output_path, f"{gL}/{keyword}")
                self.GetEnrichment(gL_hold, keyword, background_genes, max_workers = 6, outpath = hold_OutPutPath)
