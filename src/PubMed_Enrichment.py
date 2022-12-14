from Bio import Entrez
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import concurrent.futures
from itertools import repeat
from helper_functions import CreateDir

# Define Class
class PubMed_Enrichment():
    def __init__(self, df, df_name, interest_list, keywords, output_path):
        self.email = "kwilhelm95@gmail.com"
        self.api_key = '3a82b96dc21a79d573de046812f2e1187508'
        self.df = df
        self.df_name = df_name
        self.interest_list = interest_list
        self.key_word_list = keywords.split(",")
        self.background_path = f"{os.path.abspath(os.path.join(os.getcwd(), '../'))}/refs/background_genes_n=19819.csv"
        self.output_path = output_path
        self.main()

    def LoadBackgroundGenes(self):
        self.background_genes = pd.read_csv(self.background_path, sep ='\t', )

    # API Call to PubMed
    def search(gene, disease, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        new_query = f'"{gene}" AND ("{disease}") AND ("gene" or "protein")'
        handle = Entrez.esearch(db = 'pubmed',
                                sort = 'relevance',
                                retmax = '100000',
                                retmode = 'xml',
                                term = new_query,
                                field = 'title/abstract')
        results = Entrez.read(handle)
        return results

    def GetEnrichment(self, query, disease_query, max_workers, outpath):
        ## Perform PubMed query on our query genes
        print('Pulling Publications for Query Genes')
        # Create hold df
        df = pd.DataFrame(columns=['Gene and paper', 'Gene and disease'], index=query)
        # Create the output df
        out_df = pd.DataFrame(columns=['Count', 'PMID for Gene + ' + str(disease_query)], index=query)

        # Thread api calls to Pubmed
        with concurrent.futures.ThreadPoolExecutor(max_workers = max_workers) as executor:
            # Using map, you do not need list iteration
            results = executor.map(self.search, query, repeat(disease_query), repeat(self.email), repeat(self.api_key))

            for result in results:
                try: gene = result['TranslationStack'][0]['Term'].split('"')[1]
                except:
                    gene = result['QueryTranslation'].split('"')[1]
                n_paper_dis = result['Count']
                df.loc[gene, 'Gene and disease'] = int(n_paper_dis)

                # Append values to out_df
                out_df.loc[gene, 'PMID for Gene + ' + disease_query] = "; ".join(result['IdList'])
                out_df.loc[gene, 'Count'] = int(n_paper_dis)

        # Write summary of query genes to table
        out_df.sort_values(by = 'Count', ascending = False).to_csv(outpath + "PMIDresults_Title-Abstract_query+" + disease_query + ".csv")

        # Perform pubmed query on random genes
        trials = 100
        randfs = []
        print('Pulling Publications for {} random gene sets of {} genes'.format(trials, len(query)))
        for i in range(trials):
            if i % 10 == 0:
                print(" Random Trial : ", i)
            randgenes = np.random.choice(self.background_genes, len(query), replace = False)
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
        thrshlds = [0, 4, 9] # The counts are offset by one - if you are looking for 5 or more papers, put 4
        for paper_thrshld in thrshlds:
            observation = df[df['Gene and disease'] > paper_thrshld].shape[0]
            background = [tmp[tmp['Gene and disease'] > paper_thrshld].shape[0] for tmp in randfs]
            back_mean = np.mean(background)
            back_std = np.std(background)
            z = (observation - back_mean)/back_std

            # Print results
            print('Observation is that {} out of {} genes had at least {} publications with "{}"'.format(observation, len(query), (paper_thrshld+1), disease_query))
            print('Z-score of this was {}'.format(z))

            # Plot Observation and Random Tests
            fig, ax = plt.subplots(figsize = (6,3.5), facecolor = 'white', frameon = False)

            y, x, _ = plt.hist(background, color = 'dimgrey', edgecolor = 'white')
            plt.axvline(x = observation, ymax = 0.5, linestyle='dotted', color = 'red')
            plt.text(x = observation*0.99, y = (y.max()/1.95), s = '{}/{} (Z = {})'.format(observation, len(query), round(z, 2)), color = 'red', ha = 'right')
            plt.xlabel(f'# of genes with {paper_thrshld+1} papers co-mentioning "{disease_query}"'.format(disease_query), fontsize = 15)
            plt.ylabel('Count', fontsize=15)
            plt.show()
            plt.close()
            plt.savefig(outpath + f"{disease_query}_Plot.png")
        
    def main(self):
        # Load background genes
        self.LoadBackgroundGenes

        # Loop through gene lists
        for gL in self.interest_list:
            gL_hold = gL.dropna()
            # Loop through keyword queries
            for keyword in self.key_word_list:
                # Create directory
                hold_OutPutPath = CreateDir(self.output_path, f"{gL}/{keyword}")
                self.GetEnrichment(gL_hold, keyword, max_workers = 16, outpath = hold_OutPutPath)