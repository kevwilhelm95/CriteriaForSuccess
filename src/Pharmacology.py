import pandas as pd
import os
import requests
import io
import datetime
import multiprocessing as mp

# Create GetGeneDrugInteractions class
class GetGeneDrugInteractions():
    # Define init
    def __init__(self, df, df_name, interst_list, experiment_name, outpath, cores):
        self.consDf = df
        self.consDfName = df_name
        self.methodList = interst_list
        self.expName = experiment_name
        self.oPath = outpath
        self.Cores = cores

        self.main()

    # Get uniprot IDs
    def GetUniprotIDs(self, col):
        '''
        Input: 1) Self, 2) a column from a dataframe (df.EAML)
        '''
        # Parse gene names from col into list
        genes = [x for x in col]

        # Set holding directory
        accessionDict = {}

        # Loop through genes and request accession number from UniProt
        for gene in genes:
            response = requests.get("https://www.uniprot.org/uniprot/?query=gene_exact:" + gene +
                                    "+AND+organism:9606&columns=id,genes(PREFERRED),genes(ALTERNATIVE)&format=tab")
            if response.status_code != 200:
                print("Error in UniProt API call for " + gene)
                continue

            # Parse into pandas dataframe
            response_io = io.StringIO(str((response.text)))
            try:
                response_df = pd.read_table(
                    response_io, sep='\t', dtype='string')
            except pd.errors.EmptyDataError:
                print("No data for " + gene)
                accessionDict[gene] = "No data"
                continue

            # Filter gene names for best matches
            if response_df.iloc[:, 1].str.contains(gene, case=False).any():
                response_df_clean = response_df[response_df.iloc[:, 1].str.contains(
                    gene, case=False)].reset_index(drop=True)
            else:
                response_df_clean = response_df[(response_df.iloc[:, 1].str.contains(gene, case=False)) | (
                    response_df.iloc[:, 2].str.contains(gene, case=False))].reset_index(drop=True)

            # Append accession ID starting with "P" or second-best options
            if response_df_clean.Entry.str.startswith('P').any():
                p_hold = response_df_clean[response_df_clean.Entry.str.startswith(
                    'P')].reset_index(drop=True)
                accessionDict[gene] = p_hold.Entry[0]

            else:
                try:
                    accessionDict[gene] = response_df_clean.Entry[0]
                except KeyError:
                    print("No IDs for " + gene)
                    accessionDict[gene] = "No IDs"

        return accessionDict


    # Get and format API call for DGIdb
    def GetDGIdb(self, method, output_path, date):
        # Genes to analyze
        input_genes = self.consDf[method].dropna()

        # Prep parameters for API call
        # Genes/Drugs (required for call) - comma delimited list of gene or drug names/symbols
        dgidb_genes = ",".join(input_genes)

        # API Call
        DGIdb_response = requests.get(
            "https://dgidb.org/api/v2/interactions.json?genes=" + dgidb_genes)

        # Check for server error
        if DGIdb_response.status_code == 200:
            print("Response OK")
        elif DGIdb_response.status_code != 200:
            print("Error in DGIdb request: ", DGIdb_response.status_code)

        # Parse response into interactions by gene name
        # matchedTerms
        matchedTerms_DGIdb = pd.json_normalize(DGIdb_response.json(), record_path=["matchedTerms", "interactions"], meta=[
            ["matchedTerms", "geneName"], ["matchedTerms", "geneLongName"]])

        # ambiguousTerms
        ambiguousTerms_DGIdb = pd.json_normalize(
            DGIdb_response.json(), record_path=["ambiguousTerms"])

        # unmatchedTerms
        unmatchedTerms_DGIdb = pd.json_normalize(
            DGIdb_response.json(), record_path=["unmatchedTerms"])

        # Check for missing terms and write to output files
        matchedTerms_DGIdb.to_csv(output_path + "/" + 
                                    self.expName + "_DGIdb_matchedTerms_DrugInteractions_" + str(date) + ".csv", index=False)

        if ambiguousTerms_DGIdb.shape[0] == 0 and unmatchedTerms_DGIdb.shape[0] == 0:
            print("All search genes were matched in DGIdb")
        elif ambiguousTerms_DGIdb.shape[0] != 0 or unmatchedTerms_DGIdb.shape[0] != 0:
            print("Some search terms were not found in DGIdb")
            ambiguousTerms_DGIdb.to_csv(output_path + "/" +
                                        self.expName + "_DGIdb_ambiguousTerms_" + str(date) + ".csv", index=False)
            unmatchedTerms_DGIdb.to_csv(output_path + "/" +
                                        self.expName + "_DGIdb_unmatchedTerms_" + str(date) + ".csv", index=False)
        return matchedTerms_DGIdb


    # Get and format API call to BindingDB
    def GetBindingDB(self, method, output_path, date):
        # Get genes to search Uniprot_ids
        input_genes = self.consDf[method].dropna()

        # Get UniProtIDs for genes
        accessionDict = self.GetUniprotIDs(input_genes)

        # Prep accession list and make API Call
        BindingDB_accession = []
        for key, value in accessionDict.items():
            BindingDB_accession.append(value)
        BindingDB_accession_str = ",".join(BindingDB_accession)

        # Our filter - Genes, only interactions with < 200 nM, only compounds that are FDA approved
        BindingDB_response = requests.get("http://bindingdb.org/axis2/services/BDBService/getLigandsByUniprots?uniprot=" +
                                            BindingDB_accession_str + "&cutoff=200&code=2&response=application/json")
        if BindingDB_response.status_code != 200:
            print("Error in BindingDB API Call")

        # Format output for return
        matchedTerms_BindingDB = pd.json_normalize(BindingDB_response.json(), record_path=[
            'getLigandsByUniprotsResponse', 'affinities'])

        if matchedTerms_BindingDB.shape[0] == 0:
            print("No search terms found for ", method)

        else:
            # Convert accessionDict to df
            accessionDict_df = pd.DataFrame.from_dict(
                accessionDict, orient='index').reset_index()
            accessionDict_df.rename(
                columns={'index': 'gene', 0: 'accession_id'}, inplace=True)

            # Merge with matchedTerms
            try: allTerms_BindingDB = pd.merge(accessionDict_df, matchedTerms_BindingDB,
                                            left_on='accession_id', right_on='query', how='outer')
            except KeyError:
                print('No drug interactions found')

            # Write to .csv output and return
            allTerms_BindingDB.to_csv(output_path + "/" +
                                        self.expName + "_BindingDB_matchedTerms_DrugInteractions_" + str(date) + ".csv", index=False)

            return matchedTerms_BindingDB


    def main(self):
        # Call Functions
        date = datetime.date.today()
        pool = mp.Pool(self.Cores)
        for method in self.methodList:
            # Create new output
            os.makedirs(self.oPath + method, exist_ok= True)

            # Get drug gene interactions
            pool.apply(self.GetDGIdb, args = (method, self.oPath + method + "/", date))
            pool.apply(self.GetBindingDB, args = (method, self.oPath + method + "/", date))
        pool.close()
        pool.join()


