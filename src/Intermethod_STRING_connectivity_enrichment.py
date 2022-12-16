"""
@author: Jenn Asmussen
    - Adapted for pyCFS by Kevin Wilhelm

Script to perform intermethod STRING gene connectivity analysis with degree matched bootstrapping

"""

import sys
import os
import networkx as nx
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from scipy.stats import norm, hypergeom
from venn import venn
from helper_functions import CreateDir

class InterMethod_Connectivity():
    def __init__(self, gold_standards, evidences, confidence, df, df_name, interstList, outpath):
        self.GS = gold_standards
        self.evidences = evidences.split(',')
        self.confidence = confidence
        self.df = df
        self.df_name = df_name
        self.interstList = interstList
        self.outpath = outpath
        self.main()

    
    def selectEvidences(self, evidence_lst, network):
        return network[['node1', 'node2'] + evidence_lst]

    def getEvidenceTypes(self, evidence_lst):
        evidence_lst = [x for x in evidence_lst if x != 'None']
        if 'all' in evidence_lst:
            evidence_lst = ['neighborhood','fusion','cooccurence','coexpression','experimental','database','textmining']
        else:
            pass
        return evidence_lst


    def getCombinedScore(self, net_df): #as calculated by STRING network
        cols = net_df.columns.values.tolist()
        cols = cols[2:]
        p = 0.041
        for col in cols:
            net_df[col] = 1 - ((net_df[col]/1000) - p) / (1 - p)
        for col in cols:
            net_df[col] = np.where(net_df[col] == 1 - ((0 - p) / (1 - p)), 1, net_df[col])
        net_df['score'] = 1 - np.product([net_df[cols[i]] for i in range(len(cols))], axis=0)
        net_df['score'] = net_df['score'] + p * (1 - net_df['score'])
        return net_df


    def getEdgeWeight(self, edge_confidence):
        if edge_confidence == 'low':
            weight = 0.2
        elif edge_confidence == 'high':
            weight = 0.7
        elif edge_confidence == 'highest':
            weight = 0.9
        elif edge_confidence == 'all':
            weight = 0.0
        else:
            weight = 0.4
        return weight


    def getGeneSources(self, set_dict, savepath, status):
        geneSource = {}
        for k, v in set_dict.items():
            for protein in set_dict[k]:
                if protein in geneSource:
                    geneSource[protein].append(k)
                else:
                    geneSource[protein] = [k]

        if status == 'randomSets':
            pass
        else:
            df = pd.DataFrame.from_dict(geneSource, orient = 'index')
            df.to_csv(savepath + 'GeneSet_Sources.csv', header = False)
        return geneSource

    def getUniqueGenes(self, source_dict):
        uniqueGenes = {}
        for k, v in source_dict.items():
            if len(v) > 1:
                pass
            else:
                uniqueGenes[k] = v
        return uniqueGenes

    def checkConnection(self, x, y, gene_lst):
        if x in gene_lst and y in gene_lst:
            out = 'yes'
        else:
            out = 'no'
        return out

    def getNodePair(self, x, y):
        pair = [x, y]
        pair.sort()
        return str(pair)

    def getUniqueGeneNetwork(self, gene_lst, network):
        n_df = network.copy()
        n_df1 = n_df[n_df['node1'].isin(gene_lst)]
        n_df2 = n_df[n_df['node2'].isin(gene_lst)]
        n_df_final = pd.concat([n_df1, n_df2])
        n_df_final['bwSetConnection'] = n_df_final.apply(lambda x: checkConnection(x['node1'], x['node2'], gene_lst),
                                                        axis=1)
        n_df_final = n_df_final[n_df_final['bwSetConnection'] == 'yes']
        n_df_final.drop_duplicates(inplace=True)
        n_df_final['pair'] = n_df_final.apply(lambda x: getNodePair(x['node1'], x['node2']), axis=1)
        n_df_final.sort_values(by='node1', inplace=True)
        n_df_final.drop_duplicates(subset=['pair'], inplace=True, keep='first')
        return n_df_final

    def getUniqueGeneNetwork_BwMethodConnections(self, network, unique_dict):
        net_genes = list(set(network['node1'].tolist() + network['node2'].tolist()))
        bw_method_edges = []

        for g in net_genes:
            df1 = network[network['node1'] == g]
            df1_genes = df1['node2'].tolist()

            for p in df1_genes:
                if unique_dict[p][0] == unique_dict[g][0]:
                    pass
                else:
                    bw_method_edges.append([g, p])

            df2 = network[network['node2'] == g]
            df2_genes = df2['node1'].tolist()

            for p in df2_genes:
                if unique_dict[p][0] == unique_dict[g][0]:
                    pass
                else:
                    bw_method_edges.append([g, p])

        bw_method_edges = [sorted(x) for x in bw_method_edges]
        bw_method_edges = [str(x) for x in bw_method_edges]
        bw_method_edges = list(set(bw_method_edges))

        return bw_method_edges

    def getUniqueGeneCounts(self, uniqueGeneDict):
        uniqueGenesPerSet = {}
        for k, v in uniqueGeneDict.items():
            if v[0] in uniqueGenesPerSet:
                uniqueGenesPerSet[v[0]].append(k)
            else:
                uniqueGenesPerSet[v[0]] = [k]
        return uniqueGenesPerSet

    def getNodeDegree(self, node, degree_dict):
        return degree_dict[node]

    def getNodeDegreeDict(self, uniqueGenes, degree_df):
        df = degree_df.copy()
        df = df.loc[uniqueGenes]
        return df

    def createRandomDegreeMatchedSet(self, uniqueGeneSets, stringNet_allGenes, stringNet_degree_df):
        randomSets = {}
        for k, v in uniqueGeneSets.items():
            uniqueMappedGenes = v
            uniqueMappedGenes = [x for x in uniqueMappedGenes if x in stringNet_allGenes] #need to filter for genes that are mapped to STRING appropriately
            uniqueMappedGenes_degree_df = getNodeDegreeDict(uniqueMappedGenes, stringNet_degree_df)
            uniqueMappedGenes_degree_df = pd.DataFrame(uniqueMappedGenes_degree_df.groupby('degree_rounded')['degree'].count())
            #in this dictionary: key is degree, value is count of genes with that degree
            uniqueMappedGenes_degree_dict = dict(zip(uniqueMappedGenes_degree_df.index.tolist(), uniqueMappedGenes_degree_df['degree'].tolist()))

            randomGenes = []
            for k1, v1 in uniqueMappedGenes_degree_dict.items():
                degree_df = stringNet_degree_df.copy()
                degree_df = degree_df[degree_df['degree_rounded']==k1]

                degreeMatchedGenes = degree_df.index.tolist()
                x = np.random.choice(degreeMatchedGenes, v1, replace=False).tolist()
                randomGenes.extend(x)
            #print(len(uniqueMappedGenes), len(randomGenes)) #these numbers should match
            randomSets[k] = randomGenes
        return randomSets

    def plotResultsNorm(self, randomDist, trueConnections, savepath):
        fig, ax = plt.subplots(tight_layout = True)
        ax.hist(randomDist, color='grey', density = True)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.axvline(np.mean(randomDist), linestyle='dashed', color='black')
        plt.axvline(len(trueConnections), linestyle='dashed', color='red')
        plt.xlabel('Number Intermethod Edges', size=14)
        plt.ylabel('Random Analysis Count', size=14)

        mu, sigma = norm.fit(randomDist)
        z_fit = (len(trueConnections) - mu)/sigma
        xmin, xmax = plt.xlim()
        plt.xlim(xmin, xmax)
        plt.text(xmin + 5, 0.002, 'Z-score: ' + str(round(z_fit, 2)) + '\nNumber True Connections: ' + str(len(trueConnections)),
                bbox={'facecolor': 'blue', 'alpha': 1.0, 'pad': 2}, color = 'white')

        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, sigma)
        ax.plot(x, p, linewidth = 2)
        plt.savefig(savepath + 'Result_Plot.pdf', transparent=True)
        plt.show()
        plt.close()

    def outputRandomSets(self, randomSets_Connections):
        with open(args.savepath + 'RandomSetConnections.txt','w') as f:
            for i in randomSets_Connections:
                f.writelines(str(i) + '\n')

    def outputUniqueGeneNetwork(self, uniqueGeneNetwork, uniqueGenes, savepath):
        uniqueGeneNetwork['node1_source'] = uniqueGeneNetwork['node1'].map(uniqueGenes)
        uniqueGeneNetwork['node2_source'] = uniqueGeneNetwork['node2'].map(uniqueGenes)
        uniqueGeneNetwork = uniqueGeneNetwork[uniqueGeneNetwork['node1_source'] != uniqueGeneNetwork['node2_source']]
        uniqueGeneNetwork = uniqueGeneNetwork[['node1','node2','score','node1_source','node2_source']]
        uniqueGeneNetwork.to_csv(savepath + 'UniqueGeneNetwork.csv', index = False)


    def hypergeoOverlapGSgenes(self, method_genes, gs_genes, overlap, bkgd_size):
        M = bkgd_size  # bkgd population size
        N = method_genes  # num draws from population
        n = gs_genes  # num successes in population
        k = overlap  # num success in draw
        pval = hypergeom.sf(k - 1, M, n, N)
        return pval

    def geneSetOverlap(self, geneSets, savepath):
        geneSets_ = {}
        for k, v in geneSets.items():
            geneSets_[k] = set(v)
        venn(geneSets_, fontsize=8, legend_loc="upper left")
        plt.savefig(savepath + 'GeneSet_VennDiagram.pdf', transparent = True)


    def main(self):
        #load and customize STRINGv11 network for analysis (evidence types, edge weight)
        cwd1_path = os.path.join(os.getcwd(), "..")
        stringNet = pd.read_csv(cwd1_path + "/refs/STRINGv11_ProteinNames_DetailedEdges_07062022.csv")
        stringNet_allGenes = list(set(stringNet['node1'].unique().tolist() + stringNet['node2'].unique().tolist()))
        evidence_lst = self.getEvidenceTypes(self.evidences)
        stringNet = self.selectEvidences(evidence_lst, stringNet)
        stringNet = self.getCombinedScore(stringNet)

        #Filtering network for edge weight
        edgeWeight = self.getEdgeWeight(self.confidence)
        stringNet = stringNet[stringNet['score']>= edgeWeight]

        #get degree connectivity after edgeweight filtering
        G_stringNet = nx.from_pandas_edgelist(stringNet[['node1', 'node2']], 'node1', 'node2') #network is already edge weight filtered
        G_stringNet_degree = dict(G_stringNet.degree)
        stringNet_degree_df = pd.DataFrame(index = stringNet_allGenes)
        stringNet_degree_df['degree'] = stringNet_degree_df.index.map(G_stringNet_degree)
        stringNet_degree_df.fillna(0, inplace = True) #fillna with zeros for genes w/o appropriate edge weights
        stringNet_degree_df['degree_rounded'] = stringNet_degree_df['degree'].apply(lambda x : round(x/10)*10) #round each degree to nearest 10

        #true gene sets b/w set unique gene connectivity
        # Here can we add between methods???
        # Need to loop through dfs and gs's
        for gL in self.df.columns:
            # Get query genes
            query = self.df[gL].dropna()
            for gS in self.GS.columns:
                # Get gold standard gene set
                ref_genes = self.GS[gS].dropna()

                # Create output path 
                newOutPutPath = CreateDir(self.outpath, f"{gL}/{gS}_Evidences-{','.join(evidence_lst)}_Confidence-{self.confidence}/")

                geneSets = {gL:list(query), gS:list(ref_genes)}
                geneSources = self.getGeneSources(geneSets, newOutPutPath, 'trueSets')
                uniqueGenes = self.getUniqueGenes(geneSources)
                uniqueGeneNetwork = self.getUniqueGeneNetwork(list(uniqueGenes.keys()), stringNet)
                trueConnections = self.getUniqueGeneNetwork_BwMethodConnections(uniqueGeneNetwork, uniqueGenes)
                #print('Number of b/w set connects in real gene sets:', len(trueConnections))

                #ouput true gene set b/w set network and venn diagram of overlapping genes
                self.outputUniqueGeneNetwork(uniqueGeneNetwork, uniqueGenes, newOutPutPath)
                self.geneSetOverlap(geneSets, newOutPutPath)

                #random gene sets b/w set unique gene connectivity w/ degree matching
                uniqueGeneSets = self.getUniqueGeneCounts(uniqueGenes) #dictionary with unique genes (values) per each set (keys) in true gene lists

                randomSets_Connections = []
                for i in range(250):
                    randomSets = self.createRandomDegreeMatchedSet(uniqueGeneSets, stringNet_allGenes, stringNet_degree_df)
                    geneSources = self.getGeneSources(randomSets, newOutPutPath, 'randomSets')
                    uniqueGenes = self.getUniqueGenes(geneSources)
                    uniqueGeneNetwork = self.getUniqueGeneNetwork(list(uniqueGenes.keys()), stringNet)
                    randomConnections = self.getUniqueGeneNetwork_BwMethodConnections(uniqueGeneNetwork, uniqueGenes)
                    randomSets_Connections.append(len(randomConnections))

                mu, sigma = norm.fit(randomSets_Connections)
                z_fit = (len(trueConnections) - mu)/sigma

                self.plotResultsNorm(randomSets_Connections, trueConnections, newOutPutPath)
                self.outputRandomSets(randomSets_Connections)