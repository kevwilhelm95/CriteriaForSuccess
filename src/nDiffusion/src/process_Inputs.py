'''@author: minhpham'''

import networkx as nx
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, coo_matrix, csgraph, identity
from collections import Counter
import random

def getGraph (net_lst = '../data/networks/STRINGv10.txt'):
    G = nx.read_edgelist(open(net_lst, 'r'), data=(('weight',float),))
    graph_node = list(G.nodes())
    adjMatrix = nx.to_scipy_sparse_matrix(G)
    node_degree = dict(nx.degree(G))
    G_degree = node_degree.values()
    return G, graph_node, adjMatrix, node_degree, G_degree

def getDiffusionParam (adjMatrix):
    L = csgraph.laplacian(adjMatrix, normed=True)
    n = adjMatrix.shape[0]
    I = identity(n, dtype='int8', format='csr')
    axisSum = coo_matrix.sum(np.abs(L), axis=0)
    sumMax = np.max(axisSum)
    diffusionParameter = (1 / float(sumMax))
    ps = (I + (diffusionParameter * L))
    return ps

def readInput (fl):
    lst = []
    for line in open(fl).readlines():
        line = line.strip('\n').strip('\r').split('\t')
        lst.append(line[0])
    return lst

def getIndexdict (graph_node):
    graph_node_index = {}
    for i in range(len(graph_node)):
        graph_node_index[graph_node[i]] = i
    return graph_node_index

def getIndex (lst, graph_node_index):
    index = []
    for i in lst:
        ind = graph_node_index[i]
        index.append(ind)
    return index

def getDegree (pred_node, node_degree):
    pred_degree = []
    for i in pred_node:
        pred_degree.append(node_degree[i])
    pred_degree_count = dict(Counter(pred_degree))
    return pred_degree_count

def parseGeneInput (fl1, fl2, graph_node, graph_node_index, node_degree, graph_gene=[]):
    ### Parsing input files
    group1 = set(readInput(fl1))
    group2 = set(readInput(fl2))
    fl1_name = fl1.split('/')[-1].split('.')[0]
    fl2_name = fl2.split('/')[-1].split('.')[0]
    overlap = list(set(group1).intersection(group2))
    group1_only = list(set(group1)-set(overlap))
    group2_only = list(set(group2)-set(overlap))
    ### Mapping genes into the network
    group1_node = list(set(group1).intersection(graph_node))
    group2_node = list(set(group2).intersection(graph_node))
    overlap_node = list(set(overlap).intersection(graph_node))
    if graph_gene == []:
        other = list(set(graph_node) - set(group1_node) - set(group2_node))
    else:
        other = graph_gene
    group1_only_node = list(set(group1_node)-set(overlap_node))
    group2_only_node = list(set(group2_node)-set(overlap_node))
    print("{} genes are mapped (out of {}) in {}\n {} genes are mapped (out of {}) in {}\n {} are overlapped and mapped (out of {})\n".format(len(group1_node), len(group1), fl1_name, len(group2_node), len(group2), fl2_name, len(overlap_node), len(overlap)))
    ### Getting indexes of the genes in the network node list
    group1_only_index = getIndex(group1_only_node, graph_node_index)
    group2_only_index = getIndex(group2_only_node, graph_node_index)
    overlap_index = getIndex(overlap_node, graph_node_index)
    other_index = list(set(range(len(graph_node))) - set(group1_only_index) - set(group2_only_index)-set(overlap_index))
    ### Getting counter dictionaries for the connectivity degrees of the genes
    group1_only_degree_count = getDegree(group1_only_node, node_degree)
    group2_only_degree_count = getDegree(group2_only_node, node_degree)
    overlap_degree_count = getDegree(overlap_node, node_degree)
    ### Combining these features into dictionaries
    GP1_only_dict={'orig': group1_only, 'node':group1_only_node, 'index':group1_only_index, 'degree': group1_only_degree_count}
    GP2_only_dict={'orig': group2_only,'node':group2_only_node, 'index':group2_only_index, 'degree': group2_only_degree_count}
    overlap_dict={'orig': overlap, 'node':overlap_node, 'index':overlap_index, 'degree': overlap_degree_count}
    other_dict={'node':other, 'index':other_index} 
    
    return GP1_only_dict, GP2_only_dict, overlap_dict, other_dict 

def combineGroup (gp1_dict, gp2_dict):
    combine_dict = {}
    combine_dict['orig'] = gp1_dict['orig']+gp2_dict['orig']
    combine_dict['node'] = gp1_dict['node']+gp2_dict['node']
    combine_dict['index'] = gp1_dict['index']+gp2_dict['index']
    combine_dict['degree'] = mergeDegreeDict(gp1_dict['degree'], gp2_dict['degree'])
    return combine_dict

def getDegreeNode(G_degree, node_degree, other):
    degree_nodes = {}
    for i in set(G_degree):
        degree_nodes[i] = []
        for y in node_degree:
            if node_degree[y] == i and y in other:
                degree_nodes[i].append(y)
        degree_nodes[i] = list(set(degree_nodes[i]))
        random.shuffle(degree_nodes[i])
    return degree_nodes

def mergeDegreeDict (dict1, dict2):
    merge_dict = {}
    for k in dict1:
        try:
            merge_dict[k] = dict1[k] + dict2[k]
        except:
            merge_dict[k] = dict1[k]
    for k in dict2:
        try:
            n = dict1[k]
        except:
            merge_dict[k] = dict2[k]
    return merge_dict
