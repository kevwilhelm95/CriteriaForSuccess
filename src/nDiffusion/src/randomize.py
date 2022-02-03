'''@author: minhpham'''
from nDiffusion.src.process_Inputs import getIndex
from sklearn.metrics import roc_curve, precision_recall_curve, auc
import numpy as np
import random, os, sys
#sys.path.append('src/nDiffusion/srcß/multimodal-network-diffusion/')
from nDiffusion.src.MultimodalNetworkDiffusion.GraphBasedDiffusion import *

### Computing AUROC and AUPRC for each run 
def performance_run (from_index, to_index, graph_node, ps, exclude = [], diffuseMatrix = False):
    results = {}
    if exclude == []:
        exclude = from_index
    try:
        if diffuseMatrix == False:
            label = np.zeros(len(graph_node))
            for i in from_index:
                label[i] = 1
            diffuseMatrix = diffuse(label, ps)
    except:
        pass
    score, classify, scoreTP, gene_write = [], [], [], []
    for i in range(len(graph_node)):
        if i not in exclude:
            gene_write.append(graph_node[i])
            score.append(diffuseMatrix[i])
            if i in to_index:
                classify.append(1)
                scoreTP.append(diffuseMatrix[i])
            else:
                classify.append(0)
    results['classify'], results['score'], results['scoreTP'], results['genes'] = classify, score, scoreTP, gene_write
    results['diffuseMatrix'] = diffuseMatrix
    results['fpr'], results['tpr'], thresholds = roc_curve(classify, score, pos_label=1)
    results['auROC']= auc(results['fpr'], results['tpr'])
    results['precision'], results['recall'], thresholds = precision_recall_curve(classify, score, pos_label=1)
    results['auPRC'] = auc(results['recall'], results['precision'])
    return results 

###  Selecting random genes that have similar connectivity degrees in the network
def getRand_degree(pred_degree_count, degree_nodes, iteration=1):
    rand_node, rand_degree = [], {}
    for i in pred_degree_count:
        rand_degree[i] = []
        count = pred_degree_count[i]*iteration
        lst = []
        modifier = 1
        cnt = 0
        if float(i) <= 100:
            increment = 1
        elif float(i) <= 500:
            increment = 5
        else:
            increment = 10
        while len(lst) < count and modifier <= float(i)/10 and cnt <= 500:
            degree_select = [n for n in degree_nodes.keys() if n <= i+modifier and n >= i-modifier]
            node_select = []
            for m in degree_select:
                node_select += degree_nodes[m]
            node_select = list(set(node_select))
            random.shuffle(node_select)
            try:
                lst += node_select[0:(count-len(lst))]
            except:
                pass
            modifier += increment
            cnt += 1
            overlap = set(rand_node).intersection(lst)
            for item in overlap:
                lst.remove(item)
        rand_node += lst
        rand_degree[i] += lst
    return rand_node

### Selecting random genes uniformly
def getRand_uniform(pred_degree_count, other):
    number_rand = sum(pred_degree_count.values())
    rand_node = random.sample(other, number_rand)
    return rand_node

### Performing randomization
def runRand(node_degree_count, node_index, degree_nodes, other, graph_node_index, graph_node, ps, rand_type, node_type, diffuseMatrix=False, repeat=100):
    AUROCs, AUPRCs, scoreTPs = [], [], []
    if rand_type == 'uniform':
        getRand = getRand_uniform
        var2 = other
    elif rand_type == 'degree':
        getRand = getRand_degree
        var2 = degree_nodes
    for i in range(repeat):
        rand_node = getRand(node_degree_count, var2)
        rand_index = getIndex(rand_node, graph_node_index)
        if node_type == 'TO':
            results = performance_run(node_index, rand_index, graph_node, ps, diffuseMatrix=diffuseMatrix)
        elif node_type == 'FROM':
            results = performance_run(rand_index, node_index, graph_node, ps)
        AUROCs.append(results['auROC'])
        AUPRCs.append(results['auPRC'])
        scoreTPs += results['scoreTP']
    return AUROCs, AUPRCs, scoreTPs
