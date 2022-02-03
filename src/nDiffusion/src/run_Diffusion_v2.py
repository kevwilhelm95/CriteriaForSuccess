'''@author: minhpham'''
from process_Inputs import *
from utils Íimport *
from writeSummary import *
from IPython import embed
import sys, os

# Parse arguments for handling 
def parse_args():
    """
    Parses the Big pipeline arguments.
    """
    parser = argparse.ArgumentParser(description="nDiffusion arguments")
    parser.add_argument('--network', nargs='?', default = './src/nDiffusion/data/networks/toy_network.txt')## need to figure out the network part here
    parser.add_argument('--gl_1', nargs = '?', default = './', help='first input gene list')
    parser.add_argument('--gl_2', nargs = '?', default = './', help = 'second input gene list')
    parser.add_argument('--output', nargs = '?', default = './', help = 'output path for results')

    return parser.parse_args()

### CHANGE HERE: Assigning network, gene inputs, and result directory ###
network_fl = args.network
geneList1_fl = args.gl_1
geneList2_fl = args.gl_2
result_fl = args.output
group1_name = geneList1_fl.split('/')[-1].split('.')[0]
group2_name = geneList2_fl.split('/')[-1].split('.')[0]
repeat = 100 #number of randomization iterations

### For a multimodal network, specify graph genes
graph_gene = []
if network_fl == '../data/networks/MeTeOR.txt':
    for line in open(network_fl).readlines():
        line = line.strip('\n').split('\t')
        for i in line[:2]:
            if ':' not in i:
                graph_gene.append(i)

if __name__ == "__main__":
    ### Directory of the result folder
    result_fl_figure = result_fl + 'figures/'
    result_fl_raw = result_fl + 'raw_data/'
    result_fl_ranking = result_fl + 'ranking/'
    if not os.path.exists(result_fl):
        os.makedirs(result_fl)
    if not os.path.exists(result_fl_figure):
        os.makedirs(result_fl_figure)
    if not os.path.exists(result_fl_raw):
        os.makedirs(result_fl_raw)   
    if not os.path.exists(result_fl_ranking):
        os.makedirs(result_fl_ranking)      
    print('Running ...')

    ### Getting network and diffusion parameters
    G, graph_node, adjMatrix, node_degree, G_degree = getGraph(network_fl)
    ps = getDiffusionParam(adjMatrix)
    graph_node_index = getIndexdict(graph_node)
    GP1_only_dict, GP2_only_dict, overlap_dict, other_dict = parseGeneInput(geneList1_fl, geneList2_fl, graph_node, graph_node_index, node_degree, graph_gene)
    degree_nodes = getDegreeNode(G_degree, node_degree, other_dict['node'])

    # Combine exclusive genes and overlapped genes in each group, if there is an overlap
    if overlap_dict['node'] != []:
        GP1_all_dict = combineGroup(GP1_only_dict, overlap_dict)
        GP2_all_dict = combineGroup(GP2_only_dict, overlap_dict)
        Exclusives_dict = combineGroup(GP1_only_dict, GP2_only_dict)
   
    ### Diffusion experiments
    def getResults(gp1, gp2, result_fl, gp1_name, gp2_name, show = '', exclude=[]):
        auroc, z_auc, auprc, z_prc, pval = runrun(gp1, gp2, result_fl, gp1_name, gp2_name, show, degree_nodes, other_dict['node'], graph_node_index, graph_node, ps, exclude=exclude, repeat = repeat)
        return auroc, z_auc, auprc, z_prc, pval
    
    #### auroc, z-scores for auc, auprc, z-scores for auprc, KS pvals
    #### z-scores: from_degree, to_degree, from_uniform, to_uniform

    if overlap_dict['node'] != [] and GP1_only_dict['node'] != [] and GP2_only_dict['node'] != []: 
        # From group 1 exclusive to group 2 all:
        R_gp1o_gp2 = getResults(GP1_only_dict, GP2_all_dict, result_fl, group1_name+'Excl', group2_name, show = '__SHOW_1_')
        # From group 2 exclusive to group 1 all:
        R_gp2o_gp1 = getResults(GP2_only_dict, GP1_all_dict, result_fl, group2_name+'Excl', group1_name, show = '__SHOW_2_')     
        # From group 1 exclusive to group 2 exclusive:
        R_gp1o_gp2o = getResults(GP1_only_dict, GP2_only_dict, result_fl, group1_name+'Excl', group2_name+'Excl')
        # From group 2 exclusive to group 1 exclusive:
        R_gp2o_gp1o = getResults(GP2_only_dict, GP1_only_dict, result_fl, group2_name+'Excl', group1_name+'Excl')
        # From group 1 exclusive to the overlap
        R_gp1o_overlap = getResults(GP1_only_dict, overlap_dict, result_fl, group1_name+'Excl', 'Overlap')
        # From group 2 exclusive to the overlap
        R_gp2o_overlap = getResults(GP2_only_dict, overlap_dict, result_fl, group2_name+'Excl', 'Overlap')
        # From overlap to (group 1 exclusive and group 2 exlusive)
        R_overlap_exclusives = getResults(overlap_dict, Exclusives_dict, result_fl,'Overlap', 'Exclus')
        ### Write output
        writeSumTxt (result_fl, group1_name, group2_name, GP1_only_dict, GP2_only_dict, overlap_dict, R_gp1o_gp2, R_gp2o_gp1, R_gp1o_gp2o, R_gp2o_gp1o, R_gp1o_overlap, R_gp2o_overlap, R_overlap_exclusives)
    elif overlap_dict['node'] != [] and GP2_only_dict['node'] == []: #when group 2 is entirely part of group 1
        # From group 1 exclusive to overlap/group 2
        R_gp1o_overlap = getResults(GP1_only_dict, overlap_dict, result_fl, group1_name+'Excl', 'Overlap or'+group2_name)
        # From overlap/group 2 to group 1 exclusive
        R_overlap_gp1o = getResults(overlap_dict, GP1_only_dict, result_fl,'Overlap or'+group2_name, group1_name+'Excl')
        writeSumTxt (result_fl, group1_name, group2_name, GP1_only_dict, GP2_only_dict, overlap_dict, R_gp1o_overlap=R_gp1o_overlap, R_overlap_gp1o=R_overlap_gp1o)
    elif overlap_dict['node'] != [] and GP1_only_dict['node'] == []: #when group 1 is entirely part of group 2
        # From group 2 exclusive to overlap/group 1
        R_gp2o_overlap = getResults(GP2_only_dict, overlap_dict, result_fl, group2_name+'Excl', 'Overlap or '+group1_name)
        # From overlap/group 1 to group 2 exclusive
        R_overlap_gp2o = getResults(overlap_dict, GP2_only_dict, result_fl,'Overlap or'+group1_name, group2_name+'Excl')
        writeSumTxt (result_fl, group1_name, group2_name, GP1_only_dict, GP2_only_dict, overlap_dict, R_gp2o_overlap=R_gp2o_overlap, R_overlap_gp2o=R_overlap_gp2o)
    else: #when there is no overlap between two groups
        # From group 1 to group 2:
        R_gp1o_gp2o = getResults(GP1_only_dict, GP2_only_dict, result_fl, group1_name, group2_name, show = 'SHOW1')
        # From group 2 to group 1:
        R_gp2o_gp1o = getResults(GP2_only_dict, GP1_only_dict, result_fl, group2_name, group1_name, show = 'SHOW2')
        ### Write output
        writeSumTxt (result_fl, group1_name, group2_name, GP1_only_dict, GP2_only_dict, overlap_dict, R_gp1o_gp2o=R_gp1o_gp2o, R_gp2o_gp1o=R_gp2o_gp1o)

