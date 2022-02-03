def writeSumTxt (result_fl, group1_name, group2_name, GP1_only_dict, GP2_only_dict, overlap_dict, R_gp1o_gp2=[], R_gp2o_gp1=[], R_gp1o_gp2o=[], R_gp2o_gp1o=[], R_gp1o_overlap=[], R_gp2o_overlap=[], R_overlap_exclusives=[], R_gp1_gp2=[], R_gp2_gp1=[], R_overlap_gp1o=[], R_overlap_gp2o=[]):
          #### KS test results
          ks_result = []
          ks_result.append(['Seeds','Recipients','Randomize Seeds (degree-matched)','Randomize Recipients (degree-matched)','Randomize Seeds (uniform)','Randomize Recipients (uniform)'])
          if overlap_dict['node'] != [] and GP1_only_dict['node'] != [] and GP2_only_dict['node'] != []: 
                    ks_result.append([group1_name+' Exclusive', group2_name, R_gp1o_gp2[4][0], R_gp1o_gp2[4][1], R_gp1o_gp2[4][2], R_gp1o_gp2[4][3]])
                    ks_result.append([group2_name+' Exclusive', group1_name, R_gp2o_gp1[4][0], R_gp2o_gp1[4][1], R_gp2o_gp1[4][2], R_gp2o_gp1[4][3]])
                    ks_result.append([group1_name+' Exclusive', group2_name+' Exclusive', R_gp1o_gp2o[4][0], R_gp1o_gp2o[4][1], R_gp1o_gp2o[4][2], R_gp1o_gp2o[4][3]])
                    ks_result.append([group2_name+' Exclusive', group1_name+' Exclusive', R_gp2o_gp1o[4][0], R_gp2o_gp1o[4][1], R_gp2o_gp1o[4][2], R_gp2o_gp1o[4][3]])                    
                    ks_result.append([group1_name+' Exclusive', 'Overlap', R_gp1o_overlap[4][0], R_gp1o_overlap[4][1], R_gp1o_overlap[4][2], R_gp1o_overlap[4][3]])
                    ks_result.append([group2_name+' Exclusive', 'Overlap', R_gp2o_overlap[4][0], R_gp2o_overlap[4][1], R_gp2o_overlap[4][2], R_gp2o_overlap[4][3]])
                    ks_result.append(['Overlap', 'Exclusives', R_overlap_exclusives[4][0], R_overlap_exclusives[4][1], R_overlap_exclusives[4][2], R_overlap_exclusives[4][3]])
                    try:
                              ks_result.append([group1_name, group1_name, R_gp1_gp2[4][0], R_gp1_gp2[4][1], R_gp1_gp2[4][2], R_gp1_gp2[4][3]])
                              ks_result.append([group2_name, group1_name, R_gp2_gp1[4][0], R_gp2_gp1[4][1], R_gp2_gp1[4][2], R_gp2_gp1[4][3]])
                    except:
                              pass
          elif overlap_dict['node'] != [] and GP2_only_dict['node'] == []:
                    ks_result.append([group1_name+' Exclusive', 'Overlap or '+group2_name, R_gp1o_overlap[4][0], R_gp1o_overlap[4][1], R_gp1o_overlap[4][2], R_gp1o_overlap[4][3]])
                    ks_result.append(['Overlap or '+group2_name, group1_name+' Exclusive', R_overlap_gp1o[4][0], R_overlap_gp1o[4][1], R_overlap_gp1o[4][2], R_overlap_gp1o[4][3]]) 
          elif overlap_dict['node'] != [] and GP1_only_dict['node'] == []:    
                    ks_result.append([group2_name+' Exclusive', 'Overlap or '+group1_name, R_gp2o_overlap[4][0], R_gp2o_overlap[4][1], R_gp2o_overlap[4][2], R_gp2o_overlap[4][3]])
                    ks_result.append(['Overlap or '+group1_name, group2_name+' Exclusive', R_overlap_gp2o[4][0], R_overlap_gp2o[4][1], R_overlap_gp2o[4][2], R_overlap_gp2o[4][3]]) 
          else:
                    ks_result.append([group1_name, group2_name, R_gp1o_gp2o[4][0], R_gp1o_gp2o[4][1], R_gp1o_gp2o[4][2], R_gp1o_gp2o[4][3]])
                    ks_result.append([group2_name, group1_name, R_gp2o_gp1o[4][0], R_gp2o_gp1o[4][1], R_gp2o_gp1o[4][2], R_gp2o_gp1o[4][3]]) 

          ### ROC results
          roc_result=[]
          roc_result.append(['Seeds','Recipients','AUROC','Z-score for Random Seeds (degree-matched)','Z-score for Random Recipients (degree-matched)','Z-score for Random Seeds (uniform)','Z-score for Random Recipients (uniform)'])
          if overlap_dict['node'] != [] and GP1_only_dict['node'] != [] and GP2_only_dict['node'] != []: 
                    roc_result.append([group1_name+' Exclusive', group2_name, R_gp1o_gp2[0], R_gp1o_gp2[1][0], R_gp1o_gp2[1][1], R_gp1o_gp2[1][2], R_gp1o_gp2[1][3]])
                    roc_result.append([group2_name+' Exclusive', group1_name, R_gp2o_gp1[0], R_gp2o_gp1[1][0], R_gp2o_gp1[1][1], R_gp2o_gp1[1][2], R_gp2o_gp1[1][3]])
                    roc_result.append([group1_name+' Exclusive', group2_name+' Exclusive', R_gp1o_gp2o[0], R_gp1o_gp2o[1][0], R_gp1o_gp2o[1][1], R_gp1o_gp2o[1][2], R_gp1o_gp2o[1][3]])
                    roc_result.append([group2_name+' Exclusive', group1_name+' Exclusive', R_gp2o_gp1o[0], R_gp2o_gp1o[1][0], R_gp2o_gp1o[1][1], R_gp2o_gp1o[1][2], R_gp2o_gp1o[1][3]])                    
                    roc_result.append([group1_name+' Exclusive', 'Overlap', R_gp1o_overlap[0], R_gp1o_overlap[1][0], R_gp1o_overlap[1][1], R_gp1o_overlap[1][2], R_gp1o_overlap[1][3]])
                    roc_result.append([group2_name+' Exclusive', 'Overlap', R_gp2o_overlap[0], R_gp2o_overlap[1][0], R_gp2o_overlap[1][1], R_gp2o_overlap[1][2], R_gp2o_overlap[1][3]])
                    roc_result.append(['Overlap', 'Exclusives', R_overlap_exclusives[0], R_overlap_exclusives[1][0], R_overlap_exclusives[1][1], R_overlap_exclusives[1][2], R_overlap_exclusives[1][3]])
                    try:
                              roc_result.append([group1_name, group2_name, R_gp1_gp2[0], R_gp1_gp2[1][0], R_gp1_gp2[1][1], R_gp1_gp2[1][2], R_gp1_gp2[1][3]])
                              roc_result.append([group2_name, group1_name, R_gp2_gp1[0], R_gp2_gp1[1][0], R_gp2_gp1[1][1], R_gp2_gp1[1][2], R_gp2_gp1[1][3]])
                    except:
                              pass
          elif overlap_dict['node'] != [] and GP2_only_dict['node'] == []:
                    roc_result.append([group1_name+' Exclusive', 'Overlap or '+group2_name, R_gp1o_overlap[0], R_gp1o_overlap[1][0], R_gp1o_overlap[1][1], R_gp1o_overlap[1][2], R_gp1o_overlap[1][3]])
                    roc_result.append(['Overlap or '+group2_name, group1_name+' Exclusive', R_overlap_gp1o[0], R_overlap_gp1o[1][0], R_overlap_gp1o[1][1], R_overlap_gp1o[1][2], R_overlap_gp1o[1][3]]) 
          elif overlap_dict['node'] != [] and GP1_only_dict['node'] == []:    
                    roc_result.append([group2_name+' Exclusive', 'Overlap or '+group1_name, R_gp2o_overlap[0], R_gp2o_overlap[1][0], R_gp2o_overlap[1][1], R_gp2o_overlap[1][2], R_gp2o_overlap[1][3]])
                    roc_result.append(['Overlap or '+group1_name, group2_name+' Exclusive', R_overlap_gp2o[0], R_overlap_gp2o[1][0], R_overlap_gp2o[1][1], R_overlap_gp2o[1][2], R_overlap_gp2o[1][3]]) 
          else:
                    roc_result.append([group1_name, group2_name, R_gp1o_gp2o[0], R_gp1o_gp2o[1][0], R_gp1o_gp2o[1][1], R_gp1o_gp2o[1][2], R_gp1o_gp2o[1][3]])
                    roc_result.append([group2_name, group1_name, R_gp2o_gp1o[0], R_gp2o_gp1o[1][0], R_gp2o_gp1o[1][1], R_gp2o_gp1o[1][2], R_gp2o_gp1o[1][3]]) 

          ### PRC results
          prc_result=[]
          prc_result.append(['Seeds','Recipients','AUPRC','Z-score for Random Seeds (degree-matched)','Z-score for Random Recipients (degree-matched)','Z-score for Random Seeds (uniform)','Z-score for Random Recipients (uniform)'])
          if overlap_dict['node'] != [] and GP1_only_dict['node'] != [] and GP2_only_dict['node'] != []: 
                    prc_result.append([group1_name+' Exclusive', group2_name, R_gp1o_gp2[2], R_gp1o_gp2[3][0], R_gp1o_gp2[3][1], R_gp1o_gp2[3][2], R_gp1o_gp2[3][3]])
                    prc_result.append([group2_name+' Exclusive', group1_name, R_gp2o_gp1[2], R_gp2o_gp1[3][0], R_gp2o_gp1[3][1], R_gp2o_gp1[3][2], R_gp2o_gp1[3][3]])
                    prc_result.append([group1_name+' Exclusive', group2_name+' Exclusive', R_gp1o_gp2o[2], R_gp1o_gp2o[3][0], R_gp1o_gp2o[3][1], R_gp1o_gp2o[3][2], R_gp1o_gp2o[3][3]])
                    prc_result.append([group2_name+' Exclusive', group1_name+' Exclusive', R_gp2o_gp1o[2], R_gp2o_gp1o[3][0], R_gp2o_gp1o[3][1], R_gp2o_gp1o[3][2], R_gp2o_gp1o[3][3]])                    
                    prc_result.append([group1_name+' Exclusive', 'Overlap', R_gp1o_overlap[2], R_gp1o_overlap[3][0], R_gp1o_overlap[3][1], R_gp1o_overlap[3][2], R_gp1o_overlap[3][3]])
                    prc_result.append([group2_name+' Exclusive', 'Overlap', R_gp2o_overlap[2], R_gp2o_overlap[3][0], R_gp2o_overlap[3][1], R_gp2o_overlap[3][2], R_gp2o_overlap[3][3]])
                    prc_result.append(['Overlap', 'Exclusives', R_overlap_exclusives[2], R_overlap_exclusives[3][0], R_overlap_exclusives[3][1], R_overlap_exclusives[3][2], R_overlap_exclusives[3][3]])
                    try:
                              prc_result.append([group1_name, group2_name, R_gp1_gp2[2], R_gp1_gp2[3][0], R_gp1_gp2[3][1], R_gp1_gp2[3][2], R_gp1o_gp2[3][3]])
                              prc_result.append([group2_name, group1_name, R_gp2_gp1[2], R_gp2_gp1[3][0], R_gp2_gp1[3][1], R_gp2_gp1[3][2], R_gp2o_gp1[3][3]])  
                    except:
                              pass
          elif overlap_dict['node'] != [] and GP2_only_dict['node'] == []:
                    prc_result.append([group1_name+' Exclusive', 'Overlap or '+group2_name, R_gp1o_overlap[2], R_gp1o_overlap[3][0], R_gp1o_overlap[3][1], R_gp1o_overlap[3][2], R_gp1o_overlap[3][3]])
                    prc_result.append(['Overlap or '+group2_name, group1_name+' Exclusive', R_overlap_gp1o[2], R_overlap_gp1o[3][0], R_overlap_gp1o[3][1], R_overlap_gp1o[3][2], R_overlap_gp1o[3][3]]) 
          elif overlap_dict['node'] != [] and GP1_only_dict['node'] == []:    
                    prc_result.append([group2_name+' Exclusive', 'Overlap or '+group1_name, R_gp2o_overlap[2], R_gp2o_overlap[3][0], R_gp2o_overlap[3][1], R_gp2o_overlap[3][2], R_gp2o_overlap[3][3]])
                    prc_result.append(['Overlap or '+group1_name, group2_name+' Exclusive', R_overlap_gp2o[2], R_overlap_gp2o[3][0], R_overlap_gp2o[3][1], R_overlap_gp2o[3][2], R_overlap_gp2o[3][3]]) 
          else:
                    prc_result.append([group1_name, group2_name, R_gp1o_gp2o[2], R_gp1o_gp2o[3][0], R_gp1o_gp2o[3][1], R_gp1o_gp2o[3][2], R_gp1o_gp2o[3][3]])
                    prc_result.append([group2_name, group1_name, R_gp2o_gp1o[2], R_gp2o_gp1o[3][0], R_gp2o_gp1o[3][1], R_gp2o_gp1o[3][2], R_gp2o_gp1o[3][3]]) 

          ### Mapping results
          gene_result=[]
          gene_result.append(['**','#Total', '# Mapped in the network','Not mapped genes',' '])
          if overlap_dict['node'] != []: 
                    gene_result.append([group1_name+' Exclusive', len(GP1_only_dict['orig']), len(GP1_only_dict['node']),
                                        ';'.join(str(x) for x in set(GP1_only_dict['orig']).difference(GP1_only_dict['node']))])
                    gene_result.append([group2_name+' Exclusive', len(GP2_only_dict['orig']), len(GP2_only_dict['node']),
                                        ';'.join(str(x) for x in set(GP2_only_dict['orig']).difference(GP2_only_dict['node']))])
                    gene_result.append(['Overlap', len(overlap_dict['orig']), len(overlap_dict['node']),
                                        ';'.join(str(x) for x in set(overlap_dict['orig']).difference(overlap_dict['node']))])
          else:
                    gene_result.append([group1_name, len(GP1_only_dict['orig']), len(GP1_only_dict['node']),
                                        ';'.join(str(x) for x in set(GP1_only_dict['orig']).difference(GP1_only_dict['node']))])
                    gene_result.append([group2_name, len(GP2_only_dict['orig']), len(GP2_only_dict['node']),
                                        ';'.join(str(x) for x in set(GP2_only_dict['orig']).difference(GP2_only_dict['node']))])

          ### Write files   
          file_hand = open('{}/diffusion_result.txt'.format(result_fl), 'w')  
          file_hand.write('Comparing distributions of experimental and random diffusion values (p-values for KS tests)\n')     
          for line in ks_result:
                    val = "\t".join(str(v) for v in line)
                    file_hand.writelines("%s\n" % val)
          file_hand.write('\n')
          file_hand.write('\n')
          file_hand.write('\n')
          file_hand.write('Evaluating how well {} genes are linked to {} genes, comparing against random\n'.format(group1_name, group2_name))
          file_hand.write('\n')
          file_hand.write('ROC results\n') 
          for line in roc_result:
                    val = "\t".join(str(v) for v in line)
                    file_hand.writelines("%s\n" % val)
          file_hand.write('\n')
          file_hand.write('PRC results\n') 
          for line in prc_result:
                    val = "\t".join(str(v) for v in line)
                    file_hand.writelines("%s\n" % val) 
          file_hand.write('\n')         
          file_hand.write('Z-scores are computed for the experimental area under ROC or PRC based on distributions of the random areas under these curves\n')
          file_hand.write('Seeds: Genes where diffusion signal starts FROM)\nRecipients: Genes that receive the diffusion signal and that are in the other validated group\n')
          file_hand.write('Random genes are selected either uniformly or degree matched\n')
          file_hand.write('\n')
          file_hand.write('\n')
          file_hand.write('\n')
          file_hand.write('Number of genes\n')
          for gene in gene_result:
                    val = "\t".join(str(v) for v in gene)
                    file_hand.writelines("%s\n" % val)
          file_hand.close()
