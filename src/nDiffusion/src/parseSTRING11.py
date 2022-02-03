from IPython import embed

def getSTRINGv11(network_path = '/Users/minhpham/Desktop/Projects/Databases/STRING10/', outfile = '../data/networks/STRINGv11.txt'):
          out = open(outfile, 'w')
          
          string_entrez = {}
          def entrezString(fl):
                    for line in open(fl).readlines():
                              if '#' not in line:
                                        line = line.strip('\n').split()
                                        string_entrez[line[1].replace('9606.','')] = line[0]
          entrezString(network_path+'entrez_gene_id.vs.string.v10.28042015.tsv')
          entrezString(network_path+'entrez_gene_id.vs.string.v9.05.28122012.txt')
          
          def human_entrez_string (fl = network_path + 'human.entrez_2_string.2018.tsv'):
                    for line in open(fl).readlines():
                              if '#' not in line:
                                        line = line.strip('\n').split()
                                        string_entrez[line[2].replace('9606.','')] = line[1]
          human_entrez_string()

          entrez_name = {}
          for line in open(network_path+'HGNC_3-20-18.txt').readlines():
                    if 'Symbol' not in line:
                              line = line.strip('\n').split()
                              try:
                                        entrez_name[line[1]] = line[0]
                              except:
                                        pass
          
          mapping_fl = network_path+'all_go_knowledge_explicit.tsv'
          id_name = {}
          for line in open(mapping_fl).readlines():
                    line = line.strip('\n').split()
                    id_name[line[1]] = line[2]
          id_name['ENSP00000365759'] = 'DDR1'
          id_name['ENSP00000407685'] = 'TNXB'
          id_name['ENSP00000370781'] = 'CDC20B'
          id_name['ENSP00000353245'] = 'C17orf97'
          id_name['ENSP00000366348'] = 'OR12D3'
          id_name['ENSP00000265299'] = 'FAM188B'
          id_name['ENSP00000381760'] = 'IQCF6'
          id_name['ENSP00000351052'] = 'ZC3HC1'
          id_name['ENSP00000322898'] = 'EBF1'
          id_name['ENSP00000349437'] = 'IGF2R'
          id_name['ENSP00000483018'] = 'ABO'
          def human_name_string (fl = network_path + 'human.name_2_string.tsv'):
                    for line in open(fl).readlines():
                              if '#' not in line:
                                   line = line.strip('\n').split()
                                   id_name[line[2].replace('9606.', '')] = line[1]
          human_name_string()
          def protein_info (fl = network_path + '9606.protein.info.v11.0.txt'):
                    for line in open(fl).readlines():
                              if 'preferred_name' not in line:
                                   line = line.strip('\n').split()
                                   id_name[line[0].replace('9606.', '')] = line[1]
          protein_info()
          # Parsing STRING network
          net_name = 'STRING11_combined'
          net_fl = network_path+'9606.protein.links.detailed.v11.0.txt'
          net_lst, unmapped, mapped = [], [], []
          for line in open(net_fl).readlines():
                    if 'combined_score' not in line:
                              line = line.strip('\n').split()
                              sums = float(line[9])
                              try:
                                        protein1 = id_name[line[0].replace('9606.', '')]
                                        mapped.append(line[0].replace('9606.', ''))
                              except:
                                        try:
                                                  protein1 = entrez_name[string_entrez[line[0].replace('9606.', '')]]
                                                  mapped.append(line[0].replace('9606.', ''))
                                        except:
                                                  protein1 = line[0].replace('9606.', '')
                                                  unmapped.append(protein1)
                              try:
                                        protein2 = id_name[line[1].replace('9606.', '')]
                                        mapped.append(line[1].replace('9606.', ''))
                              except:
                                        try:
                                                  protein2 = entrez_name[string_entrez[line[1].replace('9606.', '')]]
                                                  mapped.append(line[1].replace('9606.', ''))
                                        except:
                                                  protein2 = line[1].replace('9606.', '')
                                                  unmapped.append(protein2)
                              out.write('{}\t{}\t{}\n'.format(protein1, protein2, sums))
                              net_lst.append('{}\t{}\t{}'.format(protein1, protein2, sums))
          unmapped = list(set(unmapped))
          mapped = list(set(mapped))
          print("{} nodes in the network are mapped to genes, {} are not mapped".format(len(mapped), len(unmapped)))
          out.close()

getSTRINGv11()
