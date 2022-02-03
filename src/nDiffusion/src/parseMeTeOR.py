from IPython import embed

def parseGeneMap (fl = '/Users/minhpham/Desktop/Projects/MeTeOR_0418/data/geneMapping.txt'):
    entrezid_gene, gene_entrezid = {}, {}
    ids = []
    for line in open(fl).readlines():
        if 'MESHUI' not in line:
            line = line.strip('\n').split()
            entrezid_gene['EntrezID:{}'.format(line[1])] = line[2]
            gene_entrezid[line[2]] = 'EntrezID:{}'.format(line[1])
            ids.append('EntrezID:{}'.format(line[1]))
    return entrezid_gene, gene_entrezid, ids

def getMeTeOR(network_path = '/Users/minhpham/Desktop/Projects/MeTeOR_0418/results/network_com/MeTeORAll.txt', outfile = '../data/networks/MeTeOR.txt'):
          out = open(outfile, 'w')
          unmapped = []
          entrez_gene, gene_entrez, entrezs = parseGeneMap()
          for line in open(network_path).readlines():
                    line = line.strip('\n').split('\t')
                    line1 = []
                    for i in line:
                              if 'EntrezID:' in i:
                                        try:
                                                  mapped = entrez_gene[i]
                                                  line1.append(mapped)
                                        except:
                                                  line1.append(i)
                                                  unmapped.append(i)
                              else:
                                        line1.append(i)
                    line1_txt = '\t'.join(str(x) for x in line1)
                    out.write('{}\n'.format(line1_txt))
          out.close()
          unmapped = set(unmapped)
          return unmapped

meteor_unmapped = getMeTeOR()
embed()