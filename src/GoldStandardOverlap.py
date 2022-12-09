import pandas as pd
from scipy.stats import hypergeom as hg
from matplotlib import pyplot as plt
#from matplotlib_venn import venn2_unweighted
from matplotlib_venn import venn2
import os



# Define Class
class GoldStandardOverlap():
    ''' Call to run GoldStandardOverlap analysis and write output files
    Required Inputs:
    1) '''
    # __init__
    def __init__(self, gold_standards, df, df_name, num_genes, interst_list, experiment_name, output_path):
        self.goldStandards = gold_standards
        self.consDf = df
        self.consDfName = df_name
        self.numGenes = num_genes
        self.intrstList = interst_list
        self.expName = experiment_name
        self.oPath = output_path
        self.main()

    def main(self):
        # Sub functions
        # Calculates the number and p-value of Overlap
        def CalcOverlap(self, gL, gS):
            # Define holders
            gs_name_hold = []
            len_gs = []
            exp_name_hold = []
            len_method = []
            count_hold = []
            pval_hold = []
            over_gene_hold = []

            # Create hold dfs
            gL_hold = self.consDf[gL].dropna()
            gS_hold = self.goldStandards[gS].dropna()

            # Get overlap and calculate p-val
            overlap = list(gS_hold[gS_hold.isin(gL_hold)])
            pval = hg(M=self.numGenes, n=len(gS_hold), N=len(gL_hold)).sf(len(overlap) - 1)

            # Create return dict
            return_dict = {'GS': [gS], 
                            'Len(GS)': [len(gS_hold)],
                            'Method': [gL],
                            'Len(Method)': [len(gL_hold)],
                            'Count': [len(overlap)],
                            'P-value': [pval],
                            'Genes': [", ".join(list(overlap))]}
            return_df = pd.DataFrame(return_dict)

            # Make Venn Diagram of the results
            VennPlot(gS, gL, gS_hold, gL_hold, overlap, pval, self.oPath, self.expName) 

            # Return the dict
            return return_df


        # Plot venn diagram of overlap between lists
        def VennPlot(gs, exp, gs_hold, exp_hold, overlap, pval, outpath, experimental_name):
            # Make experimental list output directory
            outpath = outpath + exp + '/'
            os.makedirs(outpath, exist_ok = True)

            # Plot figure
            fig = plt.figure(figsize=(10, 5))
            out = venn2(subsets=((len(exp_hold) - len(overlap)),
                                            (len(gs_hold) - len(overlap)),
                                            len(overlap)),
                                set_labels=(exp, gs),
                                set_colors=('white', 'white'),
                                alpha=0.7)
            overlap1 = out.get_patch_by_id("A")
            overlap1.set_edgecolor("red")
            overlap1.set_linewidth(3)
            overlap2 = out.get_patch_by_id("B")
            overlap2.set_edgecolor('blue')
            overlap2.set_linewidth(3)

            for text in out.set_labels:
                text.set_fontsize(20)
            for text in out.subset_labels:
                text.set_fontsize(18)
            plt.text(0, -0.7,
                    str("p = " + str(round(pval, 8))),
                    horizontalalignment='center',
                    verticalalignment='top',
                    fontsize=14.0)
            plt.text(0, -0.78,
                    ", ".join(overlap),
                    horizontalalignment='center',
                    verticalalignment='top',
                    fontsize=14.0)
            fig.savefig(outpath + experimental_name + '_' + gs + '_' +
                        exp + '_VennDiagram.png', bbox_inches='tight')  # , dpi=150)
            plt.close()
        

        # Main function
        # Loop through genes lists of interest
        out_df = pd.DataFrame(columns=['GS', 'Len(GS)', 'Method', 'Len(Method)', 'Count', 'P-value', 'Genes'])
        for gL in self.intrstList:
            for gS in self.goldStandards.columns:
                # Get holding variables for function calls
                hold_out = CalcOverlap(self, gL, gS)
                out_df = pd.concat([out_df, hold_out], axis = 0, ignore_index=True)

        out_df.to_csv(self.oPath + self.expName + "_GoldStandardOverlap_Summary.csv", index=False)

