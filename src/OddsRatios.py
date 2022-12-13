import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np
from adjustText import adjust_text
import os
import scipy.stats as stats
import datetime
import matplotlib.pyplot as plt
import multiprocessing as mp
import math
from pathlib import Path
import subprocess
from helper_functions import *

# Define odds ratio class
class GetOddsRatios():
    # Define init
    def __init__(self, df, df_name, interstList, CaseControl, CaseControl_path, ref, VCF, exp_name, outpath, cores):
        self.consDf = df
        self.consDfName = df_name
        self.expName = exp_name
        self.oPath = outpath
        self.CaseControl = CaseControl
        self.CaseControl_path = CaseControl_path
        self.ref = ref
        self.VCF_path = VCF
        self.methodList = interstList
        self.Cores = cores

        self.main()

    # Run ExactTest script to get variant counts in cases and controls
    def RunExactTest(self):
        print("---Preparing intermediate files for ExactTest---")
        # Create intermediate files needed for ExactTest script
        main_outpath = os.path.abspath(os.path.join(self.oPath, "../.."))
        intermediate_outpath = CreateDir(main_outpath, "IntermediateFiles")
        print(intermediate_outpath)
        CreateSampleOnlyFile(self.CaseControl, intermediate_outpath)
        CreateSampleFamFile(self.CaseControl, intermediate_outpath)
        if "Genes" in self.consDf.columns:
            CreateGeneRegionFile(self.consDf.Genes, self.ref, intermediate_outpath)
        else: 
            CreateGeneRegionFile(self.consDf.AllUnique, self.ref, intermediate_outpath)

        exactTest_sh = f"{os.getcwd()}/ExactTest.sh"
        # Run ExactTest.sh script and wait for output
        cmd = [exactTest_sh, 
                self.VCF_path, 
                intermediate_outpath + "AllUniqueGenesLocationFile.txt",
                self.CaseControl_path,
                f"{intermediate_outpath}/CaseControl_SampleOnly.txt",
                f"{intermediate_outpath}/CaseControl_fam.fam",
                intermediate_outpath]
        #proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = False)
        print("---ExactTest Script Started---")
        subprocess.Popen(cmd, shell = False, check = True)

        #stdout, stderr = proc.communicate()
        #print(stdout)
        #print(stderr)

        # Load outputback in as self.exactTest
        print("---ExactTest Script Finished ---")
        self.exactTest = pd.read_csv(f"{intermediate_outpath}/CaseControl.Variants.OR.txt",
                                    sep='\t', header=None,
                                     names=['chrom', 'pos', 'ref', 'alt', 'gene', 'ENSP', 'Consequence','HGVSp', 'EA', 'AN_0', 'AN_1', 'Cases', 'Controls', 'AC_1', 'AC_Het_1', 'AC_Hom_1', 'AC_0', 'AC_Het_0', 'AC_Hom_0', 'CC_ALL', 'CC_DOM', 'CC_REC', 'AF', 'AC'], low_memory=False)
         # Clean the file now
        self.exactTest = self.exactTest[self.exactTest.Consequence.str.contains("frameshift_variant|missense_variant|stop_gained|stop_lost", case = False)].reset_index() # Need to include splice sites
        self.exactTest['EA-Clean'] = [x.split(',')[0] for x in self.exactTest['EA']]
        self.exactTest['EA-Clean'] = [-1 if x == '.' else x for x in self.exactTest['EA-Clean']]
        self.exactTest['EA-Clean'] = self.exactTest['EA-Clean'].astype(float)


    # Clean input files to format for program - self.totalCases, self.totalControls
    def CleanInputs(self, method):
        # Get total number of cases and controls
        self.totalCases = self.CaseControl[self.CaseControl[1]== 1].shape[0]
        self.totalControls = self.CaseControl[self.CaseControl[0]
                                                == 0].shape[0]

        # Further subset exactTest into gene-list specific tables
        hold_genes = list(self.consDf[method].dropna())
        exactTest_method = self.exactTest[self.exactTest.gene.isin(
            hold_genes)]
        #exactTest_method['EA-Clean'].mask(exactTest_method['EA-Clean'] == '.', np.nan, inplace = True)

        # Split exactTest into rare and common variants
        variants_rare = exactTest_method[exactTest_method.AF <= 0.01]
        variants_common = exactTest_method[exactTest_method.AF >= 0.01]

        # Print how many variants are in the analysis
        print("Gene list: ", method)
        print("Total Variants: ", exactTest_method.shape[0])
        print("Common Variants: ", variants_common.shape[0])
        print("Rare Variants: ", variants_rare.shape[0])

        return variants_rare, variants_common


    # Plotting function for Odds Ratio S-curves
    def PlotOR(self, matrix2, matrix1, outpath, experiment_name, analysis, date):
        # Plotting sorted Odds Ratio data
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(matrix2.index, matrix2.OR, color='blue')
        ax.scatter(matrix1.index, matrix1.OR, color='red',
                    label='Significant Gene (q-val<0.1)')
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.axhline(y=1, color='black', linestyle='--')
        ax.legend()
        ax.set_ylabel('OR', fontsize=16)
        ax.set_xlabel("Index")
        ax.set_title(str("Odds Ratio-" + experiment_name + "_" + analysis))

        # Defining and adjusting the text labels for points of interest
        texts = []
        
        # Further filter matrix1 to get only the top significant genes (Cleans up the plot)
        # OR > 1
        mat1_geq1 = matrix1[matrix1.OR >= 1].sort_values(by = ['OR'], ascending = False)
        mat1_geq1 = mat1_geq1.iloc[0:5, :]
        # OR < 1
        mat1_leq1 = matrix1[matrix1.OR <= 1].sort_values(by = ['OR'], ascending = True)
        mat1_leq1 = mat1_leq1.iloc[0:5, :]
        # Concatenate them together to annotate the highest and lowest ORs
        mat1_filt = pd.concat([mat1_geq1, mat1_leq1], axis = 0)#, ignore_index=True)

        # Get the index of a gene in the significant gene matrix
        #for gene in mat1_filt.gene:
            #hold_df = mat1_filt[mat1_filt.gene == gene]
            #x, y = hold_df.index.tolist()[0], np.float64(hold_df.OR)
            #texts.append(ax.text(x, y, gene, fontsize=14))

        #adjust_text(texts, expand_text=(0.3,2), force_text = (2.2,2.2), force_points=(0.2,0.2), arrowprops=dict(arrowstyle="-", lw=1), ax=ax)

        plt.savefig(str(outpath + experiment_name + '_' +
                        analysis + '_' + str(date) + '.png'))


    # Plotting the number of alleles observed at rare variants
    def PlotEABar(self, df, gene_name, outpath, experiment_name, analysis, date):
        def SortDf(df):
            # Get sum of variants
            df.loc[:, 'Sum'] = df.Cases + \
                df.Controls

            # Dual sort for total variants and the number of alleles in cases
            df = df.sort_values(
                by=["Sum", "Cases"], ascending=True)

            # Split dataframe into variants with cases and controls and those that appear in only one group
            non0_df = df[(df.Cases >= 1) & (df.Controls >= 1)]
            zero_df = df[~df.Variants.isin(non0_df.Variants)]
            # Sort zero df by EA score
            zero_df = zero_df.sort_values(by = 'EA', ascending = True)

            # Combine them together
            df = pd.concat([zero_df, non0_df], axis = 0)
            df = df[['Variants', 'Cases', 'Controls']]

            # Get the max number of patients with variant for plotting
            y_max = df[["Cases", "Controls"]].max().max()

            return df, y_max

        # Define holding lists for gene parsing
        vars_hold = []
        cases = []
        controls = []
        ea_hold = []

        # Loop and get total alleles in cases and controls for rare variants
        for i in df.index.tolist():
            try: var = str(df.HGVSp[i].split(':')[1])
            except IndexError:
                continue
            if 'Ter' in var:
                ea = 100
                ea_hold += [float(ea)]
            else:
                ea = df.EA[i].split(',')[0]
                if ea == '.':
                    ea = np.nan
                    ea_hold += [ea]
                else:
                    ea_hold += [float(ea)]
            #ea_hold += [float(ea)]
            vars_hold += [var + "- EA(" + str(ea) + ")"]
            #vars_hold += [str(df.HGVSp[i].split(':')[1])]
            cases += [int(df.Cases[i].split(',')[2])]
            controls += [int(df.Controls[i].split(',')[2])]

        # Create and sort dataframe
        vars_df_hold = pd.DataFrame({"Variants": vars_hold, "Cases": cases, "Controls": controls, "EA": ea_hold})
        vars_df_hold, y_max = SortDf(vars_df_hold)

        # Construct plot
        fig, axes = plt.subplots(ncols=2, sharey=True)
        fig.set_size_inches(8, len(vars_hold)/3)
        # Subplot - Cases
        axes[0].barh(vars_df_hold.Variants, vars_df_hold.Cases,
                        align='center', color='black', zorder=10)
        axes[0].set(title='# Alleles - Cases')
        axes[0].set_xlim((0, y_max+10))
        for i, v in enumerate(vars_df_hold.Cases):
            axes[0].text(v + (y_max/10), i-0.25, str(v),
                            color='black', fontweight='bold')

        # Subplot - Controls
        axes[1].barh(vars_df_hold.Variants, vars_df_hold.Controls,
                        align='center', color='gray', zorder=10)
        axes[1].set(title='# Alleles - Controls')
        axes[1].set_xlim((0, y_max+10))
        for i, v in enumerate(vars_df_hold.Controls):
            axes[1].text(v + (y_max/10), i-0.25, str(v),
                            color='black', fontweight='bold')

        # Mirror cases subplot across x-axis
        axes[0].invert_xaxis()
        # Get the ticks
        tick_loc = range(len(vars_df_hold.Variants))
        # set the ticks
        axes[0].yaxis.set_ticks(tick_loc)
        axes[0].yaxis.set_ticklabels(vars_df_hold.Variants)
        axes[0].yaxis.tick_right()

        for ax in axes.flat:
            ax.margins(0.03)
            ax.grid(True)

        fig.tight_layout()
        #fig.subplots_adjust(top=0.95)
        fig.suptitle(gene_name + "-Rare Variants (AF < 0.01)", y = 0.98)
        fig.subplots_adjust(wspace=0.65)
        plt.savefig(outpath + experiment_name + analysis +
                    gene_name + "_VariantsByCount_" + str(date) + '.png')
        plt.close()


    # Calculate the FDR for observed ORs
    def ORFdr(self, ORmatrix):
        pvals = list(ORmatrix['pvalue'])
        if len(pvals) != 0:
            qvals = multipletests(
                pvals, alpha=0.05, method='fdr_bh', is_sorted=True)[1]
            #qvals = multipletests(pvals, alpha = 0.05, method='fdr_bh', is_sorted = False)[1]#bonferroni, fdr_bh
            ORmatrix['qvalue'] = qvals
            return qvals, ORmatrix


    # Calculate and plot ORs for variants grouped by gene (Only including variants with AF < 0.01)
    def ORGene(self, variants_df, outpath, analysis, date):
        Case_Allele = []
        Cont_Allele = []
        OR = []
        pval = []
        a = []
        b = []
        CI_upper = []
        CI_lower = []
        n_cases = []
        n_conts = []
        genes = list(np.unique((variants_df.gene.to_list())))
        dmatrix = pd.DataFrame(np.zeros((len(genes), 9)), columns=[
            'gene', 'n_cases', 'n_controls', 'OR', 'pvalue', 'qvalue', 'LowerCI', 'UpperCI', 'pEA'])

        for g in genes:
            a = []
            b = []
            a1 = []
            b1 = []
            A = variants_df[variants_df.gene == g]
            for i in A.index.tolist():
                # Get Allele counts in cases and controls
                a += [int(A.AC_1[i])]
                b += [int(A.AC_0[i])]
                # Get allele numbers in cases and controls (Genotyping rate)
                a1 += [int(A.AN_1[i])]
                b1 += [int(A.AN_0[i])]

            # Calculate sum of alleles across the gene
            n_cases += [np.sum(a)]
            n_conts += [np.sum(b)]
            AN_Cases = math.floor(np.mean(a1))
            AN_Conts = math.floor(np.mean(b1))
            c = AN_Cases-np.sum(a)
            d = AN_Conts-np.sum(b)
            gene_counts = [np.sum(a), np.sum(b)]
            no_gene_counts = [c, d]
            oddsratio, pvalue = stats.fisher_exact(
                [gene_counts, no_gene_counts])
            if np.sum(a) == 0 or np.sum(b) == 0 or c == 0 or d == 0:
                upper_CI = np.nan
                lower_CI = np.nan
            else:
                upper_CI = np.exp(
                    np.log(oddsratio) + 1.96*np.sqrt((1/np.sum(a)) + (1/np.sum(b)) + (1/c) + (1/d)))
                lower_CI = np.exp(
                    np.log(oddsratio) - 1.96*np.sqrt((1/np.sum(a)) + (1/np.sum(b)) + (1/c) + (1/d)))
                upper_CI = round(upper_CI, 3)
                lower_CI = round(lower_CI, 3)

            OR += [oddsratio]
            CI_upper += [upper_CI]
            CI_lower += [lower_CI]
            pval += [pvalue]

        dmatrix['gene'] = genes
        dmatrix['OR'] = OR
        dmatrix['pvalue'] = pval
        dmatrix['n_cases'] = n_cases
        dmatrix['n_controls'] = n_conts
        dmatrix['LowerCI'] = CI_lower
        dmatrix['UpperCI'] = CI_upper

        # Sort and filter out dmatrix for plotting and saving to summary table
        dmatrix1 = dmatrix.sort_values(
            by='pvalue', axis=0, ascending=True, inplace=False)
        if dmatrix1.shape[0] != 0:
            qval, dmatrix2 = self.ORFdr(dmatrix1)
        else:
            dmatrix2 = dmatrix1
        dmatrix3 = dmatrix2.sort_values(
                by='OR', axis=0, ascending=True, inplace=False).reset_index()
            
        dmatrix3.to_csv(outpath + self.expName + '_' + analysis +
                    '_Gene-Based_' + str(date) + '.csv', index=False)
        matrix1 = dmatrix3[dmatrix3['qvalue'] < 0.1] # Matrix 1 = significant genes
        matrix1 = matrix1[((matrix1['UpperCI'] > 1) & (matrix1['LowerCI'] > 1)) | ((matrix1['UpperCI'] < 1) & (matrix1['LowerCI'] < 1))]
        matrix2 = dmatrix3[dmatrix3['qvalue'] > 0.1]

        # Plotting case and control counts by gene
        for gene in matrix1.gene.unique():
            variants_df_hold = variants_df[variants_df.gene == gene].reset_index(drop = True)
            self.PlotEABar(variants_df_hold, gene, outpath,
                        self.expName, analysis, date)

        self.PlotOR(matrix2, matrix1, outpath, self.expName, analysis, date)
    

    # Calculate and plot ORs for variants grouped by gene (Only including variants with AF < 0.01)
    def ORVariant(self, variants_df, outpath, analysis, date):
        # Create holding variables for loop
        Case_Allele = []
        Cont_Allele = []
        eas = []
        afs = []
        OR = []
        pval = []
        a = []
        b = []
        CI_upper = []
        CI_lower = []
        n_cases = []
        n_conts = []

        # Get the genes of input file and create the output matrix
        genes = list((variants_df.gene.to_list()))
        hgvsp = list(variants_df.HGVSp)
        p_vars = []
        for x in hgvsp:
            try: p_vars.append(x.split(":")[1])
            except IndexError: 
                p_vars.append(x)
                continue
        #p_vars = [x.split(":")[1] for x in hgvsp]
        variants = []
        for i in range(len(hgvsp)):
            hold_gene = genes[i]
            hold_vars = p_vars[i]
            hold_variant = ":".join([hold_gene, hold_vars])
            variants.append(hold_variant)
        variants_df.loc[:, 'variants'] = variants
        dmatrix = pd.DataFrame(np.zeros((len(genes), 11)), columns=[
            'variant', 'gene', 'n_cases', 'n_controls', 'OR', 'pvalue', 'qvalue', 'LowerCI', 'UpperCI', 'EA', 'AF'])

        # Loop through each variant and calculate ORs
        for v in variants:
            #if v == 'GCSAML-AS1:p.Thr20Ala':
                #continue

            a = []
            b = []
            a1 = []
            b1 = []
            AF_hold = []
            A = variants_df[variants_df.variants == v]
            ea_hold = A.EA.values[0]
            AF_hold = A.AF.values
            eas.append(ea_hold)
            afs.append(AF_hold)
            for i in A.index.tolist():
                a += [int(A.AC_1[i])]
                b += [int(A.AC_0[i])]
                a1 += [int(A.AN_1[i])]
                b1 += [int(A.AN_0[i])]
                #AF_hold += [int(A.AF[i])]
                #print(AF_hold)

            n_cases += [np.sum(a)]
            n_conts += [np.sum(b)]
            AN_Cases = np.sum(a1)
            AN_Conts = np.sum(b1)
            c = AN_Cases-np.sum(a)
            d = AN_Conts-np.sum(b)

            gene_counts = [np.sum(a), np.sum(b)]
            no_gene_counts = [c, d]
            oddsratio, pvalue = stats.fisher_exact(
                [gene_counts, no_gene_counts])
            if AF_hold[0] >= 0.5:
                try: oddsratio = 1/oddsratio
                except ZeroDivisionError: oddsratio = oddsratio
            # Calculate CIs
            if np.sum(a) == 0 or np.sum(b) == 0 or c == 0 or d == 0:
                upper_CI = np.nan
                lower_CI = np.nan
            else:
                upper_CI = np.exp(
                    np.log(oddsratio) + 1.96*np.sqrt((1/np.sum(a)) + (1/np.sum(b)) + (1/c) + (1/d)))
                lower_CI = np.exp(
                    np.log(oddsratio) - 1.96*np.sqrt((1/np.sum(a)) + (1/np.sum(b)) + (1/c) + (1/d)))
                upper_CI = round(upper_CI, 3)
                lower_CI = round(lower_CI, 3)

            OR += [oddsratio]
            CI_upper += [upper_CI]
            CI_lower += [lower_CI]
            pval += [pvalue]

        dmatrix['variant'] = variants
        dmatrix['gene'] = variants
        dmatrix['OR'] = OR
        dmatrix['pvalue'] = pval
        dmatrix['n_cases'] = n_cases
        dmatrix['n_controls'] = n_conts
        dmatrix['LowerCI'] = CI_lower
        dmatrix['UpperCI'] = CI_upper
        dmatrix['AF'] = afs
        dmatrix['EA'] = eas

        dmatrix1 = dmatrix.sort_values(
            by='pvalue', axis=0, ascending=True, inplace=False)
        if dmatrix1.shape[0] != 0:
            qval, dmatrix2 = self.ORFdr(dmatrix1)
        else:
            dmatrix2 = dmatrix1
        dmatrix3 = dmatrix2.sort_values(
            by='OR', axis=0, ascending=True, inplace=False).reset_index(drop=True)
        dmatrix3.to_csv(outpath + self.expName + '_' + analysis +
                        '_Variant-Based' + str(date) + '.csv', index=False)
        matrix1 = dmatrix3[dmatrix3['qvalue'] < 0.1]
        matrix1 = matrix1[((matrix1['UpperCI'] > 1) & (matrix1['LowerCI'] > 1)) | (
            (matrix1['UpperCI'] < 1) & (matrix1['LowerCI'] < 1))]
        matrix2 = dmatrix3[dmatrix3['qvalue'] > 0.1]

        self.PlotOR(matrix2, matrix1, outpath, self.expName, analysis, date)


    def main(self):
        # Run ExactTest script
        nothing = self.RunExactTest()

        # Set up pool for parallelizing
        pool = mp.Pool(self.Cores)
        # Define the thresholds to analyze
        thresholds_2analyze = ['All Variants', 'EA 0-30', 'EA 30-70', 'EA 70-100']
        date = datetime.date.today()
        for method in self.methodList:
            # Create new output path
            newOROutPath = CreateDir(self.oPath, method)
            
            # Clean files
            variants_rare, variants_common = self.CleanInputs(method)

            # Check thresholds and run each method after filtering lists to thresholds
            if 'All Variants' in thresholds_2analyze:
                os.makedirs(newOROutPath + "/All Variants/", exist_ok = True)
                # Run async pooling for ORs
                pool.apply(self.ORGene, args=(variants_rare, newOROutPath+"/All Variants/", "All Variants-Rare(Gene)", date))
                pool.apply(self.ORVariant, args=(variants_common, newOROutPath+"/All Variants/", "All Variants-Common(Variant)", date))


            if 'EA 0-30' in thresholds_2analyze:
                os.makedirs(newOROutPath + "/EA 0-30/", exist_ok = True)
                # Filter variants
                variants_common_0EA30 = variants_common[(variants_common['EA-Clean'] < 30)
                                                        & (variants_common['EA-Clean'] >= 0)].reset_index(drop=True)
                variants_rare_0EA30 = variants_rare[(variants_rare['EA-Clean'] < 30)
                                                     & (variants_rare['EA-Clean'] >= 0)].reset_index(drop=True)
                # Run ORs
                pool.apply(self.ORGene, args=(
                    variants_rare_0EA30, newOROutPath+"/EA 0-30/", "EA 0-30-Rare(Gene)", date))
                pool.apply(self.ORVariant, args=(variants_common_0EA30,
                                                 newOROutPath+"/EA 0-30/", "EA 0-30-Common(Variant)", date))

            if 'EA 30-70' in thresholds_2analyze:
                os.makedirs(newOROutPath + "/EA 30-70/", exist_ok = True)
                # Filter variants
                variants_common_30EA70 = variants_common[(variants_common['EA-Clean'] >= 30) &
                                                         (variants_common['EA-Clean'] <= 70)].reset_index(drop=True)
                variants_rare_30EA70 = variants_rare[(variants_rare['EA-Clean'] >= 30) &
                                                    (variants_rare['EA-Clean'] <= 70)].reset_index(drop=True)
                # Run ORs
                pool.apply(self.ORGene, args=(
                    variants_rare_30EA70, newOROutPath+"/EA 30-70/", "EA 30-70-Rare(Gene)", date))
                pool.apply(self.ORVariant, args=(variants_common_30EA70,
                                                  newOROutPath + "/EA 30-70/", 'EA 30-70-Common(Variant)', date))

            if 'EA 70-100' in thresholds_2analyze:
                os.makedirs(newOROutPath + "/EA 70-100/", exist_ok = True)
                # Filter variants
                variants_common_70EA100 = variants_common[variants_common['EA-Clean'] > 70].reset_index(
                    drop=True)
                variants_rare_70EA100 = variants_rare[variants_rare['EA-Clean'] > 70].reset_index(
                    drop=True)

                # Run ORs
                pool.apply(self.ORGene, args=(variants_rare_70EA100,
                                               newOROutPath+"/EA 70-100/", "EA 70-100-Rare(Gene)", date))
                pool.apply(self.ORVariant, args=(variants_common_70EA100,
                                                  newOROutPath+"/EA 70-100/", "EA 70-100-Common(Variant)", date))

        pool.close()
        pool.join()

