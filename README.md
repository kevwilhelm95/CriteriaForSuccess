# CriteriaForSuccess
The aggregation of Lichtarge Lab genotype-phenotype validation experiments

Installation
1. git clone https://github.com/kevwilhelm95/CriteriaForSuccess.git
2. conda env create -f environment.yml
3. conda activate pyCFS

# Required arguments:

| Argument | Description | Options |
|--------- | ----------- | ------- |
|--ExperimentName | Name of disease and/or cohort | NA|
|--InputPath or --InputList | Path to BigPipeline result directory or single-column .txt file with genes of interest | NA |
|--PickExperiments| No-space comma-separated list of experiments to run | All (default), GS Overlap, nDiffusion, MGI, OR, PubMed Enrichment, Variants By Sample, Intermethod Connectivity, Pharmacology |
|--OutPutPath | Path to output directory | NA |
|--cores | Number of cores used to run | NA | 

# Optional arguments:

| Argument | Description | Options | Required for |
|----------|-------------|---------|--------------|
|--GSPath | Path to CSV of Gold Standard Lists | NA | All, GS Overlap, nDiffusion, Intermethod Connectivity |
|--VCF | Path to VCF | NA | All, OR |
|--PickNetwork | Network to use for nDiffusion/Intermethod Connectivity | STRINGv10, STRINGv11, MeTEOR (default: STRINGv11) | All, nDiffusion, Intermethod Connectivity |
|--PubMedKeywords | No-space comma-separated list of key words to query for co-mention (e.g. "Type 2 Diabetes,Insulin") | NA | All, PubMed Enrichment |
|--InterConnectivity_Evidences | No-space comma-separated list of evidence types | neighborhood, fusion, cooccurence, coexpression, experimental, database, textmining | All, Intermethod Connectivity |
|--InterConnectivity_Confidence | Edge weight to test for connections | all, medium, high, highest | All, Intermethod Connectivity |
|--CaseControlPath | Path to CSV; Col. 1- IDs, Col.2- 1,0 (Case,Control); No header | NA | All, OR, Variants By Sample |
|--AC_Threshold | Threshold of Pathways to incorporate into Consensus2 (Default = 5) | NA | Only used if --InputPath used |
|--ref | Genomic reference locations | GRCh37, GRCh38 (default), hg19, hg38 | All, OR, Variants By Sample |

# Example Call
python CriteriaForSuccess.py \
--ExperimentName Test_run \
--InputPath /path/to/dir/ \
--PickExperiments All \
--CaseControlPath /path/to/sample.csv \
--OutPutPath /path/to/output/dir \
--cores 1 \
--GSPath /path/to/GS.csv \
--VCF /path/to/vcf.gz \
--PickNetwork STRINGv11 \
--PubMedKeywords "Type 2 Diabetes,Insulin" \
--InterConnectivity_Evidences "neighborhood,fusion,coexpression" \
--InterConnectivity_Confidence highest \
--CaseControlPath /path/to/samples.csv \
--AC_Threshold 5 \
--ref GRCh38

# Experiment Descriptions
| Experiment | Description |
|-----------| -------------|
| GS Overlap | Number of Gold Standard genes recovered by experimental list. Significance tested by hypergeometric test |
| nDiffusion | Network connectivity against Gold Standard genes using nDiffusion algorithm (PMID: 31797167). Significance tested by Z score against 100 random, degree-matched gene sets |
| InterMethod Connectivity | First-neighbor enrichment to Gold Standard genes. Significance tested by Z score against 100 random, degree-matched gene sets |
| PubMed Enrichment | Co-mentions in title or abstract by query gene for defined keywords. Significance tested by Z score using 100 random gene sets |
| OR (Odds Ratio) | Allele frequency differences between cases/controls for by variants (Allele frequency > 1%) or by gene (Allele frequency < 1% variants collapsed by gene). Significance tested by Fisher's Exact Test |
| Variants By Sample | Parser to obtain .txt file containing information on which samples harbor which variants |
| MGI | Enrichment test for high-level phenotypic abnormalities observed in mouse-models obtained from the Mouse Genome Informatics resource. Significance tested by Fisher's Exact Test |
| Pharmacology | Pulls drug-protein interactions for query genes from DGIdb.org |