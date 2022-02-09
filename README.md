# CriteriaForSuccess
The integration of Lichtarge Lab genotype-phenotype validation experiments

Installation
1. git clone https://github.com/kevwilhelm95/CriteriaForSuccess.git
2. conda env create -f environment.yml
3. conda activate pyCFS


Required arguments:

| Argument | Description |
|--------- | ----------- |
|--ExperimentName | Name of disease and or cohort |
|--InputPath | Path to BigPipeline Results |
|--GSPath | Path to CSV of Gold Standard Lists |
|--CaseControlPath | Path to CSV; Col. 1- IDs, Col.2- 1,0 (Case,Control); No header |
|--ExactTestPath | Path to .txt output from ExactTest.sh |
|--OutPutPath | Path to output directory |


Optional arguments:

| Argument | Description |
|----------|-------------|
|--nDiffusionGraph | Network to use for nDiffusion; Choices = "STRINGv11", "STRINGv10", "MeTeOR" |
|--AC_Threshold | Threshold of Pathways to incorporate into Consensus2 (Default = 5) |
|--cores | Number of cores used to run |



Required reference files not in repo:
- STRINGv11.txt
- STRINGv10.txt
- MeTeOR.txt
