# Network Diffusion (nDiffusion) to Validate Gene Connectedness

Repository for [*nDiffusion website*](http://ndiffusion.lichtargelab.org) 


# Publications 

- Pham M, Govindarajan H, Lichtarge O. nDiffusion: a network diffusion application to evaluate gene functionality (submitted)

- [*Pham M, Lichtarge O. Graph-based information diffusion method for prioritizing functionally related genes in protein-protein interaction networks. Pac Symp Biocomput. 2020;25:439-450*](https://www.worldscientific.com/doi/10.1142/9789811215636_0039)



If you have any questions or comments, feel free to contact Minh Pham (minh.pham@bcm.edu), Harikumar Govindarajan (Harikumar.Govindarajan@bcm.edu), or Olivier Lichtarge (lichtarge@bcm.edu).

--------
## Content
 - [Download code](#download-code)
 - [Installation and download network data](#installation-and-download-network-data)
 - [Run tutorial](#run-tutorial)

--------
## Download code
```bash
git clone https://github.com/LichtargeLab/nDiffusion.git 
```

--------
## Installation

### Install Environment
- Requirement: python=3.5.2
```bash
conda create -n nDiffusion python=3.5.2
source activate nDiffusion
pip install -r requirements.txt
```
--------
## Run tutorial

### 1. Activate environment (skip it if you already activated)
```bash
source activate nDiffusion
```
### 2. Open src/run_Diffusion.py and edit the appropriate variables in the CHANGE HERE section. For the tutorial, keep the default
#### (1) A network file: network_fl. Format: Entity1\tEntity2\tEdge_weight. Default: toy_network.txt
#### (2) Input gene lists: geneList1_fl and geneList2_fl. Format: a column of gene names
#### (3) A result folder: result_fl. Default: '../results/'

### 3. Run run_Diffusion.py
```bash
cd src/
python run_Diffusion.py
```
### 4. Check for the outputs in the result folder
=======

