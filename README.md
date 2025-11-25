# BEACON

## Overview

**BEACON** is a computational framework designed to perform Bayesian analysis on expression and gene dependency data across different cell lineages. The model integrates expression as the dependent variable and gene dependency data as the independent variable.

The framework uses **JAGS** (Just Another Gibbs Sampler) for Bayesian inference and MCMC (Markov Chain Monte Carlo) simulations to estimate parameters of interest, including the correlation between gene expression and dependency.

## Citation

Elmas A, Layden HM, Ellis JD, Bartlett LN, Zhao X, Kawabata-Iwakawa R, Obinata H, Hiebert SW, Huang KL. Expression-Driven Genetic Dependency Reveals Targets for Precision Medicine. bioRxiv [Preprint]. 2024 Oct 21:2024.10.17.618926. doi: 10.1101/2024.10.17.618926. PMID: 39484404; PMCID: PMC11527036.
[![DOI](https://zenodo.org/badge/DOI/doi.org/10.1101/2024.10.17.618926.svg)](https://doi.org/10.1101/2024.10.17.618926)

## Features

- **Flexible Data Input:** Supports multiple data types such as mRNA, Protein, and RNA transcripts.
- **Customizable Parameters:** Users can adjust parameters like number of iterations, adaptation steps, and lineages of interest.
- **Reproducibility:** The code can reproduce results or calculate false discovery rates (FDR) based on user inputs; Generates detailed output files with Bayesian analysis results for each lineage.

## Installation

1. Clone the repository:

2. Install the necessary R packages:
   ```R
   install.packages(c("openxlsx", "rjags"))
   ```
   
3. Install JAGS:
   - JAGS can be downloaded and installed from [JAGS official site](https://sourceforge.net/projects/mcmc-jags/).

4. Download required files and setup folders:

   - sample_info.csv     @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=35020903
   - CCLE_expression.csv (further gzipped)     @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34989919
   - Supplementary data (Table S2: normalized protein expressions) from Nusinow et al. paper (doi.org/10.1016/j.cell.2019.12.023)     @ https://www.cell.com/cms/10.1016/j.cell.2019.12.023/attachment/3709dedc-3a01-4e1d-ab4c-82597295c5d2
   - CRISPR_gene_effect.csv (further gzipped)     @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34990036

   - Folder structure:  
   ```
   BEACON-main/
   ├── LineageMCMC.R
   ├── PanLineageMCMC.R
   ├── DepMap_data/
   │   ├── sample_info_22Q2.csv  
   │   ├── CCLE_expression_22Q2.csv.gz
   │   └── CRISPR_gene_effect_22Q2.csv.gz
   ├── QuantProtCCLE_Nusinow_Cell2020/  
   │   └── mmc2.xlsx  
   └── out/
   ```

## Computational Requirements

**Runtime:** Approximately 20 hours on a 16-core processor with 128 GB memory. Actual runtime may vary depending on system specifications and the number of genes analyzed. For subset analyses or testing, consider running on a limited set of genes initially.

**System Requirements:**
- R version 4.4.3 or later
- See `requirements.txt` for complete package versions and dependencies

## Usage

To run the analysis, modify the R scriptS according to your data and parameters. The primary script performs the following steps:

1. **Data Preparation:**
   - Load mRNA/protein expression and CRISPR dependency data.
   - Compress the expression and the dependency data (gzip CCLE_expression.csv | gzip CRISPR_gene_effect.csv).
   - Map and filter data based on lineage and gene selection.
   
2. **Model Initialization:**
   - Initialize the Bayesian model with uninformative priors.
   
3. **Run MCMC:**
   - Perform MCMC simulations for each lineage.
   - Save results to an output directory.
   
4. **Reproducibility:**
   - Optionally reproduce previous results by loading existing data and recalculating FDR.
   
### Example

```R
# Example of running the analysis
n.adapt = 100
n.update = 100
n.iter = 500
reproduce.results = TRUE

# Load and prepare mRNA data
data = 'mRNA'
cell.type = 'All'
panel = ''

# Run the Bayesian analysis for a specific lineage
lineage='SOFT.TISSUE'
# Modify other parameters as necessary
# ...

# Run the analysis
source('LineageMCMC.R')
```

## Output

The analysis generates the following output files:
- `Table.<data>.dependency.Bayesian.lineage.<lineage>.<panel>.xlsx`: Summary of Bayesian analysis for each lineage.
- Log files and intermediate results saved in the specified output directory.
