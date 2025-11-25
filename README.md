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

## Usage

To run the analysis, modify the R scriptS according to your data and parameters. The primary script performs the following steps:

1. **Data Preparation:**
   - Load mRNA/protein expression and CRISPR dependency data.
   - Gzip the expression and the dependency data.
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
