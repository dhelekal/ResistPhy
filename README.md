# ResistPhy

## Introduction

ResistPhy is a R package that can estimate the cost and benefit of antimicrobial resistance from dated phylogenies and antimicrobial usage data.

For more information, see the preprint [_Estimating the fitness cost and benefit of antimicrobial resistance from pathogen genomic data_](https://www.biorxiv.org/content/10.1101/2022.12.02.518824v1).

Scripts and `.rmd` notebooks used for reproducing the data analysis and figures can be found in the directory `./run`.

## Installation

This package can be installed using 

``` r
devtools::install_github('dhelekal/ResistPhy')
```

Dependent packages available from CRAN will be automatically installed. 
However, there are three dependent packages not available on CRAN that may need manual installation beforehand:

1. [CmdStanR](https://mc-stan.org/cmdstanr/)
2. [treeio](https://bioconductor.org/packages/release/bioc/html/treeio.html)
3. [phylodyn](https://github.com/mdkarcher/phylodyn)

## Usage

The main function is `infer_costs2` which infers parameters including the fitness cost and benefit of resistance given dated phylogenies of susceptible and resistance lineages, plus the antibiotic usage function over time. Full details of how to use this and other functions can be obtained using `help(package='ResistPhy')`.

For examples of how to use ResistPhy see the applications to simulated and real datasets contained in the R markdown files in the directory `./run`. 
