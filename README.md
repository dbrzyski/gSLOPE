# gSLOPE

## group SLOPE &mdash; a statistical method for the selection of explanatory variables with gFDR control.

This repository contains code for the simulation studies and figures presented in Brzyski, D., Gossmann, A., Su, W., & Bogdan, M. (2016). *Group SLOPE-adaptive selection of groups of predictors*, [arXiv:1610.04960](http://arxiv.org/abs/1610.04960).

### Matlab version

gSLOPE uses the C code for the fast prox evaluation and the interface to this C code must be compiled prior to use. In Matlab this is done by running the makemex script, available in "\gSLOPE_matlab\SLOPE_code", from the Matlab command line.

To start using the software we suggest to run the file [`WorkingExample.m`](https://github.com/dbrzyski/gSLOPE/blob/master/WorkingExample.m).

### R version

The subdirectory [`R/grpSLOPE`](https://github.com/dbrzyski/gSLOPE/tree/master/R/grpSLOPE) contains the R versions of the simulations. The code uses the R package `grpSLOPE` ([on CRAN](https://CRAN.R-project.org/package=grpSLOPE), [on github](https://github.com/agisga/grpSLOPE)).

### gSLOPE paper: figures and data

According to dbGaP policy we could not provide any part of Human Genetic Data, which were used do get Figure 5 and Figure 6

All scripts needed to generate Figure 1, Figure 2, Figure 3 and Figure 4 are available. The user needs to update the paths inside scripts to the data and to the directory where the results should be saved.

Figures 2-6 were generated with the usage of ggplot2 R package. We used knitr - dynamic report generation with R - to get pdf plots. These plots are stored in corresponding 'Figures#_files', available in 'Figures\Rscripts_for_figures\'
