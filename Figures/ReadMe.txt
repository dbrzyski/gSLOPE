gSLOPE uses the C code for the fast prox evaluation and the interface to this C code must be compiled prior to use. In Matlab this is done by running the makemex script, available in "\gSLOPE_matlab\SLOPE_code", from the Matlab command line.

According to dbGaP policy we could not provide any part of Human Genetic Data, which were used do get Figure 5 and Figure 6

All scripts needed to generate Figure 1, Figure 2, Figure 3 and Figure 4 are available. The user needs to update the paths inside scripts to the data and to the directory where the results should be saved.

Figures 2-6 were generated with the usage of ggplot2 R package. We used knitr - dynamic report generation with R - to get pdf plot. These plots are stored in corresponding 'Figures#_files', available in 'Figures\Rscripts_for_figures\'