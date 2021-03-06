# ICcalibReproduce
Reproducibility of results from interval-censoring calibration paper:

A novel calibration framework for survival analysis when a binary covariate is measured at sparse time points, 
Nevo et al., Biostatistics, in press
https://doi.org/10.1093/biostatistics/kxy063

R code and simulation results to reproduce results in ICcalib paper. This repository includes allows to fully reproduce Table 2, and Tables A.1-A.3. The two main folder are

1. DataAnalysisCMV: Implementing our method via the ICcalib R package on the CMV data that was previously analyzed by Goggins et al. (1999). 
2. Simulations: Fully reproduce all simulation results from the paper. Auxiliary functions that carry the simulations and simulation scripts, values for `set.seed` and actual results obtained from simulations in the form of .txt files
