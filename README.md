# krr-target-align
This repository contains code for simulation studies in the following manuscript:

Arash A. Amini, Richard Baumgartner, Dai Feng. "Target alignment in truncated kernel ridge regression". 2022.

## Basic info
Need following **R** packages: *Matrix, ggplot2, dplyr, tidyverse, plotly, reticulate*

## Brief explaination of all R files

[krr_funcs.R:] Creat Gaussian Kernels

[simul_pattern3_PhaseTransitonAndDoubleDecent.R:] Generate results and corresponding figures demonstratng the multiple-descent and phase transition of for $\lambda$-regularization curve

[simul_twoSegment_truncation r_mseVSr_RatioSigmaToN.R:] Generate results and corresponding figures demonstratng the double-descent and phase transition for $r$-regularization curves

[test_ta.R:] Generate results and corresponding figures demonstratng emprical vs. theroetical results

## Brief explaination of src code
[krr_calc.cpp:] C++ code for computation of MSEs
