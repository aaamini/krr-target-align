# krr-target-align
This repository contains code for simulation studies in the following manuscript:

Arash A. Amini, Richard Baumgartner, Dai Feng. "Target alignment in truncated kernel ridge regression". 2022.

## Basic info
Need following **R** packages: *Matrix, ggplot2, dplyr, tidyverse, plotly, reticulate*

## Brief explanation of all R files

[krr_funcs.R:] Create Gaussian kernels

[simul_pattern3_PhaseTransitonAndDoubleDecent.R:] Generate results and corresponding figures demonstrating the multiple-descent and phase transition of $\lambda$-regularization curve

[simul_twoSegment_truncation r_mseVSr_RatioSigmaToN.R:] Generate results and corresponding figures demonstrating the double-descent and phase transition of $r$-regularization curves

[test_ta.R:] Generate results and corresponding figures demonstrating empirical vs. theoretical results

## Brief explanation of src code
[krr_calc.cpp:] C++ code for computation of MSEs
