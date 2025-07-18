# Improving Genomic Prediction using High-dimensional Secondary Phenotypes: The Genetic Latent Factor Approach

This repository contains all scripts required to generate the simulated, hyperspectral, and timing data. It also contains scripts to analyze the datasets using different methods and visualize the results, as in the paper. The individual scripts can all be run conveniently from four main files that source all scripts, [`run_all_p800.R`](run_all_p800.R), [`run_all_hyper.R`](run_all_hyper.R), [`run_all_timing.R`](run_all_timing.R), and [`run_all_misc.R`](run_all_misc.R).

Source all scripts in the order presented in these main files. No manual changes to any scripts are required.

## Notes:

-   All analyses together take a significant amount of time. Things can be sped up considerably if run in Windows subsystem for Linux (WSL) using Intel's oneMKL BLAS and LAPACK libraries. If using oneMKL, set MKL_DYNAMIC to FALSE and MKL_NUM_THREADS to 1 to avoid issues with parallelization.
-   There are five scripts that allow for checking of reproducibility of some intermediate results. These can be sourced from [`run_all_rep_checks.R`](run_all_rep_checks.R).
-   Deep learning was run on the CPU in Windows.
-   It is safest to run MegaLMM on Linux or in WSL as I've had it crash randomly in Windows.
-   Info on the hyperspectral dataset can be found in the paper by [Krause et al. 2019](https://doi.org/10.1534/g3.118.200856).
-   The [glfBLUP package](glfBLUP_1.0.0.tar.gz) is included in this repository, but can be installed from [this github repository](https://github.com/KillianMelsen/glfBLUP) as well.
-   There is an [old version of the glfBLUP package](gfBLUPold_1.3.1.tar.gz) included in this repository that is required to run the lsBLUP and siBLUP models.

## Folders:
-   [genotypes](genotypes) contains marker data and kinships for the analyses.
-   [helper_functions](helper_functions) contains functions required for the analyses and data generation.
-   [hyper_1415B5IR](hyper_1415B5IR) contains all files for the hyperspectral B5IR analyses.
-   [hyper_1415HEAT](hyper_1415HEAT) contains all files for the hyperspectral HEAT analyses.
-   [hyper_datafiles](hyper_datafiles) contains the original hyperspectral data.
-   [logs](logs) would contain log files from analyses that make use of parallelization.
-   [p800](p800) contains all files for the simulations.
-   [plots](plots) contains all plots.
-   [rep_check_scripts](rep_check_scripts) contains the scripts that check reproducibility of some intermediate results.
-   [timing](timing) contains files that measure runtimes of glfBLUP for different numbers of features.
-   [gfBLUPold_1.3.1.tar.gz](gfBLUPold_1.3.1.tar.gz) is an old version of the glfBLUP package required to run lsBLUP and siBLUP.
-   [glfBLUP_1.0.0.tar.gz](glfBLUP_1.0.0.tar.gz) is the version of the glfBLUP package used in the glfBLUP analyses.
-   [plot_hyper.R](plot_hyper.R) is the script that generates the combined plot of the B5IR and HEAT hyperspectral results.
-   [run_all_hyper.R](run_all_hyper.R) is the masterscript for all hyperspectral analyses.
-   [run_all_misc.R](run_all_misc.R) is the masterscript for some additional plotting.
-   [run_all_p800.R](run_all_p800.R) is the masterscript for all simulations.
-   [run_all_p800_v2.R](run_all_p800_V2.R) is the masterscript for all simulations using a low-rank residual covariance structure.
-   [run_all_rep_checks.R](run_all_rep_checks.R) is the masterscript for reproducibility checks.
-   [run_all_timing.R](run_all_timing.R) is the masterscript for the glfBLUP timing.

## Figures:
-   Figure 2 corresponds to [p800_main.png](plots/p800_main.png).
-   Figure 3 corresponds to [hyper.png](plots/hyper.png).
-   Figure 4 corresponds to [gfBLUP_hyper_1415B5IR_single_date.png](plots/gfBLUP_hyper_1415B5IR_single_date.png).
-   Figure 5 corresponds to [p800.png](plots/p800.png).
-   Figure S1 corresponds to [siBLUP_hyper_1415B5IR_single_date.png](plots/siBLUP_hyper_1415B5IR_single_date.png) and [lsBLUP_hyper_1415B5IR_single_date.png](plots/lsBLUP_hyper_1415B5IR_single_date.png).
-   Figure S2 corresponds to [MegaLMM_hyper_1415B5IR_single_date_M3.png](plots/MegaLMM_hyper_1415B5IR_single_date_M3.png), [MegaLMM_hyper_1415B5IR_single_date_M5.png](plots/MegaLMM_hyper_1415B5IR_single_date_M5.png), and [MegaLMM_hyper_1415B5IR_single_date_M10.png](plots/MegaLMM_hyper_1415B5IR_single_date_M10.png).
-   Figure S3 corresponds to [MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M3.png](plots/MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M3.png).
-   Figure S4 corresponds to [MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M5.png](plots/MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M5.png).
-   Figure S5 corresponds to [MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M10_part1.png](plots/MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M10_part1.png).
-   Figure S6 corresponds to [MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M10_part2.png](plots/MegaLMM_hyper_1415B5IR_single_date_traceplot_WL_M10_part2.png).
-   Figure S7 corresponds to [MegaLMM_hyper_1415B5IR_single_date_traceplot_Y_M3.png](plots/MegaLMM_hyper_1415B5IR_single_date_traceplot_Y_M3.png), [MegaLMM_hyper_1415B5IR_single_date_traceplot_Y_M5.png](plots/MegaLMM_hyper_1415B5IR_single_date_traceplot_Y_M5.png), and [MegaLMM_hyper_1415B5IR_single_date_traceplot_Y_M10.png](plots/MegaLMM_hyper_1415B5IR_single_date_traceplot_Y_M10.png).
-   Figure S8 corresponds to [timing.png](plots/timing.png).
-   Figure S9 corresponds to [hyper_RF_1415B5IR.png](plots/hyper_RF_1415B5IR.png).
-   Figure S10 corresponds to [hyper_RF_1415HEAT.png](plots/hyper_RF_1415HEAT.png).

## sessionInfo() output:
```R
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_Netherlands.utf8  LC_CTYPE=English_Netherlands.utf8    LC_MONETARY=English_Netherlands.utf8 LC_NUMERIC=C                        
[5] LC_TIME=English_Netherlands.utf8    

time zone: Europe/Amsterdam
tzcode source: internal

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggforce_0.4.2      statgenGWAS_1.0.9  statgenHTP_1.0.6.1 LMMsolver_1.0.8    lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1     
 [8] dplyr_1.1.4        purrr_1.0.2        readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       tidyverse_2.0.0    ggplot2_3.5.1     
[15] rrBLUP_4.6.3       tensorflow_2.16.0  keras_2.15.0       MegaLMM_0.9.5      MCMCglmm_2.36      ape_5.8-1          coda_0.19-4.1     
[22] Matrix_1.7-0       glfBLUP_1.0.0      gfBLUPold_1.3.1    doParallel_1.0.17  iterators_1.0.14   foreach_1.5.2      tictoc_1.2.1      
[29] rlist_0.4.6.2     

loaded via a namespace (and not attached):
 [1] dotCall64_1.1-1   spam_2.10-0       gtable_0.3.5      tensorA_0.36.2.1  lattice_0.22-6    tzdb_0.4.0        vctrs_0.6.5       tools_4.4.1      
 [9] tfruns_1.5.3      generics_0.1.3    fansi_1.0.6       SpATS_1.0-18      pkgconfig_2.0.3   data.table_1.16.0 lifecycle_1.0.4   cubature_2.1.1   
[17] farver_2.1.2      compiler_4.4.1    munsell_0.5.1     codetools_0.2-20  pillar_1.9.0      nloptr_2.1.1      whisker_0.4.1     MASS_7.3-61      
[25] boot_1.3-31       abind_1.4-8       nlme_3.1-166      tidyselect_1.2.1  digest_0.6.37     stringi_1.8.4     reshape2_1.4.4    splines_4.4.1    
[33] polyclip_1.10-7   rprojroot_2.0.4   here_1.0.1        colorspace_2.1-1  cli_3.6.3         magrittr_2.0.3    base64enc_0.1-3   utf8_1.2.4       
[41] corpcor_1.6.10    withr_3.0.1       scales_1.3.0      rappdirs_0.3.3    timechange_0.3.0  lme4_1.1-35.5     reticulate_1.39.0 hms_1.1.3        
[49] png_0.1-8         rlang_1.1.4       Rcpp_1.0.13       zeallot_0.1.0     glue_1.7.0        tweenr_2.0.3      rstudioapi_0.16.0 minqa_1.2.8      
[57] jsonlite_1.8.9    R6_2.5.1          plyr_1.8.9
```
