# Improving Genomic Prediction using High-dimensional Secondary Phenotypes: The Genetic Latent Factor Approach

This repository contains all scripts required to generate the simulated, hyperspectral, and timing data. It also contains scripts to analyze the datasets using different methods and visualize the results, as in the paper. The individual scripts can all be run conveniently from four main files that source all scripts, [`run_all_p800.R`](run_all_p800.R), [`run_all_hyper.R`](run_all_hyper.R), [`run_all_timing.R`](run_all_timing.R), and [`run_all_misc.R`](run_all_misc.R).

## Notes:

-   All analyses together take a significant amount of time. Things can be sped up considerably if run in Windows subsystem for Linux (WSL) using Intel's oneMKL BLAS and LAPACK libraries. If using oneMKL, set MKL_DYNAMIC to FALSE and MKL_NUM_THREADS to 1 to avoid issues with parallelization.
-   Deep learning was run on the CPU in Windows.
-   It is safest to run MegaLMM on Linux or in WSL as I've had it crash randomly in Windows.
