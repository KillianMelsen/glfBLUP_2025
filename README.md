# gfBLUP_analyses

Contains all scripts required to generate the simulated data and analyze it using different methods.

## **Part 1:** data simulation
Two scripts are used to generate that data. The first one simulates the actual data based on the factor model while the second one randomly divides each dataset into a training and test set. Note that these scripts can be run on WSL using Intel's oneMKL BLAS and LAPACK libraries. In that case MKL_DYNAMIC and MKL_NUM_THREADS must be set to FALSE and 3, respectively. This assumes a 20 thread system (5 parallel processes each using 3 threads for MKL). For an 8 thread system MKL_NUM_THREADS should be set to 1.

## **Part 2:** data analyses
Placeholder

## **Part 3:** result merging
Placeholder

## **Part 4:** result visualization
Placeholder
