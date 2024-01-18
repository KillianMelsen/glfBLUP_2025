# gfBLUP 2024

Contains all scripts required to generate the simulated, hyperspectral, and timing data. Also contains scripts to analyze the datasets using different methods and visualize the results.

## 1 - High-dimensional p800 simulated data

This section provides information on the generation and analysis of the high-dimensional simulated data.

### 1.1 - Data simulation

Two scripts are used to generate that data. The first one simulates the actual data based on the factor model while the second one randomly divides each dataset into a training and test set. Note that these scripts can be run on WSL using Intel's oneMKL BLAS and LAPACK libraries. In that case MKL_DYNAMIC and MKL_NUM_THREADS must be set to FALSE and 3, respectively. This assumes a 20 thread system (5 parallel processes each using 3 threads for MKL). For an 8 thread system MKL_NUM_THREADS should be set to 1.

### 1.2 - Data analysis

Placeholder

### 1.3 - Result visualization

Placeholder

## 2 - Hyperspectral data

This section provides information on the pre-processing and analysis of the CIMMYT hyperspectral data.

### 2.1 - Data pre-processing

Placeholder

### 2.2 - Data analysis

Placeholder

### 2.3 - Result visualization

Placeholder

### 2.4 - Interpretation of a single date

Placeholder

## 3 - Timing data

This section provides information on the generation of the data used to visualize the computational times of gfBLUP as p grows.

### 3.1 - Data simulation

Placeholder

### 3.2 - Timing

Placeholder

### 3.3 - Result visualization

Placeholder

# Supplementary Material

Some extra material related to MegaLMM.

## S1 - MegaLMM

Some basic results related to a few curious characteristics of MegaLMM with regards to the number of factors and the posterior factor loading and yield prediction samples.

### S1.1 - Example data

Placeholder

### S1.2 - Sampling

Placeholder

### S1.3 - Traceplotting

Placeholder
