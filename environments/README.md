# Analyses environments

### 1. Docker Environment: `koalive/mibiprofiling`
**Used for:** Preprocessing steps and notebooks requiring Julia packages.  
Available from [DockerHub](https://hub.docker.com/r/koalive/mibiprofiling).

### 2. Kasumi R Environment
**Used for:** R-based statistical analysis and visualization.  
Create conda environment from `kasumi.yml` and install R packages listed in `kasumi.csv`.

### 3. scverse Environment
**Used for:** Factor analysis and enrichment analyses.  
Create conda environment from `scverse.yml`.

### 4. xgboost Environment
**Used for:** Machine learning predictions and PyTorch-based analyses.  
Create conda environment from `xgboost.yml`. Optimized for NVIDIA GPU and CUDA 12.3.  

*Note: to ensure accurate reproduction of the manuscript results, these are not minimal environments but the original ones used when designing and running the analyses.*