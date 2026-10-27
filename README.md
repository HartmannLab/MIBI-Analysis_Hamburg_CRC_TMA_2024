# Robust multicellular programs dissect the complex tumor microenvironment and track disease progression in colorectal cancer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17008987.svg)](https://zenodo.org/badge/latestdoi/17008987)

**Authors:** Loan Vulliard, Teresa Glauner, Sven Truxa, Miray Cetin, Yu-Le Wu, Ronald Simon, Laura Behm, Jovan Tanevski, Julio Saez-Rodriguez, Guido Sauter, Felix J. Hartmann

This repository contains the computational analysis pipeline for our study on multicellular programs in colorectal adenocarcinomas using MIBI (Multiplexed Ion Beam Imaging) technology on tissue microarrays (TMAs) from CRC biopsies collected at the Institute of Pathology, University Medical Center Hamburg-Eppendorf. You can find the full preprint here: https://arxiv.org/abs/2510.05083

## Overview

This project employs spatial proteomics to characterize the tumor microenvironment in colorectal cancer, identifying robust multicellular programs that track disease progression. The analysis pipeline includes preprocessing, cell type annotation, multicellular program identification, and disease stage prediction using machine learning approaches.

## Repository Structure

```
├── notebooks/                   # Jupyter notebooks for analysis
│   ├── .utils/                  # Helpers and scripts shared by notebooks
│   ├── preprocessing/           # Data preprocessing and filtering
│   ├── lineage/                 # Cell type annotation workflows
│   ├── metabolism/              # Metabolic marker analysis
│   ├── metadata/                # Clinical metadata processing
│   ├── multicellular/           # Multicellular program analysis
│   ├── qc/                      # Data quality checks
│   ├── spatial/                 # Spatial analyses
│   ├── validation/              # Transcriptomics reanalyses
│   └── viz/                     # Visualization notebooks
├── data/                        # Data files (not included in repo, see DATA.md)
├── environments/                # Environment specification files, see environments/README.md
├── DATA.md                      # Description of the shared data files
└── README.md                    # This file
```

## Data

The cell-level measurements, cell type annotations and clinical metadata are deposited on Zenodo: [10.5281/zenodo.18504904](https://doi.org/10.5281/zenodo.18504904). The supplementary tables are deposited separately: [10.5281/zenodo.17280877](https://doi.org/10.5281/zenodo.17280877). Multiplexed images and segmentation masks will be provided on the BioImage Archive.

[DATA.md](DATA.md) describes the columns of each file, how the files connect, and how to load them.

## Computational Environments

The analysis steps use six environments:

| Environment | Used for | Specification |
|---|---|---|
| Docker image `koalive/mibiprofiling` | Preprocessing and notebooks requiring Julia packages | [DockerHub](https://hub.docker.com/r/koalive/mibiprofiling) |
| Kasumi | R-based statistics and visualization, MISTy and Kasumi | `environments/kasumi.yml`, `environments/kasumi.csv` |
| scverse | Factor analysis and enrichment analyses | `environments/scverse.yml` |
| xgboost | Cell type annotation, machine learning predictions and PyTorch-based analyses | `environments/xgboost.yml` |
| scverse-squidpy | Factor refits with other spatial encodings, including squidpy | `environments/scverse_squidpy.yml` |
| scimap_cozi | COZI enrichment between cell types | `environments/scimap_cozi.yml` |

[environments/README.md](environments/README.md) lists the notebooks that use each environment and explains how to create it. These are not minimal environments, but the original ones used to design and run the analyses, to reproduce the manuscript results accurately.

## Getting Started

1. **Set up environments**: follow [environments/README.md](environments/README.md).
2. **Data**: download the files from Zenodo into `data/`, as described in [DATA.md](DATA.md).
3. **Preprocessing**: start with notebooks in `notebooks/preprocessing/` using the Docker environment.
4. **Cell annotation**: proceed to `notebooks/lineage/` for cell type identification.
5. **Metadata preparation**: harmonize patient annotations for stratified analyses with notebooks in `notebooks/metadata/`.
6. **Multicellular analysis**: use notebooks in `notebooks/multicellular/` for higher-order analysis.
7. **Visualization**: generate additional figures using notebooks in `notebooks/viz/`.

The deposited files are the outputs of steps 3 and 4 and the input of step 5, so analyses can also start from them directly.

## MuVIcell

The multicellular analysis we performed can be adapted to your use cases using the dedicated [MuVIcell package](https://github.com/HartmannLab/MuVIcell).

## Contact

For questions about the analysis pipeline or code, please open an issue in this repository.

## Citation

If you use this code or analysis pipeline, please cite our manuscript: https://arxiv.org/abs/2510.05083


