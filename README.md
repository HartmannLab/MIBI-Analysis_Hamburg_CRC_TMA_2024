# Robust multicellular programs dissect the complex tumor microenvironment in colorectal adenocarcinomas and track disease progression

[![DOI](https://zenodo.org/badge/750937733.svg)](https://doi.org/10.5281/zenodo.17008987)

**Authors:** Loan Vulliard, Teresa Glauner, Sven Truxa, Miray Cetin, Yu-Le Wu, Laura Behm, Jovan Tanevski, Julio Saez-Rodriguez, Guido Sauter, Felix J. Hartmann

This repository contains the computational analysis pipeline for our study on multicellular programs in colorectal adenocarcinomas using MIBI (Multiplexed Ion Beam Imaging) technology on tissue microarrays (TMAs) from CRC biopsies collected at the Institute of Pathology, University Medical Center Hamburg-Eppendorf.

## Overview

This project employs spatial proteomics to characterize the tumor microenvironment in colorectal cancer, identifying robust multicellular programs that track disease progression. The analysis pipeline includes preprocessing, cell type annotation, multicellular program identification, and disease stage prediction using machine learning approaches.

## Repository Structure

```
├── notebooks/                   # Jupyter notebooks for analysis
│   ├── preprocessing/           # Data preprocessing and filtering
│   ├── lineage/                 # Cell type annotation workflows
│   ├── metabolism/              # Metabolic marker analysis
│   ├── metadata/                # Clinical metadata processing
│   ├── multicellular/           # Multicellular program analysis
│   └── qc/                      # Data quality checks
│   └── spatial/                 # Spatial analyses
│   └── validation/              # Transcriptomics reanalyses
│   └── viz/                     # Visualization notebooks
├── mibi_analysis/               # Python package for MIBI analysis
│   ├── __init__.py             # Package initialization
│   ├── preprocessing.py        # Data normalization and preparation
│   ├── mofa.py                 # MuVI/MOFA wrapper functions
│   ├── visualization.py       # Plotting and visualization
│   ├── interpretation.py      # Statistical analysis and interpretation
│   └── README.md               # Package documentation
├── examples/                    # Usage examples
│   └── mibi_analysis_example.ipynb # Complete workflow example
├── data/                        # Data files (not included in repo)
├── environments/                # Environment specification files
├── setup.py                     # Package installation script
└── README.md                    # This file
```

## MIBI Analysis Package

This repository now includes a Python package (`mibi_analysis/`) that provides helper functions for preprocessing, running MuVI factor analysis, and visualizing results from MIBI data. The package extracts and modularizes the key functionality from the `notebooks/multicellular/MOFACell.ipynb` notebook.

### Package Features

- **Preprocessing**: Normalize muon objects with `normalize_mudata_features()` and `prepare_mudata_for_mofa()`
- **MuVI Integration**: Run factor analysis with `run_mofa_analysis()` and extract results with `extract_factor_scores()`  
- **Visualization**: Generate plots with `plot_factor_scores()`, `plot_factor_loadings()`, and confidence ellipses
- **Interpretation**: Test clinical associations with `test_factor_associations()` and identify top loadings

### Quick Start

```python
import mibi_analysis as ma
import mudata as mu

# Load and preprocess data
mdata = mu.read_h5mu("data/celltype_features.h5mu")
ma.prepare_mudata_for_mofa(mdata, normalize=True)

# Run MOFA analysis  
model_path = ma.run_mofa_analysis(mdata, n_factors=10)

# Extract and visualize results
factor_scores = ma.extract_factor_scores(model_path, mdata=mdata)
plot = ma.plot_factor_scores(factor_scores, 'Factor 1', 'Factor 2', color_by='Stage')
```

See `examples/mibi_analysis_example.ipynb` for a complete usage example.

## Computational Environment Requirements

Due to the diverse computational requirements of different analysis steps, this project uses **four distinct environments**:

### 1. Docker Environment: `koalive/mibiprofiling`
**Used for:** Preprocessing steps and notebooks requiring Julia packages
- **Notebooks:**
  - `notebooks/preprocessing/FilterAndViz_Hamburg.ipynb` - Data filtering using BioProfiling.jl
  - `notebooks/preprocessing/FilterSegmentationMask.ipynb` - Segmentation mask filtering
  - `notebooks/viz/ExampleImages.ipynb` - Pipeline visualization examples

### 2. Kasumi R Environment
**Used for:** R-based statistical analysis and visualization
- **Notebooks:**
  - `notebooks/metadata/CohortPlotR.ipynb` - Clinical cohort visualization and statistics
  - `notebooks/multicellular/BiplotFactors.ipynb` - Visualize factors associated with clinical variables
  - `notebooks/spatial/KasumiColorectal.ipynb` - Spatial metabolic analysis with Kasumi
  - `notebooks/spatial/MistyLines.ipynb` – Spatial lineage analysis with MISTy

### 3. scverse Environment
**Used for:** Factor analysis and enrichment analyses
- **Notebooks:**
  - `notebooks/multicellular/PrepareMuData.ipynb` - Multi-modal data preparation
  - `notebooks/multicellular/MOFACell.ipynb` - Factor analysis

### 4. xgboost Environment
**Used for:** Machine learning predictions and PyTorch-based analyses
- **Notebooks:**
  - `notebooks/lineage/annotation_scyan.ipynb` - Cell type annotation using Scyan
  - `notebooks/lineage/annotation_consensus.ipynb` - Consensus cell type annotation
  - `notebooks/lineage/Annotation_gating_Hamburg.ipynb` - Manual gating-based annotation
  - All notebooks involving XGBoost predictions and deep learning models

*Note: to ensure accurate reproduction of the manuscript results, these are not minimal environments but the original ones used when designing and running the analyses.*

## Environment Setup

Detailed environment specification files and installation instructions can be found in the `environments/` directory. Each environment includes:
- Package lists and version specifications
- Installation commands
- Environment activation instructions

### Installing the MIBI Analysis Package

To use the helper functions provided in this repository:

```bash
# Clone the repository
git clone https://github.com/HartmannLab/MIBI-Analysis_Hamburg_CRC_TMA_2024.git
cd MIBI-Analysis_Hamburg_CRC_TMA_2024

# Install the package in development mode
pip install -e .

# Or install dependencies manually for the scverse environment
conda env create -f environments/scverse.yml
conda activate liana
pip install -e .
```

## Data Requirements

This analysis requires:
- Raw MIBI imaging data
- Segmentation masks (cellpose-generated)
- Clinical metadata
- Transformed cell tables

*Note: Due to patient data privacy and size constraints, raw data files are not included in this repository. They will be uploaded on data sharing platforms as soon as possible.*

## Getting Started

1. **Set up environments**: Follow instructions in `environments/` to install required computational environments
2. **Data preparation**: Place your data files in the appropriate `data/` subdirectories
3. **Preprocessing**: Start with notebooks in `notebooks/preprocessing/` using the Docker environment
4. **Cell annotation**: Proceed to `notebooks/lineage/` for cell type identification
5. **Multicellular analysis**: Use notebooks in `notebooks/multicellular/` for higher-order analysis
6. **Visualization**: Generate figures using notebooks in `notebooks/viz/` and `notebooks/metadata/`

## Key Analysis Steps

Metadata preparation (`notebooks/metadata`)  

Preprocessing (`notebooks/preprocessing`):

* Image extraction and channel compensation with Toffy
* Gradient correction
* CLAHE correction of segmentationg channels, image segmentation with CellPose and feature extraction
* Cell and mask filtering
* Data transformation

Quality control (`notebooks/qc`)  

Visualization (`notebooks/viz`)  

Cell type annotation (`notebooks/lineage`)  

* Major lineage marker annotation (`Annotation_gating_Hamburg`)
* Curation and agreement with PyFlowSOM (`Pyflowsome_clustering_Hamburg`)
* Scyan clustering (`annotation_scyan`)
* Consensus lineage annotations (`annotation_consensus`)
* Analysis and visualization of cell type abundance (`CellTypeViz`)

## Reproducibility

To ensure reproducibility:
- Use the specified computational environments
- Follow the notebook execution order outlined above
- Check environment specifications in `environments/` directory
- Verify data preprocessing steps before downstream analysis

## Contact

For questions about the analysis pipeline or code, please open an issue in this repository.

## Citation

If you use this code or analysis pipeline, please cite our manuscript: [Citation to be added upon publication].


