# Analyses environments

### 1. Docker Environment: `koalive/mibiprofiling`
**Used for:** Preprocessing steps and notebooks requiring Julia packages.  
Available from [DockerHub](https://hub.docker.com/r/koalive/mibiprofiling).
- `notebooks/preprocessing/FilterAndViz_Hamburg.ipynb`: data filtering with BioProfiling.jl
- `notebooks/preprocessing/FilterSegmentationMask.ipynb`: segmentation mask filtering
- `notebooks/viz/ExampleImages.ipynb`: pipeline visualization examples

### 2. Kasumi R Environment
**Used for:** R-based statistical analysis and visualization.  
Create conda environment from `kasumi.yml` and install R packages listed in `kasumi.csv`.
- `notebooks/metadata/CohortPlotR.ipynb`: clinical cohort visualization and statistics
- `notebooks/multicellular/BiplotFactors.ipynb`: factors associated with clinical variables
- `notebooks/spatial/KasumiColorectal.ipynb`: spatial metabolic analysis with Kasumi
- `notebooks/spatial/MistyLines.ipynb`: spatial lineage analysis with MISTy

### 3. scverse Environment
**Used for:** Factor analysis and enrichment analyses.  
Create conda environment from `scverse.yml`.
- `notebooks/multicellular/PrepareMuData.ipynb`: multi-modal data preparation
- `notebooks/multicellular/MOFACell.ipynb`: factor analysis
- `notebooks/multicellular/ClinicalAssociations.ipynb`: clinical annotations against pT and the factors

### 4. xgboost Environment
**Used for:** Machine learning predictions and PyTorch-based analyses.  
Create conda environment from `xgboost.yml`. Optimized for NVIDIA GPU and CUDA 12.3.  
- `notebooks/lineage/annotation_scyan.ipynb`: cell type annotation with Scyan
- `notebooks/lineage/annotation_consensus.ipynb`: consensus cell type annotation
- `notebooks/lineage/Annotation_gating_Hamburg.ipynb`: annotation by manual gating
- All notebooks with XGBoost predictions or deep learning models, such as `notebooks/spatial/StageXGB.ipynb` and `notebooks/spatial/EncodingClinicalInformation.ipynb`

### 5. scverse-squidpy Environment
**Used for:** Refitting the MuVI factors with other spatial encodings, including squidpy neighbourhood enrichment (`multicellular/SpatialFeatureContribution.ipynb`, `.utils/stage_fold_factors.py`).  
The scverse environment with squidpy added. `scverse_squidpy.yml` lists the exact versions, but conda cannot create it in one step: squidpy 1.8.2 declares zarr 3, which breaks anndata 0.11.2, so zarr stays below 3. Clone the scverse environment instead and add the two packages:
```bash
conda create --name liana-squidpy --clone liana
conda activate liana-squidpy
pip install squidpy==1.8.2
pip install "zarr<3"
```
pip then reports that squidpy, spatialdata and ome-zarr expect zarr 3. The squidpy functions the notebooks use run with zarr 2. Installing squidpy raises numpy from 1.26.4 to 2.0.2, and `SpatialFeatureContribution.ipynb` checks that the published MuVI factors still reproduce. The clone was made from the scverse environment after its torch was updated to 2.9.1 for CUDA 13, so `scverse_squidpy.yml` and `scverse.yml` differ in torch.

### 6. scimap_cozi Environment
**Used for:** COZI enrichment between cell types (`spatial/CellTypeCOZI.ipynb`).  
Create conda environment from `scimap_cozi.yml`. It installs the COZI fork of scimap, [SchapiroLabor/scimap_COZI](https://github.com/SchapiroLabor/scimap_COZI), at commit `191a1b1`, together with `natsort`, `mpmath` and `setuptools<81`, which the fork needs but does not declare.  
Packages installed in the user site directory, such as `~/.local/lib/python3.10`, would shadow those of the environment, so register its Jupyter kernel with user site packages disabled:
```bash
conda activate scimap_cozi
python -m ipykernel install --user --name scimap_cozi --display-name scimap_cozi --env PYTHONNOUSERSITE 1
```

*Note: to ensure accurate reproduction of the manuscript results, these are not minimal environments but the original ones used when designing and running the analyses.*