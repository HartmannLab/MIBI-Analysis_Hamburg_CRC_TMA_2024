# Data

The cell-level data of the study are deposited on Zenodo: [10.5281/zenodo.18504904](https://doi.org/10.5281/zenodo.18504904). The supplementary tables are deposited separately, [10.5281/zenodo.17280877](https://doi.org/10.5281/zenodo.17280877). Supplementary Table 3 lists the antibody panel, with the target, clone, metal and mass of each channel. This page describes the columns of the three deposited files, how they connect, and how to load them.

## How the files connect

| File | One row per | Key | Made by |
|---|---|---|---|
| `cell_table.csv` | cell, after artefact filtering | `fov`, `label` | [`preprocessing/FilterAndViz_Hamburg.ipynb`](notebooks/preprocessing/FilterAndViz_Hamburg.ipynb) |
| `adata_consensus_cell_types.h5ad` | cell, with cell type annotations | `obs["fov"]`, `obs["label"]` | [`lineage/annotation_consensus.ipynb`](notebooks/lineage/annotation_consensus.ipynb) |
| `summary_clinical_data_modified.csv` | TMA core, which is one patient | `TMA Localisation` | [`metadata/Making a comprehensive list.ipynb`](notebooks/metadata/Making%20a%20comprehensive%20list.ipynb) |

Each core of the tissue microarray was imaged once, as one field of view of 400 × 400 µm at 1024 × 1024 pixels, so 390 nm per pixel. The field of view is named after the position of its core on the array, such as `A1a`.

- To link cells to patients, match `fov` in the cell tables to `TMA Localisation` in the clinical table.
- To link a cell to the segmentation masks, use `fov` and `label`. Together they identify a cell uniquely, and `label` is the value of the cell in the mask of its core.
- The two cell tables hold the same cells in the same order, except that the AnnData file leaves out core `C6h`, where most cells are badly segmented. The row names of the AnnData file are the row numbers of `cell_table.csv`.

`cell_table.csv` holds 579,947 cells from 470 cores, and the AnnData file 579,011 cells from 469 cores. The clinical table lists all 522 cores of the array, including the cores that were damaged or not imaged. The 11 control cores, whose names start with `E1` to `E4`, have no clinical row. The analyses use `E3` and `E4` as normal colon, labelled `Colon-no.`, and no comparison includes `E1` and `E2`.

## `cell_table.csv`

One row per segmented cell that passed the artefact filters: an area between 120 and 1500 pixels, an integrated nuclear signal above 2000, and a mean gold signal below 1000. Cores with fewer than 300 such cells are left out. The first, unnamed column is the row number.

The channel columns are named after the image files, with the suffix `.tiff`. Each holds the mean pixel value over the mask of the cell, with no further transformation. For the marker channels, the images are gradient-corrected, as the preprint describes: [`preprocessing/gradient_correction.ipynb`](notebooks/preprocessing/gradient_correction.ipynb) divides the rescaled images that toffy writes after Rosetta compensation by a gradient profile estimated for each channel. The values are therefore neither ion counts nor the toffy output, but on an arbitrary scale, and the technical channels below sit on a different scale.

| Columns | Content |
|---|---|
| `SMA.tiff`, `CD4.tiff`, `CD31.tiff`, `CD7.tiff`, `CD163.tiff`, `CD68.tiff`, `CD8.tiff`, `CD3e.tiff`, `MPO.tiff`, `CD20.tiff`, `HLADRa.tiff`, `CD14.tiff`, `CD45.tiff`, `PanCK.tiff`, `FoxP3.tiff` | Lineage markers, used for the cell type annotation |
| `CA9.tiff`, `CD98.tiff`, `CytC.tiff`, `MCT1.tiff`, `ASCT2.tiff`, `LDH.tiff`, `GS.tiff`, `GLS.tiff`, `ATP5A.tiff`, `CS.tiff`, `PKM2.tiff`, `GLUT1.tiff`, `ARG1.tiff`, `CPT1A.tiff`, `Ki67.tiff` | Metabolic markers and Ki67 |
| `MSH2.tiff`, `MSH6.tiff`, `PDL1.tiff`, `PD1.tiff`, `STING1.tiff`, `DCN.tiff`, `Vimentin.tiff` | Other functional and structural markers |
| `nuclei.tiff` | Nuclear channel, mass 115, which carries both the dsDNA and the histone H3 antibodies |
| `membrane.tiff` | Membrane channel, mass 176, which carries both the Na/K-ATPase and the HLA class I antibodies |
| `nuclei_clahe.tiff`, `membrane_clahe.tiff`, `GLUT1_clahe.tiff`, `ASCT2_clahe.tiff` | The same channels after CLAHE contrast equalisation, on a 0 to 255 scale. They served the segmentation only |
| `Au.tiff`, `Fe.tiff`, `Noodle.tiff`, `chan_39.tiff`, `chan_48.tiff`, `chan_69.tiff`, `chan_71.tiff` | Channels at masses that no antibody of the panel uses. The gold channel `Au` serves the artefact filter. The analyses use none of them |
| `centroid-0`, `centroid-1` | Position of the cell centre in pixels: `centroid-0` is the row, from the top of the image, and `centroid-1` the column, from the left |
| `area`, `convex_area` | Area of the cell mask and of its convex hull, in pixels. One pixel is 0.152 µm² |
| `perimeter`, `axis_major_length`, `axis_minor_length` | Perimeter of the mask and axis lengths of the ellipse with the same second moments, in pixels |
| `eccentricity` | Eccentricity of that ellipse, from 0 for a circle to 1 for a line |
| `label` | Value of the cell in the segmentation mask of its core |
| `fov` | Core, matching `TMA Localisation` in the clinical table |
| `AreaShape_Center_X`, `AreaShape_Center_X_1`, `AreaShape_Center_Y`, `AreaShape_Center_Y_1` | Copies of `centroid-1` and `centroid-0` in the format BioProfiling.jl expects |
| `NuclearIntensity` | `area` times `nuclei_clahe.tiff`, the integrated nuclear signal of the filter |
| `url` | Path of the nuclear image used to inspect the filters. It points to the machine that ran the filtering and has no use elsewhere |

The morphology comes from `skimage.measure.regionprops_table`, and the segmentation from Cellpose, as described in [`preprocessing/SegmentAndExtract.ipynb`](notebooks/preprocessing/SegmentAndExtract.ipynb).

## `adata_consensus_cell_types.h5ad`

An AnnData object with one observation per cell and the 15 lineage markers as variables. All marker values in it are transformed as in [`preprocessing/TransformMarkerData.ipynb`](notebooks/preprocessing/TransformMarkerData.ipynb): each channel is divided by its 99.9th percentile over all cells and multiplied by 10, then passed through arcsinh and divided by arcsinh(10). 99.9% of the values of a channel therefore lie between 0 and 1.

| Slot | Content |
|---|---|
| `X` | Transformed intensities of the 15 lineage markers, in the order of `var_names`: SMA, CD4, CD31, CD7, CD163, CD68, CD8, CD3e, MPO, CD20, HLADRa, CD14, CD45, PanCK, FoxP3 |
| `obs`, marker columns | Transformed intensities of the other 22 markers, under the names of `cell_table.csv` without `.tiff` |
| `obs["nuclei"]`, `obs["membrane"]` | The nuclear and membrane channels, not transformed |
| `obs`, morphology columns | `centroid-0`, `centroid-1`, `eccentricity`, `perimeter`, `convex_area`, `area`, `axis_major_length`, `axis_minor_length`, `label` and `fov`, as in `cell_table.csv` |
| `obs["annotation_t"]` | Cell type from manual gating, in [`lineage/Annotation_gating_Hamburg.ipynb`](notebooks/lineage/Annotation_gating_Hamburg.ipynb), 17 classes |
| `obs["annotation_scyan"]` | Cell type from Scyan, in [`lineage/annotation_scyan.ipynb`](notebooks/lineage/annotation_scyan.ipynb), 6 broad classes |
| `obs["annotation_consensus"]` | Final cell type, combining both annotations. The notebook assigns a cell type to each compatible pair of gating and Scyan labels, and every other pair becomes `Unclear`. This is the annotation of all analyses, with 15 classes: APC, B_cell, CAF, CD163_Macrophage, CD4_Tcell, CD68_Macrophage, CD8_Tcell, Cancer_cell, Endothelial_cell, Monocyte, NK_cell, Neutrophil, Other_immune_cell, T_reg_cell, Unclear |
| `obs["consensus"]` | An identical copy of `annotation_consensus` |
| `obsm["X_umap"]`, `obsp`, `uns` | UMAP embedding and neighbour graph computed by scanpy on `X` with 10 neighbours, and the colours of the annotation plots |

## `summary_clinical_data_modified.csv`

One row per core of the array, encoded in UTF-8. Column names and most values are in German, as recorded by the Institute of Pathology. The columns with `RED` in their name hold fewer, grouped categories than the raw field next to them.

The notebooks translate and group these values before using them, and they are the reference for how to use the clinical annotations:
- [`metadata/CohortPlotR.ipynb`](notebooks/metadata/CohortPlotR.ipynb) recodes the cohort description and writes the cohort summary, `clinical_summary_by_pT.csv`, deposited as `Supp_Table_2.csv`. It bins age and tumour diameter into thirds instead of using `Age Group` and `Tumor Size Group`.
- [`multicellular/ClinicalAssociations.ipynb`](notebooks/multicellular/ClinicalAssociations.ipynb) maps every annotation it tests to English categories in one place, the `annotations` dictionary in the cell that reads the clinical table.

| Column | Meaning |
|---|---|
| `TMA name`, `TMA sector` | Tissue microarray block and sector |
| `TMA Localisation` | Position of the core, the key to `fov` |
| `Tumorlokalisation COL10`, `Tumorlokalisation alle RED` | Tumour site, as recorded and in six groups |
| `Tumor Seite alle` | Side of the colon |
| `Tumortyp` | Tumour type |
| `pT raw`, `pN raw` | pT and pN categories, as reported |
| `pT group`, `pN group` | pT category without subdivisions, which is the stage of the analyses, and nodal status, `pN0` or `pN+` |
| `Grade raw`, `pM group` | Empty |
| `Tumordurchmesser (mm)` | Tumour diameter in mm |
| `Muzinöse Tumorkomponente` | Mucinous component of the tumour |
| `Infiltrationstyp` | Growth pattern at the invasion front |
| `Tumorbudding` | Tumour budding |
| `Peritumorale lymphozytäre Infiltration` | Peritumoral lymphocytic infiltration, free text |
| `Angioinvasion V RED`, `Lymphbahninvasion L` | Venous and lymphatic invasion |
| `MSI im Befund`, `MSI gesamt RED` | Microsatellite instability in the pathology report, and overall microsatellite status |
| `RAS Mutation RED`, `RAS Mutation RED 2` | Mutated RAS codon, and RAS status |
| `BRAF Mutation RED`, `HER2 Amplifikation` | BRAF status and HER2 amplification |
| `Age`, `Sex` | Age in years, and sex |
| `MIBI` | Whether the core was imaged, `acquired` or `no` |
| `Tumor Size Group`, `Age Group` | Tumour diameter in four bins of equal width, named by their lower edge in mm, and age below or from 70 |

## Loading the data

In Python, with `pandas` and `anndata`:

```python
import anndata as ad
import pandas as pd

adata = ad.read_h5ad("adata_consensus_cell_types.h5ad")
clinical = pd.read_csv("summary_clinical_data_modified.csv")

# Clinical annotation of every cell; the control cores have none
cells = adata.obs.merge(clinical, left_on="fov", right_on="TMA Localisation", how="left")

# Cell type composition of each core
composition = pd.crosstab(adata.obs["fov"], adata.obs["annotation_consensus"], normalize="index")

# Untransformed intensities and the technical columns, for the same cells
cell_table = pd.read_csv("cell_table.csv", index_col=0)
raw = cell_table.loc[adata.obs_names.astype(int)]
```

The notebooks expect the files in `data/` at the root of the repository. [`metabolism/MetabViz.ipynb`](notebooks/metabolism/MetabViz.ipynb) joins the AnnData file with the clinical table and writes `data/cell_table_with_types_stage.csv`, the cell table with stage labels that most later notebooks read.
