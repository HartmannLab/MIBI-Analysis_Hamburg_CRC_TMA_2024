# MIBI Analysis – Hamburg CRC TMA 2024

This repository compiles the scripts used to process and analyze the data for the spatial metabolic protein atlas of colorectal cancers.

## Study authors

Teresa Glauner, Loan Vulliard, Felix Hartmann and friends.

## Data availability

Data will be made available at a later point.

## Analysis workflow

Metadata preparation (`notebooks/metadata`)  

Preprocessing (`notebooks/preprocessing`):
* Image extraction and channel compensation with Toffy
* Gradient correction
* CLAHE correction of segmentationg channels, image segmentation with CellPose and feature extraction
* Cell and mask filtering
* Data transformation

Quality control (`notebooks/qc`)  

Visualization (`notebooks/viz`)  

Cell type annotation (ongoing)  


## Code overview

Jupyter notebooks can be found in the `notebooks` folder.

## Environment

Most scripts can be run using the [koalive/mibiprofiling](https://hub.docker.com/r/koalive/mibiprofiling/tags) Docker image.



