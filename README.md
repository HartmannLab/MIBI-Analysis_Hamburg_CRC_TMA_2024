# MIBI Analysis – Hamburg CRC TMA 2024

This repository compiles the scripts used to process and analyze the data from the *insert a fancy project name here*.

## Study authors

Teresa Glauner, Loan Vulliard, Felix Hartmann and friends.

## Data availability

Data will be made available at a later point

## Analysis workflow

* Image extraction and channel compensation with Toffy
* Gradient correction
* CLAHE correction of segmentationg channels, image segmentation with CellPose and feature extraction
* Cell and mask filtering
* Cell type annotation (ongoing)

## Code overview

Jupyter notebooks can be found in the `notebooks` folder.

## Environment

Several scripts can be run using the [koalive/mibiprofiling](https://hub.docker.com/r/koalive/mibiprofiling/tags) Docker image.



