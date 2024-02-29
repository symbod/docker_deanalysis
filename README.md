
# README for Differential Expression (DE) Analysis Pipeline

## Overview
This repository contains a Nextflow pipeline for performing differential expression (DE) analysis. It's designed to compare various conditions and timepoints in proteomics data, identifying differentially expressed elements between these groups.

## Container
The corresponding docker container, which is called in the Nextflow file, can be found here: https://hub.docker.com/repository/docker/kadam0/deanalysis/general

## Usage
To run the pipeline, you need to specify three flags:

1. `--meta_file`: This flag should be followed by the path to your metadata file. This file contains essential information about your samples, such as experimental conditions and timepoints.

2. `--count_file`: This flag requires the path to your count file, which contains the raw counts data for your analysis.

3. `--output`: Use this flag to specify the output directory where the results will be stored.

Example command:
```
./nextflow main.nf --meta_file path/to/meta_file.csv --count_file path/to/count_file.csv --output output_directory
```

## Output Description
The pipeline generates several output files and directories:

- **`de_data` Directory**: Contains the raw results of the DE analysis for each pairwise comparison. Files are named following the pattern `DEdata.TX_CX-TY_CY.txt`, where `TX` and `TY` represent specific timepoints, and `CX` and `CY` represent conditions. The underscore `_` connects attributes of the same type (e.g., timepoint or condition), while the hyphen `-` separates the two types being compared.

- **Metadata Encoding Files**: 
  - `Metadata_encoding.txt.json`: Describes the encoding used in the analysis, explaining what codes like `T1`, `C1`, etc., represent.
  - `Metadata_reverse_encoding.txt.json`: Provides reverse mapping for the encoding, aiding in understanding the raw data and results.

- **Overview Files**:
  - `DE_overview_XXX_sum.csv`: A matrix summarizing the total number of differentially expressed elements for each pairwise comparison.
  - `DE_overview_XXX_updown.csv`: In this matrix, the upper right triangle shows the number of upregulated elements, and the lower left triangle shows the number of downregulated elements for each comparison.

## Getting Started
To get started, clone this repository and ensure you have Nextflow installed. Prepare your metadata and count files according to the required format and run the command with the appropriate paths to your files. Check the output directory for results and detailed analysis.

---

**Note**: This README provides a general guide. Users should adjust paths and file names according to their specific project structure and requirements.
