# BGlab Smart-Seq2 preprocessing toolkit

## Introduction

smqpp is the preprocessing pipeline for Smart-Seq2 data, specifically for datasets generated from Gottgens lab. The QC part of code was adpoted from [bglab](https://github.com/wjawaid/bglab) package developed by Wajid Jawaid.

The package contains the following steps:
1. **<ins>Preanalysis</ins>**
	- **generate_feature_table**: If gene feature table not available then this can be generate using this function
	- **read_in_files**: Read in count and QC inputs and format them into anndata object
	- **reformat_meta**: Reformat metatable to keep all versions consistent (Due to different versions of metadata spread sheet from google drive)
	- **smartseq_qc**: bglab equivalent quality control
	- **normalise_data**: Data normalisation using DESeq2 method
	- **quantile_norm**: Quantile normalisation
	- **quantile_norm_log**: Log quantile normalisation
	- **downsampling_norm**: Downsampling normalisation (Not recommanded for TenX as it will shrink more counts to 0)
	- **tech_var**: Highly variable gene (HVG) calculation using [Brennecke et. al
](https://www.nature.com/articles/nmeth.2645?proof=trueInJun) method
	- **plot_tech_var**: Plot the HVG prediction
	- **detect_outlier_cells**: filter out cells that effect the selection of HVGs

2. **<ins>Differential expression analysis</ins>**
	- **plot_ma**: MAplot for rank_genes_group from [Scanpy](https://github.com/theislab/scanpy) and select significant genes with high confidence

3. **<ins>Pseudotime time analysis</ins>**
	- **GeneExp_LLR_test**: Likelihood ratio test to select genes that differentially expressed along pseudotime. Linear models were fitted between log norm exp and smoothed PT by applying natural spline.
	- **plot_genes_along_pt**: Plotting out gene expression pattern along PT. Gene exp was smoothed using Guassian filter.

4. **<ins>Projection</ins>**
	- **quick_neighbors**: Neighbors calculation adpoted from [scanpy](https://github.com/theislab/scanpy). Two constraints applied: 1) reference cells only allow neighbors between themselves; 2) new cells only allow neighbors with reference cells
	- **quick_umap**: Similar to the ingest function in [scanpy](https://github.com/theislab/scanpy). Umap was calculated using [umap](https://github.com/lmcinnes/umap) python package. Parameters used as [scanpy](https://github.com/theislab/scanpy) defaults.
	- **quick_umap_proj**: Projection of new data onto reference data

5. **<ins>3d plots</ins>**
	- **plot_3d**: Generate 3d plots from anndata object as the projection='3d' function does not work properly in the latest scanpy due to matplotlib issues
	
6. **<ins>Pathway analysis</ins>**:
	- **pathway_score_cal**: Calculate geometric mean for each terms in the databse for each cell, which can be used to color the defined layout
	- **pathway_analysis**: Calculate the enriched database terms for a given gene set using hypergeometric test.

## Installation

smqpp depends on numpy, matplotlib, pandas, anndata, scipy and statsmodels. The package is available on pip and can be easily installed as follows:

	pip install smqpp
or
	
	download the file from github using git clone
	tar zxvf smqpp
	cd smqpp
	pip install .

## Usage and Documentation

The smqpp should be fairly simple to use and it is based on Scanpy's AnnData object:

	import smqpp
	smqpp.read_in_files(...)

## Example Notebooks

Examples can be found in the following folders: 
1. **<ins>[GSK analysis](https://github.com/SharonWang/GSK_analysis)</ins>**
2. **<ins>[Patel study](https://github.com/SharonWang/Patel_Study)</ins>**

## Contact

If there are any issues, please contact xw251@cam.ac.uk.



