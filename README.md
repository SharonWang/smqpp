# BGlab Smart-Seq2 preprocessing toolkit

## Introduction

smqpp is the preprocessing pipeline for Smart-Seq2 data, specifically for datasets generated from Gottgens lab. The code was adpoted from [bglab](https://github.com/wjawaid/bglab) package developed by Wajid Jawaid.

The package contains the following steps:
- read_in_files: Read in count and QC inputs and format them into anndata object
- reformat_meta: Reformat metatable to keep all versions consistent (Due to different versions of metadata spread sheet from google drive)
- smartseq_qc: bglab equivalent quality control
- normalise_data: Data normalisation using DESeq2 method
- tech_var: Highly variable gene (HVG) calculation using [Brennecke et. al
](https://www.nature.com/articles/nmeth.2645?proof=trueInJun) method
- plot_tech_var: Plot the HVG prediction
- plot_ma: MAplot for rank_genes_group from [Scanpy](https://github.com/theislab/scanpy) and select significant genes with high confidence

## Installation

smqpp depends on numpy, matplotlib, pandas, anndata, scipy and statsmodels. The package is available on pip and can be easily installed as follows:

	pip install smqpp
or
	download the file from github 
	tar zxvf smqpp
	cd smqpp
	pip install .

## Usage and Documentation

The smqpp should be fairly simple to use and it is based on Scanpy's AnnData object:

	import smqpp
	smqpp.read_in_files(...)

## Example Notebooks

To be provided

## Contact

If there are any issues, please contact xw251@cam.ac.uk.



