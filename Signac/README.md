# Signac Tutorial

This repo contains code for peak-gene linking using Cicero (Pliner 2018 Mol Cell).

## Step 1: Preprocessing for Signac

The script `Signac_preprocessing.R` takes as input an assembled Seurat object and runs preprocessing steps needed for Signac peak-gene linking.

Example command: 

`Rscript Signac_preprocessing.R --seurat_object <seurat_object> --signac_output_dir <signac_output_dir>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of `../seurat_object_preprocessing.R`).\
<signac_output_dir>: Path to folder where outputs will be saved.
                    
Outputs: 

1) Processed Seurat object at <signac_output_dir>/seurat_object_signac_preprocessed.rds

## Step 2: Running Signac

The script `run_Signac.sh` takes as input a pre-processed Seurat object and runs Signac peak-gene linking.

Example command: 

`bash run_Signac.sh <seurat_object> <signac_output_dir>`

<seurat_object>: A preprocessed Seurat object (output of `Signac_preprocessing.R`).\
<signac_output_dir>: Path to folder where outputs will be saved (one file per chromosome).
                    
Outputs: 

1) Signac peak-gene link predictions (one file per chromosome) in <signac_output_dir>/chr*.tsv

## Step 3: Postprocessing per-chromosome predictions

The script `Signac_postprocessing.sh` takes as input a folder containing per-chromosome Signac peak-gene link predictions and generates a single file of peak-gene link predictions.

Example command: 

`Rscript Signac_postprocessing.R --input_folder <input_folder>`

<input_folder>: Path to a folder containing per-chromosome Signac peak-gene link predictions in the format chr*.tsv.
                    
Outputs: 

1) Concatenated Signac peak-gene link predictions at <input_folder>/signac_peak_gene_links.tsv.

## Step 3: Postprocessing for IGVF portal (see main E2G_Method_Tutorials folder)
