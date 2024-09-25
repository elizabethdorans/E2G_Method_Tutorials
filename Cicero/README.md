# Cicero Tutorial

This repo contains code for peak-gene linking using Cicero (Pliner 2018 Mol Cell).

## Step 1: Run Cicero co-accessibility.

The script `run_Cicero.R` takes as input an assembled Seurat object and computes Cicero peak-peak co-accessibilities.

Example command: 

`Rscript run_Cicero.R --seurat_object <seurat_object> --cicero_output_dir <cicero_output_dir>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of `../seurat_object_preprocessing.R`).\
<cicero_output_dir>: Path to folder where outputs will be saved.
                    
Outputs: 

1) Peak-peak co-accessibility predictions at <cicero_output_dir>/cicero_connections.tsv

## Step 2: Postprocessing for IGVF portal (see main E2G_Method_Tutorials folder)
