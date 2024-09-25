# ArchR Tutorial

This repo contains code for peak-gene linking using ArchR (Granja 2021 Nature Genetics).

## Step 1: Run ArchR peak-gene linking.

The script `run_ArchR.R` takes as input ATAC fragment files, an assembled Seurat object, and runs ArchR peak2gene.

Example command: 

`Rscript run_ArchR.R --seurat_object <> --atac_fragments <> --archr_output_dir <>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of `../seurat_object_preprocessing.R`). Must also contain a 'barcode' column in metadata.
<atac_fragments>: ATAC fragment file (.tbi file with same file basename must be in the same folder).
<archr_output_dir>: Path to folder where outputs will be saved.
                    
Outputs: 

1) Peak-gene link predictions at <archr_output_dir>/archr_peak_gene_links.tsv
2) ArchR project at <output_dir>/project

## Step 2: Postprocessing for IGVF portal (see main E2G_Method_Tutorials folder)
