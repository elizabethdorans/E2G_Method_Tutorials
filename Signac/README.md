parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object with RNA, ATAC, and peak assays [.rds]")
parser$add_argument("--signac_output_dir", default = ".",
                    help = "Path to directory for output files")

# Signac

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

`bash Signac_postprocessing.sh <input_file> <macs2_peaks_bedfile>`

<input_file>: A .tsv file containing eak-peak co-accessibility predictions.\
<macs2_peaks_bedfile>: A bedfile with ATAC peaks (output of `../seurat_object_preprocessing.R`).\
<promoter_bedfile>: [DEFAULT ../TSS_plusminus500bp.bed] A bedfile specifying gene promoter regions (columns chr, start, end, gene). The default file contains the regions +/- 500 bp from the gene TSS specified in ../CollapsedGeneBounds.hg38.bed.
                    
Outputs: 

1) Peak-gene link predictions at cicero_peak_gene_links.tsv (in the same folder as <input_file>).

## Step 3: Postprocessing for IGVF portal (see main E2G_Method_Tutorials folder)
