# Cicero Tutorial

This repo contains code for peak-gene linking using Cicero (Pliner 2018 Mol Cell).

## Step 1: Run Cicero co-accessibility

The script `run_Cicero.R` takes as input an assembled Seurat object and computes Cicero peak-peak co-accessibilities.

**NOTE: In order to use the included `hg38.chrom.sizes` file, run `run_Cicero.R` from this folder (`Cicero`) or specify the path to this file from your location as an argument to `--genome_size_file`.**

Example command: [~8 hours, ~30G]

`Rscript run_Cicero.R --seurat_object <seurat_object> --cicero_output_dir <cicero_output_dir>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of `../seurat_object_preprocessing.R`).\
<cicero_output_dir>: Path to folder where outputs will be saved.
                    
Outputs: 

1) Peak-peak co-accessibility predictions at <cicero_output_dir>/cicero_connections.tsv

## Step 2: Postprocessing into peak-gene links

The script `Cicero_postprocessing.sh` takes as input Cicero peak-peak co-accessibilities, ATAC peaks, and promoters overlapping ATAC peaks, and generates peak-gene link predictions.

Example command: 

`bash Cicero_postprocessing.sh <input_file> <macs2_peaks_bedfile>`

<input_file>: A .tsv file containing eak-peak co-accessibility predictions.\
<macs2_peaks_bedfile>: A bedfile with ATAC peaks (output of `../seurat_object_preprocessing.R`).\
<promoter_bedfile>: [DEFAULT ../TSS_plusminus500bp.bed] A bedfile specifying gene promoter regions (columns chr, start, end, gene). The default file contains the regions +/- 500 bp from the gene TSS specified in ../CollapsedGeneBounds.hg38.bed.
                    
Outputs: 

1) Peak-gene link predictions at cicero_peak_gene_links.tsv (in the same folder as <input_file>).

## Step 3: Postprocessing for IGVF portal (see main E2G_Method_Tutorials folder)
