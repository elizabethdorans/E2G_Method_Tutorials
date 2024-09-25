### E2G Method Tutorials
This repo contains tutorials for linking enhancers (ATAC peaks) to genes using single-cell RNA + ATAC multiome data. Data preprocessing and peak-gene linking scripts are included for the following methods: Signac (Stuart 2021 Nature Methods), Cicero (Pliner 2018 Mol Cell), ArchR (Granja 2021 Nature Genetics). A tutorial for SCENT (Sakaue 2024 Nature Genetics) is forthcoming.

Clone the repository using the following command: 

`git clone https://github.com/elizabethdorans/E2G_Method_Tutorials.git`

See below (and method-specific subfolders) for steps to run single-cell peak-gene linking methods on single-cell multiome data:

## Step 0: Pre-process single-cell multiome data.

The script `seurat_object_preprocessing.R` takes as input single-cell multiome data files (RNA count matrix, ATAC fragment files, etc.), assembles a Seurat object, and calls ATAC peaks.

Example command: 

`Rscript seurat_object_preprocessing.R --rna_matrix <rna_matrix.mtx> --rna_matrix_barcodes <rna_matrix_barcodes> --rna_matrix_genes <rna_matrix_genes> --atac_fragments <atac_fragments.tsv.gz> --filtered_barcodes <filtered_barcodes> --output_dir <output_dir>`

<rna_matrix>: A .mtx file containing a cell x gene count matrix.
<rna_matrix_barcodes>: A .txt file containing cell barcodes (one per line) corresponding to <rna_matrix>.
<rna_matrix_genes>: A .txt file containing gene names (one per line) corresponding to <rna_matrix>.
<atac_fragments>: ATAC fragment file (.tbi file with same file basename must be in the same folder).
<filtered_barcodes>: [OPTIONAL] File with subset of cell barcodes (in column labeled 'barcode') and, optionally, additional metadata columns. The Seurat object will only include these cells.
<output_dir>: Path to folder where outputs will be saved.
                    
Outputs: 

1) Seurat object at <output_dir>/seurat_object.rds
2) ATAC peak bedfile at <output_dir>/macs2_peaks.bed

## Running E2G methods

See README.md files in method-specific folders for further steps in running each method!

## Post-processing for IGVF portal:

The script `postprocessing_for_IGVF_portal.R` takes as input peak-gene link predictions, restricts to a given gene universe, and produces a file with format appropriate for the IGVF portal.

Example command: 

`Rscript postprocessing_for_IGVF_portal.R --input_file <input_file> --output_file <output_file>  --genes_file CollapsedGeneBounds.hg38.bed --cell_type <cell_type> --method <method> --version <version>`

<input_file>: A .tsv file containing E2G method predictions (includes at minimum columns 'peak', 'gene', and 'Score').
<output_file>: A .tsv file containing containing E2G method predictions restricted to gene universe and reformatted for IGVF portal.
<genes_file>: [DEFAULT CollapsedGeneBounds.hg38.bed] Gene universe file (inclues gene symbol in 'name' column and Ensembl ID in 'Ensembl_ID' column). The default gene universe file CollapsedGeneBounds.hg38.bed was obtained from https://github.com/EngreitzLab/CRISPR_comparison/blob/main/resources/genome_annotations/CollapsedGeneBounds.hg38.bed.
<cell_type>: Cell type used to generate predictions (for 'CellType' column).
<method>: E2G method used to generate predictions (for header).
<version>: Version of E2G method used to generate predictions (for header).
                    
Outputs: 

1) Reformatted predictions at <output_file>