# SCENT Tutorial

This folder contains code for peak-gene linking using SCENT (Sakaue 2024 Nature Methods). SCENT code and tutorial can be found at [https://github.com/immunogenomics/SCENT].

## Step 1: Creating candidate peak-gene link files for SCENT

SCENT requires as input a text file specifying candidate peak-gene pairs to test for association. Because the SCENT algorithm is computationally intensive, it is often wise to split candidate peak-gene pairs into multiple such files on which SCENT can be run in parallel. 

The script `create_candidate_link_files_from_Seurat_object.R` takes as input a Seurat object and generates candidate peak-gene link files.

Example command:

`Rscript create_candidate_link_files_from_Seurat_object.R --seurat_object <seurat_object> --number_candidate_links_per_file <number_candidate_links_per_file> --candidate_link_file_output_dir <candidate_link_file_output_dir>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of ../seurat_object_preprocessing.R).\
<number_candidate_links_per_file>: The maximum number of candidate peak-gene links (one linke per row) in each candidate peak-gene link file (per run of SCENT). Time to test each candidate peak-gene link can vary widely. It may be wise to try running SCENT with 100-300 candidate links per file and later reduce if runtime is too long.\
<candidate_link_file_output_dir>: Path to folder where candidate peak-gene link files will be saved.
                    
Outputs: 

1) Candidate peak-gene link files at <candidate_link_file_output_dir>/chunk_\*.txt (* will be replaced with a number).

## Step 2: Running SCENT

* If you wish to submit one candidate peak-gene link file to the SCENT algorithm *

The script `run_SCENT.sh` takes as input a Seurat object and candidate peak-gene link file and generates SCENT predictions.
**ATTENTION: Lines XX and XX contain 'sbatch' commmands to submit batch jobs to Slurm on a remote cluster. Edit these lines as appropriate.**

Example command: 

`bash run_SCENT.sh <seurat_object> <candidate_link_file> <scent_output_dir>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of ../seurat_object_preprocessing.R).\
<candidate_link_file>: A tab-delimited text file with columns "peak" and "gene" (output of create_candidate_link_files_from_Seurat_object.R).\
<scent_output_dir>: Path to folder where output file will be saved.
                    
Outputs: 

1) SCENT peak-gene link predictions in <scent_output_dir>/\*.txt (* will be replaced with the basename of the candidate peak-gene link file).

* If you wish to submit a folder of candidate peak-gene link files to the SCENT algorithm in parallel *

The script `run_SCENT_multiple.sh` takes as input a Seurat object and *a folder of* candidate peak-gene link files and generates SCENT predictions.
**ATTENTION: Lines XX and XX contain 'sbatch' commmands to submit batch jobs to Slurm on a remote cluster. Edit these lines as appropriate.**

Example command: 

`bash run_SCENT_multiple.sh <seurat_object> <candidate_link_folder> <scent_output_dir>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of ../seurat_object_preprocessing.R).\
<<candidate_link_folder>>: A folder of tab-delimited text files (with suffix "*.txt") with columns "peak" and "gene" (output of create_candidate_link_files_from_Seurat_object.R).\
<scent_output_dir>: Path to folder where output files will be saved (one file per candidate peak-gene link file).

Outputs: 

1) SCENT peak-gene link predictions in <scent_output_dir>/\*.txt (* will be replaced with the basenames of the candidate peak-gene link files).

## Step 3: Postprocessing per-chromosome predictions

The script `SCENT_postprocessing.sh` takes as input a folder containing SCENT peak-gene link predictions and concatenates into a single file of peak-gene link predictions.

Example command: 

`Rscript SCENT_postprocessing.R --scent_output_dir <scent_output_dir>`

<scent_output_dir>: Path to a folder containing output files from SCENT algorithm.
                    
Outputs: 

1) Concatenated SCENT peak-gene link predictions at <scent_output_dir>/SCENT_peak_gene_links.tsv.

## Step 3: Postprocessing for IGVF portal (see main E2G_Method_Tutorials folder)
