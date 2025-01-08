# SCENT Tutorial

This folder contains code for peak-gene linking using SCENT (Sakaue 2024 Nature Methods). SCENT code and tutorial can be found at [https://github.com/immunogenomics/SCENT].

## Step 1: Creating candidate peak-gene link files for SCENT

SCENT requires as input a text file specifying candidate peak-gene pairs to test for association. Because the SCENT algorithm is computationally intensive, it is often wise to split candidate peak-gene pairs into multiple such files on which SCENT can be run in parallel. 

The script `create_candidate_link_files_from_Seurat_object.R` takes as input a Seurat object and generates candidate peak-gene link files.

Example command:

`Rscript create_candidate_link_files_from_Seurat_object.R --seurat_object <seurat_object> --number_candidate_links_per_file <number_candidate_links_per_file> --candidate_link_file_output_dir <candidate_link_file_output_dir>`

<seurat_object>: A preprocessed Seurat object (output of `Signac_preprocessing.R`).\
<signac_output_dir>: Path to folder where outputs will be saved (one file per chromosome).


Example command: [~1 hour, ~30G]

`Rscript Signac_preprocessing.R --seurat_object <seurat_object> --signac_output_dir <signac_output_dir>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of `../seurat_object_preprocessing.R`).\
<signac_output_dir>: Path to folder where outputs will be saved.
                    
Outputs: 

1) Processed Seurat object at <signac_output_dir>/seurat_object_signac_preprocessed.rds

## Step 2: Running Signac

The script `run_Signac.sh` takes as input a pre-processed Seurat object and runs Signac peak-gene linking.
**ATTENTION: Lines 12 and 22 contain 'sbatch' commmands to submit batch jobs to Slurm on a remote cluster. Edit these lines as appropriate.**

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
