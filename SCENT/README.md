# SCENT Tutorial

This folder contains code for peak-gene linking using SCENT (Sakaue 2024 Nature Methods). SCENT code and tutorial can be found at [https://github.com/immunogenomics/SCENT].

\***If you wish to run SCENT without generating bootstrap p-values (for example, if you are generating features for [pgBoost](https://github.com/elizabethdorans/pgBoost/tree/main)), then you must:
- Download the development version of SCENT using:
```
devtools::install_github("immunogenomics/SCENT", ref="dev")
```
- Follow the additional instructions marked \*** in this README file.

## Step 1: Creating candidate peak-gene link files for SCENT

SCENT requires as input a text file specifying candidate peak-gene pairs to test for association. Because the SCENT algorithm is computationally intensive, it is often wise to split candidate peak-gene pairs into multiple such files on which SCENT can be run in parallel. 

The script `create_candidate_link_files.R` takes as input a Seurat object and generates candidate peak-gene link files.

Example command: [~3 minutes]

`Rscript create_candidate_link_files.R --seurat_object <seurat_object> --number_candidate_links_per_file <number_candidate_links_per_file> --candidate_link_file_output_dir <candidate_link_file_output_dir> --tss_coordinates_file <tss_coordinates_file>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of ../seurat_object_preprocessing.R).\
<number_candidate_links_per_file>: The maximum number of candidate peak-gene links (one linke per row) in each candidate peak-gene link file (per run of SCENT). Time to test each candidate peak-gene link can vary widely. It is reasonable to specify ~1000 candidate links per file (for ~1 day per SCENT run, see below) and adjust as needed. 
- \***If you wish to run SCENT without generating bootstrap p-values, it is reasonable to specify ~70,000 links per file (for ~1 day per SCENT run, see below) and adjust as needed.

<candidate_link_file_output_dir>: Path to folder where candidate peak-gene link files will be saved.\
<tss_coordinates_file>: [OPTIONAL] Path to file with TSS coordinates. Default is "../TSS.txt" which is included in the E2G_Method_Tutorials repo. NOTE: if you are running this script from another directory, you must specify the path, e.g. --tss_coordinates_file E2G_Method_Tutorials/TSS.txt.
                    
Outputs: 

1) Candidate peak-gene link files at <candidate_link_file_output_dir>/chunk_\*.txt (* will be replaced with a number).

## Step 2: Running SCENT

<inv>2A: Running the SCENT algorithm on on one candidate peak-gene link file<\inv>

**NOTE: Step 2A is optional, but recommended to test the SCENT code on a smaller scale before proceeding to step 2B.**

The script `run_SCENT.R` takes as input a Seurat object and candidate peak-gene link file and generates SCENT predictions.

Example command: [~1 day, ~10G; varies depending on number of cells and number of candidate peak-gene links]

`Rscript run_SCENT.R --seurat_object <seurat_object> --candidate_link_file <candidate_link_file> --scent_output_file <scent_output_file>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of `../seurat_object_preprocessing.R`).\
<candidate_link_file>: A tab-delimited text file with columns "peak" and "gene" (same as <candidate_link_file_output_dir> argument supplied to `create_candidate_link_files_from_Seurat_object.R`).\
<scent_output_file>: Path to file where SCENT output file will be saved.\

\***If you wish to run SCENT without generating bootstrap p-values, provide the `--skip_bootstrap` flag. Example command:

`Rscript run_SCENT.R --seurat_object <seurat_object> --candidate_link_file <candidate_link_file> --scent_output_file <scent_output_file> --skip_bootstrap`
                    
Outputs: 

1) SCENT peak-gene link predictions (scent_output_file).

<inv>2B: Running the SCENT algorithm on multiple candidate peak-gene link files in parallel<\inv>

The script `run_SCENT_multiple_chunks.sh` takes as input a Seurat object and *a folder of* candidate peak-gene link files and generates SCENT predictions.
**ATTENTION: Line 21 contains an 'sbatch' commmand to submit batch jobs to Slurm on a remote cluster with 10G of memory and 1 day of runtime each. Edit this line as appropriate.**

Example command: 

`bash run_SCENT_multiple_chunks.sh <seurat_object> <candidate_link_folder> <scent_output_dir> <script>`

<seurat_object>: A Seurat object containing ATAC, RNA, and peak data (output of `../seurat_object_preprocessing.R`).\
<candidate_link_folder>: A folder of tab-delimited text files with columns "gene" and "peak" (output of `create_candidate_link_files.R`).\
<scent_output_dir>: Path to file where SCENT output files will be saved.
\<script>: [OPTIONAL] Path to run_SCENT.R script (only need to supply if running run_SCENT_multiple_chunks.sh from outside the E2G_Method_Tutorials/SCENT directory).

\***If you wish to run SCENT without generating bootstrap p-values, run the the `run_SCENT_multiple_chunks_skip_bootstrap.sh` script instead. Example command:

`bash run_SCENT_multiple_chunks_skip_bootstrap.sh <seurat_object> <candidate_link_folder> <scent_output_dir> <script>`

Outputs: 

1) SCENT peak-gene link predictions in <scent_output_dir>/\*.tsv (* will be replaced with the basenames of the candidate peak-gene link files).
    
## Step 3: Concatenating SCENT predictions
   
The script `concatenate_scent_links.py` takes as input a folder of SCENT predictions (with a \*.tsv suffix) and generates a single file of peak-gene link predictions.

Example command:

`python concatenate_scent_links.py --scent_predictions_dir  <scent_output_dir> --output_file <output_file>`

<scent_predictions_dir>: A folder of tab-delimited files (\*.tsv) files containing output of the SCENT algorithm (same as <scent_output_dir> argument supplied to `run_SCENT_multiple_chunks.sh`.\
<output_file>: Path where tab-delimited file (recommend \*.tsv as the specified suffix) with concatenated SCENT predictions will be saved.
                    
Outputs: 

1) Single file of SCENT peak-gene link predictions (output_file).

## Step 4: Postprocessing for IGVF portal (see main E2G_Method_Tutorials folder)
