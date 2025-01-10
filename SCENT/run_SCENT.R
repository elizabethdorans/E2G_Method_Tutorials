suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object with RNA, ATAC, and peak assays")
parser$add_argument("--candidate_link_file",
    help = "Path to file candidate peak-gene links")
parser$add_argument("--scent_output_dir",  default = ".",
    help = "Path to directory for output files")
parser$add_argument("--scent_output_file",
    help = "Path to file to write SCENT output")

args <- parser$parse_args()

seurat_object = args$seurat_object
candidate_link_file = args$candidate_link_file
scent_output_dir = args$scent_output_dir
scent_output_file = args$scent_output_file

# Create output directory if needed
print(sprintf("Output directory: %s", scent_output_dir))
if (!dir.exists(scent_output_dir)) {
    dir.create(scent_output_dir, recursive = TRUE)
}

# Create outfile name
if (scent_output_file != NULL) {
  split = strsplit(basename(candidate_link_file), "[.]")[[1]]
  scent_output_file <- sprintf(sprintf("%s/%s.tsv", scent_output_dir, paste(split[1:length(split)-1], collapse = ".")))
}

# Read in Seurat object
print("Reading in Seurat object!")
data = readRDS(seurat_object)

