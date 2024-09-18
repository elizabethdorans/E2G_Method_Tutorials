### Inputs
# Path to Seurat object with RNA + ATAC peak counts

## Outputs
# Peak x cell matrix for running Cicero

suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="path to Seurat object with RNA, ATAC, and peak assays")

args <- parser$parse_args()

seurat_object = args$seurat_object
out_dir = dirname(seurat_object)

# Read in Seurat object
data = readRDS(seurat_object)

# Cicero input: save peak cell matrix
cicero_in_dir = sprintf("%s/cicero/input", out_dir)
dir.create(cicero_in_dir, recursive = TRUE)

accessibility.data <- GetAssayData(object = data, assay = "peaks", slot = "data")
accessibility.data <- 1 * (accessibility.data > 0)

writeMM(accessibility.data, file = sprintf("%s/matrix.mtx", cicero_in_dir))
write(rownames(accessibility.data), file = sprintf("%s/peaks.tsv", cicero_in_dir), ncolumns = 1)
write(colnames(accessibility.data), file = sprintf("%s/barcodes.tsv", cicero_in_dir), ncolumns = 1)