### Inputs
# Path to directory with Cicero inputs (see cicero_preprocessing.R)

## Outputs
# .tsv file with Cicero connections (peak-peak partial correlations)

suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(cicero))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

set.seed(1234)

genome_size_file <- "/n/groups/price/elizabeth/data/hg38.chrom.sizes"

parser <- ArgumentParser()

parser$add_argument("--input_dir",
    help="folder with cicero input subfolder, cell type tsv, etc.")
parser$add_argument("--cell_type_taskfile", 
    help="[OPTIONAL] focal cell types (will only compute links within cells of these cell types)")
parser$add_argument("--cell_types_by_barcode",
    help="file with cell barcodes and types")

args <- parser$parse_args()

input_dir = args$input_dir
cell_type_taskfile = args$cell_type_taskfile
cell_types_by_barcode = args$cell_types_by_barcode

# Define cell types and file label
if (!is.null(cell_type_taskfile)) {
    cell_types = paste(readLines(cell_type_taskfile))
    cell_types = gsub(" ", "-", cell_types)
    cell_types_label = strsplit(tail(strsplit(cell_type_taskfile, "/")[[1]], 1), "[.]")[[1]][1]

    print(cell_types_label)
    print(cell_types)
}

# create a dated subdir in results dir for output
mtx_folder = sprintf("%s/cicero/input", input_dir)
out_dir = sprintf("%s/cicero/pgl", input_dir)
dir.create(out_dir, recursive = TRUE)

outfile = sprintf("%s/cicero_connections.tsv", out_dir)
print(outfile)

# Read in matrix data using the Matrix package
inmtx <- Matrix::readMM(paste0(mtx_folder, "/matrix.mtx"))
indata <- 1*inmtx
# Binarize the matrix
indata@x[indata@x > 0] <- 1

# Format cell info
cellinfo <- read.table(paste0(mtx_folder, "/barcodes.tsv"))
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# Format peak info
peakinfo <- read.table(paste0(mtx_folder, "/peaks.tsv"))
row.names(peakinfo) <- peakinfo$V1
names(peakinfo) <- "peaks"

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))

# Subset by cell type (if applicable)
if (!is.null(cell_type_taskfile)) {
    cell_idents = read.table(cell_types_by_barcode, sep = "\t")
    colnames(cell_idents) = c("cell_type", "cell_barcode")
    focal_barcodes = subset(cell_idents, cell_type %in% cell_types)$cell_barcode
    input_cds <- input_cds[, focal_barcodes]
    print(sprintf("subsetted to %s cells", length(colnames(input_cds))))
}

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
input_cds <- monocle3::detect_genes(input_cds)

input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")

umap_coords <- reducedDims(input_cds)$UMAP

# Create Cicero CDS
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Read chromosome length information
chromosome_length <- read.table(genome_size_file)

conns <- run_cicero(cicero_cds, chromosome_length)
conns = conns[,c("Peak1", "Peak2", "coaccess")]

write.table(conns, outfile, sep = "\t", quote = F, row.names = F)