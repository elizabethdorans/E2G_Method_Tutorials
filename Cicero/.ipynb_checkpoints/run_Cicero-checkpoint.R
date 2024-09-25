suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(cicero))

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="Path to Seurat object with RNA, ATAC, and peak assays")
parser$add_argument("--cicero_output_dir",  default = ".",
    help = "Path to directory for output files")
parser$add_argument("--genome_size_file",  default = "hg38.chrom.sizes",
    help = "Path to file with chromosome sizes")

args <- parser$parse_args()

seurat_object = args$seurat_object
cicero_output_dir = args$cicero_output_dir
genome_size_file = args$genome_size_file

# Create output directory if needed
print(sprintf("Output directory: %s", cicero_output_dir))
if (!dir.exists(cicero_output_dir)) {
    dir.create(cicero_output_dir, recursive = TRUE)
}

# Create outfile name
links_outfile <- sprintf(sprintf("%s/cicero_connections.tsv", cicero_output_dir))

# Read in Seurat object
print("Reading in Seurat object!")
data = readRDS(seurat_object)

accessibility_data <- GetAssayData(object = data, assay = "peaks", slot = "data")
accessibility_data <- 1 * (accessibility_data > 0)

# Make CDS
cellinfo <- data.frame(cells = colnames(accessibility_data))
rownames(cellinfo) <- cellinfo$cells

peakinfo <- data.frame(peaks = rownames(accessibility_data))
rownames(peakinfo) <- peakinfo$peaks

input_cds <-  suppressWarnings(new_cell_data_set(accessibility_data,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
input_cds <- monocle3::detect_genes(input_cds)

# Preprocessing
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

write.table(conns, links_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
sprintf("Output to %s!", links_outfile)