### Inputs
# Chromosome for peak-gene linking
# Path to Seurat object with RNA + ATAC peak counts

## Outputs
# .tsv file with Signac peak-gene linking scores

suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))

set.seed(1234)

parser <- ArgumentParser()

parser$add_argument("chromosome",
    help="chromosome whose genes are being tested, i.e. 1")
parser$add_argument("--seurat_object",
    help="Focal cell type (data will be subsetted to only these cells for linking)")
parser$add_argument("--cell_type_taskfile",
    help="Focal cell type (data will be subsetted to only these cells for linking)")
parser$add_argument("--cell_types_by_barcode",
    help="file with cell barcodes and cell types")

args <- parser$parse_args()

chr = args$chromosome
seurat_object = args$seurat_object
cell_type_taskfile = args$cell_type_taskfile
cell_types_by_barcode <- args$cell_types_by_barcode

data_folder = dirname(seurat_object)
data_name = basename(data_folder)

# Define cell types and file label
if (!is.null(cell_type_taskfile)) {
    cell_types = paste(readLines(cell_type_taskfile))
    cell_types = gsub(" ", "-", cell_types)
    cell_types_label = strsplit(tail(strsplit(cell_type_taskfile, "/")[[1]], 1), "[.]")[[1]][1]
    print(cell_types_label)
    print(cell_types)
    data_folder = data_folder = sprintf("%s/%s", data_folder, cell_types_label)
}

out_dir <- sprintf("%s/signac/pgl", data_folder)
print(sprintf("Output directory %s", out_dir))
if (!dir.exists(out_dir)) {
    print("Creating directory")
    dir.create(out_dir, recursive = TRUE)
}

# Create outfile name
links_outfile <- sprintf("%s/chr%s.tsv", out_dir, chr)
print(links_outfile)

# Read in Seurat object
print("Reading in Seurat object!")
data <- readRDS(seurat_object)

# Subset data by cell type (if applicable)
if (!is.null(cell_type_taskfile)) {
    cell_idents = read.table(cell_types_by_barcode, sep = "\t")
    colnames(cell_idents) = c("cell_type", "cell_barcode")
    focal_barcodes = subset(cell_idents, cell_type %in% cell_types)$cell_barcode
    data <- subset(data, barcode %in% focal_barcodes)
    print(sprintf("subsetted to %s cells", 
                  nrow(data@meta.data))
         )
}

# Read in set of genes to test
print("Getting set of genes to test!")
gene_set_folder <- sprintf("%s/gene_sets_by_chr", data_folder)
gene_set = read.table(sprintf("%s/chr%s.txt", gene_set_folder, chr))[,1]
print(head(gene_set))
print(length(gene_set))

print("linking!")

DefaultAssay(data) <- "peaks"

data <- LinkPeaks(
  object = data,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = gene_set,
  pvalue_cutoff = 1,
  score_cutoff = 0
)

out_df = as.data.frame(Links(data))[,c("peak", "gene", "score", "pvalue")]
write.table(out_df, links_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
