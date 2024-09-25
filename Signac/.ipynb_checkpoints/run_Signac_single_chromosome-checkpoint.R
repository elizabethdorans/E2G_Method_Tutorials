suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))

set.seed(1234)

parser <- ArgumentParser()

parser$add_argument("chromosome",
    help="[REQUIRED] chromosome whose genes are being tested, i.e. 1")
parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object preprocessed for Signac peak-gene linking (see Signac_preprocessing.R) [.rds]")
parser$add_argument("--signac_output_dir", default = ".",
    help = "Path to directory for output files")

args <- parser$parse_args()

chromosome = args$chromosome
seurat_object = args$seurat_object
signac_output_dir = args$signac_output_dir

# Create output directory if needed
print(sprintf("Output directory: %s", signac_output_dir))
if (!dir.exists(signac_output_dir)) {
    dir.create(signac_output_dir, recursive = TRUE)
}

# Create outfile name
links_outfile <- sprintf("%s/chr%s.tsv", signac_output_dir, chromosome)

# Read in Seurat object
print("Reading in Seurat object!")
data <- readRDS(seurat_object)

# Define set of genes on chromosome for peak-gene linking
annot <- Annotation(object = data[["peaks"]])
gene_coords = as.data.table(annot)[, c("seqnames", "gene_name")]
gene_coords = unique(gene_coords)
gene_set = gene_coords[seqnames == sprintf("chr%s", chromosome), ]$gene_name

sprintf("Linking across %s possible genes on chromosome %s!", length(gene_set), chromosome)

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
sprintf("Output to %s!", links_outfile)