suppressPackageStartupMessages(library(SCENT))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(argparse))

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object with RNA, ATAC, and peak assays")
parser$add_argument("--candidate_link_file",
    help = "[REQUIRED] Path to file with candidate peak-gene links")
parser$add_argument("--scent_output_file",  default = "./scent.tsv",
    help = "Path to file for SCENT output")

args <- parser$parse_args()

seurat_object = args$seurat_object
candidate_link_file = args$candidate_link_file
scent_output_file = args$scent_output_file

# Create output directory if needed
print(sprintf("Output directory: %s", dirname(scent_output_file)))
if (!dir.exists(dirname(scent_output_file))) {
    dir.create(dirname(scent_output_file), recursive = TRUE)
}

# Read in Seurat object
print("Reading in Seurat object!")
data = readRDS(seurat_object)

# Extract expression (counts) and ATAC (binarized) matrices
rna_mtx = GetAssayData(data, assay = "RNA", layer = "counts")
atac_mtx = GetAssayData(data, assay = "peaks", layer = "counts")
atac_mtx <- 1 * (atac_mtx > 0)

# Add relevant metadata

# 1) nUMI
data@meta.data$nUMI = data@meta.data$nCount_RNA

# 2) Cell type (required by SCENT algorithm, but trivial after this preprocessing pipeline, since the Seurat object should already be restricted to cells of interest)
data@meta.data$celltype = "focal"

# 3) Percent mitochondrial reads
mito.genes <- grep(pattern = "^MT-", x = rownames(rna_mtx), value = TRUE)
percent.mito <- colSums(rna_mtx[mito.genes, ])/colSums(rna_mtx)
data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")

# Extract metadata
metadata = data@meta.data[, c("barcode", "nUMI", "percent.mito", "celltype")]
metadata$log_nUMI = log(metadata$nUMI)
metadata$cell = rownames(metadata)

# Read in candidate peak-gene links
print("Reading in candidate peak-gene links!")
candidate_pgl <- read.table(candidate_link_file)
colnames(candidate_pgl) <- c("gene", "peak")

# Run SCENT algorithm
SCENT_obj <- CreateSCENTObj(rna = rna_mtx, atac = atac_mtx, meta.data = metadata,
                            peak.info = candidate_pgl,
                            covariates = c("log_nUMI","percent.mito"), 
                            celltypes = "celltype")

print("Running SCENT!")
SCENT_obj <- SCENT_algorithm(object = SCENT_obj, celltype = "focal", ncores = 1)

write.table(SCENT_obj@SCENT.result, scent_output_file, sep = "\t", row.names = F, quote = F)
print(SCENT_obj@SCENT.result)