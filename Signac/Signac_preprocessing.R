suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(argparse))

set.seed(1234)

parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object with RNA, ATAC, and peak assays [.rds]")
parser$add_argument("--signac_output_dir", default = ".",
                    help = "Path to directory for output files")

args <- parser$parse_args()

seurat_object = args$seurat_object
signac_output_dir = args$signac_output_dir

# Create output directory if needed
print(sprintf("Output directory: %s", signac_output_dir))
if (!dir.exists(signac_output_dir)) {
    dir.create(signac_output_dir, recursive = TRUE)
}

outfile_name = sprintf("%s/seurat_object_signac_preprocessed.rds", signac_output_dir)

# Read in Seurat object
data <- readRDS(seurat_object)

# gene expression data processing
DefaultAssay(data) <- "RNA"
data <- SCTransform(data)
data <- RunPCA(data)

# DNA accessibility data processing
DefaultAssay(data) <- "peaks"
data <- FindTopFeatures(data, min.cutoff = 5)
data <- RunTFIDF(data)
data <- RunSVD(data)

# compute the GC content for each peak
data <- RegionStats(data, genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(data, outfile_name)