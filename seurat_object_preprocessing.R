suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--rna_matrix", 
                    help = "[REQUIRED] File with cell x gene RNA count matrix [.mtx]")
parser$add_argument("--rna_matrix_barcodes",
                    help = "[REQUIRED] File with cell barcodes corresponding to row names of RNA count matrix [.txt]")
parser$add_argument("--rna_matrix_genes",
                    help = "[REQUIRED] File with gene names corresponding to column names RNA count matrix [.txt]")
parser$add_argument("--atac_fragments",
                    help = "[REQUIRED] ATAC fragment file (.tsv.gz) (tsv.gz.tbi file must be in same directory)")
parser$add_argument("--filtered_barcodes", 
<<<<<<< HEAD
                    help = "File with subset of cell barcodes (in column 'barcode') to include in Seurat object (can contain additional metadata columns) [.txt]")
parser$add_argument("--macs2_folder", default = "./macs2",
                    help = "Path to macs for calling peaks")
parser$add_argument("--output_dir",
=======
                    help = "[OPTIONAL] File with subset of cell barcodes (in column 'barcode') to include in Seurat object (can contain additional metadata columns) [.txt]")
parser$add_argument("--macs2_folder", default = "./macs2",
                    help = "Path to macs2 executable for peak calling")
parser$add_argument("--output_dir", default = ".",
>>>>>>> 3c75d3e7c9c74458e36c7aa5b095c0372ac29208
                    help = "Path to directory for output files")

args <- parser$parse_args()

rna_matrix = args$rna_matrix
rna_matrix_barcodes = args$rna_matrix_barcodes
rna_matrix_genes = args$rna_matrix_genes
atac_fragments = args$atac_fragments
filtered_barcodes = args$filtered_barcodes
macs2_folder = args$macs2_folder
output_dir = args$output_dir

# Create output directory if needed
print(sprintf("Output directory: %s", output_dir))
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# Read in RNA count matrix
counts <- Matrix::readMM(rna_matrix)
rownames(counts) <- read.table(rna_matrix_barcodes)$V1
colnames(counts) <- read.table(rna_matrix_genes)$V1

# Transpose to a gene x cell matrix
counts <- t(counts)

# Load RNA count matrix into Seurat object
obj <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)

# Add cell barcode as metadata column
obj@meta.data$barcode = rownames(obj@meta.data)

# Read in filtered cell barcodes and subset Seurat object (if applicable)
if (!is.null(filtered_barcodes)) {
    meta <- read.table(filtered_barcodes, sep= "\t", header = TRUE)
    rownames(meta) = meta$barcode
    
    if (length(meta) > 1) {
        meta <- meta[complete.cases(meta),]
    }
    
    # Subset Seurat object to filtered barcodes
    focal_cells = meta$barcode
    
    obj <- subset(obj, subset = barcode %in% focal_cells)

    # If file contains additional metadata columns, add metadata to Seurat object
    if (length(meta) > 1) {
        obj <- AddMetaData(
          object = obj,
          metadata = meta
        )
    }
} else {
    focal_cells = colnames(obj)
}

sprintf("%s cells retained in Seurat object", length(focal_cells))

# Create Seurat fragment object from ATAC fragment file
frags <- CreateFragmentObject(
  atac_fragments,
  cells = focal_cells
)

# Call and refine peaks (using MACS2)
peaks <- CallPeaks(frags, 
                   macs2.path = macs2_folder
                  )
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# Define peak-cell matrix
macs2_counts <- FeatureMatrix(
  fragments = frags,
  features = peaks,
  cells = focal_cells
)

# Integrate peaks into Seurat object
obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = atac_fragments,
  annotation = annotation
)

# Load ATAC data into Seurat object
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  fragments = atac_fragments,
  annotation = annotation
)
    
# Output peaks bedfile
export.bed(peaks, sprintf("%s/macs2_peaks.bed", output_dir))  

# Output Seurat object
seurat_object_outfile <- sprintf("%s/seurat_object.rds", output_dir)
saveRDS(obj, seurat_object_outfile)