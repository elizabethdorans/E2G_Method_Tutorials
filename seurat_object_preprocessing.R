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

# Resolve bundled executables relative to this script, not the caller's working
# directory (which can differ for submitted jobs).
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
    script_dir <- getwd()
}

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--rna_matrix", 
                    help = "[REQUIRED] Gene x cell or cell x gene RNA count matrix [.mtx]")
parser$add_argument("--rna_matrix_barcodes",
                    help = "[REQUIRED] File with cell barcodes for the RNA matrix [.txt]")
parser$add_argument("--rna_matrix_genes",
                    help = "[REQUIRED] File with gene names for the RNA matrix [.txt]")
parser$add_argument("--atac_fragments",
                    help = "[REQUIRED] ATAC fragment file (.tsv.gz) (tsv.gz.tbi file must be in same directory)")
parser$add_argument("--filtered_barcodes", 
                    help = "File with subset of cell barcodes (in column 'barcode') to include in Seurat object (can contain additional metadata columns) [.txt]")
parser$add_argument("--peaks_file",
                    help = "Path to file with pre-called peaks (first 3 columns correspond to chromosome, start, end coordinates)")
parser$add_argument("--rows_are_genes",
                    action = "store_true", help = "Supply for a gene x cell matrix when separate gene and barcode files are not provided.")
parser$add_argument("--macs2_folder", default = file.path(script_dir, "macs2"),
                    help = "Path to the macs2 executable used for calling peaks")
parser$add_argument("--output_dir",
                    help = "Path to directory for output files")

args <- parser$parse_args()

rna_matrix = args$rna_matrix
rna_matrix_barcodes = args$rna_matrix_barcodes
rna_matrix_genes = args$rna_matrix_genes
atac_fragments = args$atac_fragments
filtered_barcodes = args$filtered_barcodes
peaks_file = args$peaks_file
macs2_folder = args$macs2_folder
output_dir = args$output_dir
rows_are_genes = args$rows_are_genes

# Create output directory if needed
print(sprintf("Output directory: %s", output_dir))
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# Check if RNA matrix is a MatrixMarket file
con <- file(rna_matrix, "r")
first_line <- readLines(con, n = 1)
close(con)

# Check if the first line starts with "%%MatrixMarket"
if (grepl("^%%MatrixMarket", first_line)) {
    # Read in RNA count matrix
    counts <- Matrix::readMM(rna_matrix)
} else {
    counts <- read.table(rna_matrix, sep = ',', header = TRUE, row.names = 1)
}

# Add names supplied in separate feature and barcode files. Infer the matrix
# orientation when possible because Matrix Market files do not encode it.
if (!is.null(rna_matrix_barcodes) && !is.null(rna_matrix_genes)) {
    genes <- read.table(rna_matrix_genes, header = FALSE, stringsAsFactors = FALSE)$V1
    barcodes <- read.table(rna_matrix_barcodes, header = FALSE, stringsAsFactors = FALSE)$V1

    genes_by_cells <- nrow(counts) == length(genes) &
        ncol(counts) == length(barcodes)
    cells_by_genes <- nrow(counts) == length(barcodes) &
        ncol(counts) == length(genes)

    if (genes_by_cells) {
        message("RNA matrix orientation: genes x cells")
    } else if (cells_by_genes) {
        message("RNA matrix orientation: cells x genes; transposing")
        counts <- t(counts)
    } else {
        stop(sprintf(
            paste0("RNA matrix dimensions (%s x %s) do not match ",
                   "%s genes and %s barcodes"),
            nrow(counts), ncol(counts), length(genes), length(barcodes)
        ))
    }

    rownames(counts) <- genes
    colnames(counts) <- barcodes
} else if (!rows_are_genes) {
    # Without companion name files, preserve the original command-line behavior.
    counts <- t(counts)
}

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
if (!is.null(peaks_file)) {
    print(sprintf("Reading peaks from %s", peaks_file))
    
    # Read in peak coordinates
    peaks_df <- read.table(peaks_file,
                     header = TRUE)
    peaks_df <- peaks_df[, c(1:3)]
    peaks_df <- unique(peaks_df)
    names(peaks_df) <- c('chr','start','end')
    
    # Transform to GRanges object
    peaks <- with(peaks_df, GRanges(chr, IRanges(start, end)))
} else {
    print(sprintf("Calling peaks using macs2 (%s)", macs2_folder))
    if (!file.exists(macs2_folder) || file.access(macs2_folder, mode = 1) != 0) {
        stop(sprintf(
            "macs2 executable not found or not executable: %s",
            macs2_folder
        ))
    }

    peaks <- CallPeaks(frags,
                       macs2.path = macs2_folder
                      )
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
}

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
