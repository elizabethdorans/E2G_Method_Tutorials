### Inputs
# Path to Seurat object with RNA + ATAC peak counts
# Path to ATAC fragment file (with .tbi file in same directory)
# Path to peaks bedfile
# Path to file with cell-level metadata

## Outputs
# ArchR project in /archr subfolder

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(plyranges))
addArchRGenome("hg38")

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="path to Seurat object with RNA, ATAC, and peak assays")
parser$add_argument("--fragment_file",
                    help = "path to ATAC fragment file")
parser$add_argument("--peaks_file",
                    help = "path to peaks bedfile")
parser$add_argument("--metadata_file", 
                    help = "path to file with cell-level metadata")
parser$add_argument("--result_dir",
                    help = "path to output file (will create directory if does not exist)")

args <- parser$parse_args()

seurat_object = args$seurat_object
fragment_file = args$fragment_file
peaks_file = args$peaks_file
metadata_file = args$metadata_file
result_dir = args$result_dir

data_name = basename(dirname(seurat_object))

inputFiles = c(fragment_file)
names(inputFiles) = c(data_name)
print(inputFiles)

archr_proj_dir = sprintf("%s/archr/project", result_dir)
print(archr_proj_dir)
if (!dir.exists(archr_proj_dir)) {
    print("Creating output directory")
    dir.create(archr_proj_dir, recursive = TRUE)
}

# Read in cells to consider
metadata = read.table(metadata_file, sep = "\t", header = TRUE)
cells = metadata$barcode
print(length(cells))

ArrowFiles <- sprintf("%s.arrow", data_name)
createArrowFiles(inputFiles, 
                 minTSS = 0,
                 minFrags = 0,
                 maxFrags = 1000000,
                 validBarcodes = cells,
                 addTileMat = TRUE,
                 addGeneScoreMat = FALSE,
                 force = TRUE)
proj <- ArchRProject(ArrowFiles)

# RNA data loading into Seurat object
seurat_obj <- readRDS(seurat_object)
                        
# Convert RNA data to SC Experiment and integrate into ArchR project
seRNA = as.SingleCellExperiment(seurat_obj, assay = "RNA")
logcounts(seRNA) = NULL

colnames(seRNA) = paste0(sprintf("%s#", data_name), colnames(seRNA))

gene_coords = getGenes(proj)
gene_coords = gene_coords[gene_coords$symbol %in% rownames(seRNA)]
gene_coords = sort(gene_coords, by=~symbol)

ridx <- rownames(seRNA) %in% gene_coords$symbol
seRNA = seRNA[ridx, , drop = FALSE]
seRNA = seRNA[order(row.names(seRNA)), ]

rowRanges(seRNA) <- gene_coords
rownames(seRNA) <- rowData(seRNA)$symbol

proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

proj <- addIterativeLSI(
    ArchRProj = proj, 
    clusterParams = list(resolution = c(2), 
                       sampleCells = 10000, 
                       maxClusters = 6, 
                       n.start= 10),
    saveIterations = FALSE,
    useMatrix = "TileMatrix", 
    depthCol = "nFrags",
    binarize = TRUE,
    name = "LSI_ATAC",
    force = TRUE
)

proj <- addIterativeLSI(
    ArchRProj = proj, 
    clusterParams = list(resolution = c(2), 
                       sampleCells = 10000, 
                       maxClusters = 6, 
                       n.start= 10),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA",
    force = TRUE
)

#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

peaks_granges = rtracklayer::import(peaks_file)
peaks_granges = peaks_granges %>% filter(seqnames != "chrM")
seqlevels(peaks_granges) <- seqlevelsInUse(peaks_granges)

proj <- addPeakSet(
  ArchRProj = proj,
  peakSet = peaks_granges,
  genomeAnnotation = getGenomeAnnotation(proj),
  force = TRUE
)

proj <- addPeakMatrix(
  ArchRProj = proj,
  binarize = TRUE,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addPeakMatrix")
)

saveArchRProject(ArchRProj = proj, outputDirectory = archr_proj_dir, load = FALSE)
