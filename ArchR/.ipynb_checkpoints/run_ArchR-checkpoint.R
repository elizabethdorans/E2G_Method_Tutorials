suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(plyranges))

addArchRGenome("hg38")

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object with RNA, ATAC, and peak assays")
parser$add_argument("--atac_fragments",
                    help = "[REQUIRED] ATAC fragment file (.tsv.gz) (tsv.gz.tbi file must be in same directory)")
parser$add_argument("--archr_output_dir",  default = ".",
    help = "Path to directory for output files")

args <- parser$parse_args()

seurat_object = args$seurat_object
atac_fragments = args$atac_fragments
archr_output_dir = args$archr_output_dir

# Create output directory if needed
print(sprintf("Output directory: %s", archr_output_dir))
if (!dir.exists(archr_output_dir)) {
    dir.create(archr_output_dir, recursive = TRUE)
}

# Create outfile name
links_outfile <- sprintf(sprintf("%s/archr_peak_gene_links.tsv", archr_output_dir))

# Load Seurat object
seurat_obj <- readRDS(seurat_object)

# Get focal cells from Seurat object
cells <- colnames(seurat_obj)

# Get peaks from Seurat object
peaks_granges = seurat_obj$peaks@ranges
peaks_granges = peaks_granges %>% filter(seqnames != "chrM")
seqlevels(peaks_granges) <- seqlevelsInUse(peaks_granges)

# Create ArchR project from fragment files
data_name = basename(dirname(seurat_object))
inputFiles = c(atac_fragments)
names(inputFiles) = c(data_name)
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

# Integrate RNA into ArchR project
seRNA = as.SingleCellExperiment(seurat_obj, assay = "RNA")
logcounts(seRNA) = NULL
colnames(seRNA) = paste0(sprintf("%s#", data_name), colnames(seRNA))

gene_coords = getGenes(proj)
proj_genes = gene_coords$symbol
RNA_genes = rownames(seRNA)

focal_genes = intersect(proj_genes, RNA_genes)

gene_coords = gene_coords[proj_genes %in% focal_genes]
gene_coords = sort(gene_coords, by=~symbol)

ridx = RNA_genes %in% focal_genes
seRNA = seRNA[ridx, , drop = FALSE]
seRNA = seRNA[order(row.names(seRNA)), ]

rowRanges(seRNA) <- gene_coords
rownames(seRNA) <- rowData(seRNA)$symbol

proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

# Run dimension reduction on RNA and ATAC
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

proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

# Add cell x peak matrix to ArchR project
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
  force = TRUE
)

# Run peak-gene linking
print(sprintf("Linking!"))

proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    maxDist = 500000,
    reducedDims = "LSI_Combined",
    useMatrix = "GeneExpressionMatrix",
    addEmpiricalPval = TRUE
)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = -1.5,
    FDRCutOff = 1,
    resolution = 1,
    returnLoops = FALSE
)

# Convert peak-gene links to readable output
peaks_idx = as.data.frame(metadata(p2g)[[1]])
peaks = paste0(peaks_idx$seqnames, "-", peaks_idx$start - 1, "-", peaks_idx$end)
gene_idx = as.data.frame(metadata(p2g)[[2]])
genes = gene_idx$name
p2g$peak = peaks[p2g$idxATAC]
p2g$gene = genes[p2g$idxRNA]
p2g$idxATAC = NULL
p2g$idxRNA = NULL
p2g$Score = p2g$Correlation

# Save project and peak-gene links
write.table(p2g, links_outfile, sep = "\t", row.names = FALSE, quote = FALSE)
saveArchRProject(ArchRProj = proj, outputDirectory = sprintf("%s/project", archr_output_dir), load = FALSE)
sprintf("Saved project to %s and peak-gene links to %s!", sprintf("%s/project", archr_output_dir), links_outfile)