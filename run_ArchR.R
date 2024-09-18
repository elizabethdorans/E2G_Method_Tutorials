### Inputs
# Path to ArchR project
# Path to Seurat object with RNA + ATAC peak counts

## Outputs
# .tsv file with ArchR peak-gene linking scores

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))
addArchRGenome("hg38")

set.seed(1234)

parser <- ArgumentParser()

parser$add_argument("--archr_project_path",
                    help = "path to directory containing archr project")
parser$add_argument("--cell_prefix",
                    help = "label added as prefix to all cell names in ArchR (base of RNA data directory)")
parser$add_argument("--cell_type_taskfile",
    help="Focal cell types (data will be subsetted to only these cells for linking)")
parser$add_argument("--cts_by_barcode",
    help="file with cell barcodes and types")

args <- parser$parse_args()

archr_project_path = args$archr_project_path
cell_prefix = args$cell_prefix
cell_type_taskfile = args$cell_type_taskfile
cts_by_barcode = args$cts_by_barcode

out_dir = sprintf("%s/pgl", dirname(archr_project_path))
print(out_dir)
if (!dir.exists(out_dir)) {
    print("Creating")
    dir.create(out_dir)
}
outfile = sprintf("%s/peak_gene_links_500kdist.tsv", out_dir)
print(outfile)

# Load ArchR project
proj <- loadArchRProject(path = archr_project_path)
print(length(getCellNames(proj)))

if (!file.exists(outfile)) {
    # Subset by cell type (if applicable)
    if (!is.null(cell_type_taskfile)) {
        cell_type_label = strsplit(tail(strsplit(cell_type_taskfile, "/")[[1]], 1), "[.]")[[1]][1]
        cell_types = paste(readLines(cell_type_taskfile))
        cell_types = gsub(" ", "-", cell_types)
        print(cell_type_label)
        print(cell_types)
        print(outfile)
        print(sprintf("originally %s cells", nCells(proj)))
        cell_idents = read.table(cts_by_barcode, sep = "\t", comment.char = "")
        colnames(cell_idents) = c("cell_type", "cell_barcode")
        cell_type_barcodes = subset(cell_idents, cell_type %in% cell_types)$cell_barcode
        cell_type_names = paste0(cell_prefix, cell_type_barcodes)
        print(head(cell_type_names))
        archr_proj_names = getCellNames(proj)
        print(head(archr_proj_names))
        cell_subset_names = intersect(cell_type_names, archr_proj_names)
        outfile = sprintf("%s/%s_%s", dirname(outfile), cell_types_label, basename(outfile))
    } else {
        cell_subset_names = getCellNames(proj)
    }
    print(sprintf("will link across %s cells", length(cell_subset_names)))    

    proj <- addPeak2GeneLinks(
        ArchRProj = proj,
        maxDist = 500000,
        reducedDims = "LSI_Combined",
        useMatrix = "GeneExpressionMatrix",
        cellsToUse = cell_subset_names,
        addEmpiricalPval = TRUE
    )
        
    p2g <- getPeak2GeneLinks(
        ArchRProj = proj,
        corCutOff = -1.5,
        FDRCutOff = 1,
        resolution = 1,
        returnLoops = FALSE
    )
        
    peaks_idx = as.data.frame(metadata(p2g)[[1]])
    peaks = paste0(peaks_idx$seqnames, "-", peaks_idx$start - 1, "-", peaks_idx$end)
    gene_idx = as.data.frame(metadata(p2g)[[2]])
    genes = gene_idx$name
    p2g$peak = peaks[p2g$idxATAC]
    p2g$gene = genes[p2g$idxRNA]
        
    write.table(p2g, file = outfile, sep = "\t", quote = FALSE)
        
}