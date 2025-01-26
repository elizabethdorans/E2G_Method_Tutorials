# Generate chunks of SCENT candidate peak-gene links from Seurat object

suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

set.seed(1234)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="path to Seurat object")
parser$add_argument("--number_candidate_links_per_file", 
    help="maximum number of candidate peak-gene links (one link per row) in each output file")
parser$add_argument("--candidate_link_file_output_dir", 
    help="path to folder where candidate peak-gene link files will be saved")
parser$add_argument("--tss_coordinates_file", 
    default = "../TSS.txt",
    help="path to file with TSS coordinates")

args <- parser$parse_args()

seurat_object = args$seurat_object
number_candidate_links_per_file = as.numeric(args$number_candidate_links_per_file)
candidate_link_file_output_dir = args$candidate_link_file_output_dir
tss_coordinates_file = args$tss_coordinates_file

# Create output directory if needed
print(sprintf("Output directory: %s", candidate_link_file_output_dir))
if (!dir.exists(candidate_link_file_output_dir)) {
    dir.create(candidate_link_file_output_dir, recursive = TRUE)
}

# Read in Seurat obj
data = readRDS(seurat_object)

# Read in gene universe
gene_universe = read.table(tss_coordinates_file, header = TRUE)

# Read in expression (counts) and ATAC (binarized) matrices
rna_mtx = GetAssayData(data, assay = "RNA")
atac_mtx = GetAssayData(data, assay = "peaks")
atac_mtx <- 1 * (atac_mtx > 0)

# Restrict to genes in gene universe expressed in > 5% of cells
rna_mtx = rna_mtx[intersect(gene_universe$gene, rownames(rna_mtx)),]
rna_mtx = rna_mtx[(rowSums(rna_mtx) / ncol(rna_mtx)) > 0.05,]

# Restrict to peaks active in > 5% of cells
atac_mtx = atac_mtx[(rowSums(atac_mtx) / ncol(atac_mtx)) > 0.05,]

print(sprintf("Subsetted to %s genes and %s peaks active in > 5 percent of cells", nrow(rna_mtx), nrow(atac_mtx)))

# Get peak coordinates
peaks = data.frame(peak = rownames(atac_mtx))
peaks$chr = sub("chr", "", sapply(strsplit(peaks$peak, "-"), "[[", 1))
peaks$center = (as.numeric(sapply(strsplit(peaks$peak, "-"), "[[", 2)) + as.numeric(sapply(strsplit(peaks$peak, "-"), "[[", 3))) / 2

# Get coordinates of windows +/- 500kb around gene TSS
gene_universe$tss_minus_500kb = gene_universe$tss - 500000
gene_universe$tss_plus_500kb = gene_universe$tss + 500000

# Define candidate peak-gene links
candidate_links = data.frame()

for (gene in rownames(rna_mtx)) {
    # Get coordinates of focal gene
    focal_gene_coords = gene_universe[gene_universe$gene == gene,]

    # Get peaks centered within 500kb of focal gene TSS
    focal_peaks = peaks[(peaks$chr == focal_gene_coords$chr) & 
                        (peaks$center < focal_gene_coords$tss_plus_500kb) &
                        (peaks$center > focal_gene_coords$tss_minus_500kb),]$peak
    
    if (length(focal_peaks) > 1) {
        # Define candidate links to the focal gene
        focal_candidate_links = data.frame(gene = focal_gene_coords$gene, peak = focal_peaks)

        # Add to candidate links dataframe
        candidate_links = rbind(candidate_links, focal_candidate_links)
    }
}

# Randomly split candidate peak-gene links into chunks
number_candidate_links <- nrow(candidate_links)
candidate_links_shuffled <- candidate_links[sample(1:number_candidate_links, number_candidate_links, replace = FALSE),]

chunk = 1
while (nrow(candidate_links_shuffled) > number_candidate_links_per_file) {
    outfile = sprintf("%s/chunk%s.txt", candidate_link_file_output_dir, chunk)
    write.table(tail(candidate_links_shuffled, n = number_candidate_links_per_file), outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    candidate_links_shuffled = head(candidate_links_shuffled, - number_candidate_links_per_file)
    chunk = chunk + 1
}

outfile = sprintf("%s/chunk%s.txt", candidate_link_file_output_dir, chunk)
write.table(candidate_links_shuffled, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print(sprintf("output %s chunks of candidate peak-gene links!", chunk))
