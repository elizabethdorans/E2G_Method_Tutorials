### Inputs
# Path to Seurat object with RNA + ATAC peak counts
# Name of metadata column to use for batch correction (optional)

## Outputs
# Processed Seurat object

suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(harmony))

set.seed(1234)

parser <- ArgumentParser()

parser$add_argument("--seurat_object",
    help="path to Seurat object with RNA, ATAC, and peak assays")
parser$add_argument("--batch_colname",
    help="Name of column to correct on using Harmony")
parser$add_argument("--outfile_name",
    help="Name of file to send processed/batch corrected Seurat object")

args <- parser$parse_args()

seurat_object = args$seurat_object
batch_colname = args$batch_colname
outfile_name = args$outfile_name

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

# Batch correction on RNA-Seq
if (!is.null(batch_colname)) {
    # Make factor out of batch metadata column
    data@meta.data[,batch_colname] = as.factor(data@meta.data[,batch_colname])
    data <- RunHarmony(data, batch_colname)
}

# compute the GC content for each peak
data <- RegionStats(data, genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(data, outfile_name)

# Signac input: save gene sets by chromosome
expression.data <- GetAssayData(object = data, assay = "RNA")
genecounts <- rowSums(x = expression.data > 0)
genes.keep <- genecounts > 0
expression.data <- expression.data[genes.keep, , drop = FALSE]
genes <- rownames(x = expression.data)

annot <- Annotation(object = data[["peaks"]])

gene_chr = as.data.table(annot)[, c("seqnames", "gene_name")]
gene_chr = unique(gene_chr)
gene_chr = gene_chr[!(seqnames %in% c("chrX", "chrY", "chrM")), ]

out_dir = dirname(seurat_obj)
gene_list_dir = sprintf("%s/gene_sets_by_chr", out_dir)
dir.create(gene_list_dir, recursive = TRUE)

for (chr in unique(gene_chr$seqnames)) {
    gene_lst = gene_chr[seqnames == chr, ]$gene_name
    write.table(gene_lst, sprintf("%s/%s.txt", gene_list_dir, chr), 
                quote = F, row.names = F, col.names = F)
    }