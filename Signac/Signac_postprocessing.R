suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser()

parser$add_argument("--input_folder",
    help="[REQUIRED] Path to folder containing Signac peak-gene link predictions for each chromosome")
parser$add_argument("--gene_universe_file",
    help="path to gene universe file with GeneSymbol column")

args <- parser$parse_args()

input_folder = args$input_folder
gene_universe_file = args$gene_universe_file

# Read in peak-gene links per chromosome
in_files = Sys.glob(sprintf("%s/chr*.tsv", input_folder))
pgl = data.frame()
for (file in in_files) {
    chrom_pgl = read.table(file, sep = "\t", header = TRUE)
    pgl = rbind(pgl, chrom_pgl)
}

# Restrict peak-gene links to gene universe
if (!is.null(gene_universe_file)) {
    sprintf("Restricting to gene universe in %s!", gene_universe_file)
    gene_universe <- unique(read.table(gene_universe_file, header = TRUE)$GeneSymbol)
    pgl <- pgl[pgl$gene %in% gene_universe,]
}

# Rename score column
pgl = pgl %>% rename(Score = score)

# Drop duplicates
pgl = pgl[!duplicated(pgl), ]

outfile = sprintf("%s/signac_peak_gene_links.tsv", input_folder)
write.table(pgl, outfile, sep = "\t", row.names = FALSE, quote = FALSE)