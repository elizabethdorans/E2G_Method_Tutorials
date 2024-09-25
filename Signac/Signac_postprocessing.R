suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser()

parser$add_argument("--input_folder",
    help="Path to folder containing Signac peak-gene link predictions for each chromosome")

args <- parser$parse_args()

input_folder = args$input_folder

in_files = Sys.glob(sprintf("%s/chr*.tsv", input_folder))

pgl = data.frame()
for (file in in_files) {
    chrom_pgl = read.table(file, sep = "\t", header = TRUE)
    pgl = rbind(pgl, chrom_pgl)
}

pgl = pgl %>% rename(Score = score)

outfile = sprintf("%s/signac_peak_gene_links.tsv", input_folder)
write.table(pgl, outfile, sep = "\t", row.names = FALSE, quote = FALSE)