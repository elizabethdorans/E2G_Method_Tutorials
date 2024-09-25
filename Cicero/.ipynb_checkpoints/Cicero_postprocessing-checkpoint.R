suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser()

parser$add_argument("--input_file",
    help="[REQUIRED] Path to folder containing Signac peak-gene link predictions for each chromosome")
parser$add_argument("--promoter_peaks_bedfile",
    help="[REQUIRED] Path to bedfile of peaks noting overlap with promoter(s)")

args <- parser$parse_args()

input_file = args$input_file
promoter_peaks_bedfile = args$promoter_peaks_bedfile

input_folder = dirname(input_file)

# Read in bedfile specifying peaks overlapping promoters
promoter_peaks = read.table(promoter_peaks_bedfile, sep = "\t", header = TRUE)
colnames(promoter_peaks) = c("chr", "start", "end", "promoter")
promoter_peaks$peak = paste0("chr", paste(promoter_peaks$chr, (promoter_peaks$start + 1), promoter_peaks$end, sep = "-"))
promoter_peaks = promoter_peaks %>%
    filter(promoter != ".") %>%
    select(peak, gene = promoter)

# Read in Cicero peak-peak link predictions
cicero = read.table(input_file, sep = "\t", header = TRUE)
cicero = cicero[complete.cases(cicero),]

pgl = merge(cicero, promoter_peaks, by.x = "Peak1", by.y = "peak")
pgl = pgl %>% 
    mutate(peak = Peak2, Score = coaccess) %>%
    select(peak, gene, Score)

outfile = sprintf("%s/cicero_peak_gene_links.tsv", input_folder)
write.table(pgl, outfile, sep = "\t", row.names = FALSE, quote = FALSE)