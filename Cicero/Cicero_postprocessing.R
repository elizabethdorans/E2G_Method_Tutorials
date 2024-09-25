suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser()

parser$add_argument("--input_file",
    help="Path to folder containing Signac peak-gene link predictions for each chromosome")
parser$add_argument("--promoter_peaks_bedfile",
    help="Path to bedfile of peaks noting overlap with promoter(s)")

args <- parser$parse_args()

input_file = args$input_file
promoter_peaks_bedfile = args$promoter_peaks_bedfile

# Read in bedfile specifying peaks overlapping promoters
promoter_peaks = read.table(promoter_peaks_bedfile, sep = "\t", header = TRUE)
colnames(promoter_peaks) = c("chr", "start", "end", "promoter")
promoter_peaks$peak = promoter_peaks$chr + "-" + (promoter_peaks$start + 1) + "-" + promoter_peaks$end

# Read in Cicero peak-peak link predictions
cicero = read.table(input_file, sep = "\t", header = TRUE)
cicero = cicero[complete.cases(cicero),]

pgl = cicero.merge(promoter_peaks, by.x = "Peak1", by.y = "peak", all.x = TRUE)

pgl = cicero.merge(promoter_peaks, by.x = "Peak2", by.y = "peak", all.x = TRUE)

# How many links in each processed dataframe? Need to search both peak 1 and 2?