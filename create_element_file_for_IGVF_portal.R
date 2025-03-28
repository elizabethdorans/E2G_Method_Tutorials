suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
#suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser()

parser$add_argument("--input_file",
    help="[REQUIRED] Path to bedfile with peaks called on dataset")
parser$add_argument("--peak_caller_version",
    help="version of macs2 used to call peaks")
parser$add_argument("--cell_type",
    help="version of macs2 used to call peaks")
parser$add_argument("--output_file", 
    default = "element_file.tsv",
    help="Path to for formatted element file")

args <- parser$parse_args()

input_file = args$input_file
output_file = args$output_file
peak_caller_version = args$peak_caller_version
cell_type = args$cell_type

# create header lines
header <- c(
  paste("# Source: macs2"),
  paste("# Version:", peak_caller_version),
  "# GenomeBuild: GRCh38",
  "# URL: [add url]",
  "# Assays: 10x multiome",
  "# SampleAgnostic: False",
  paste("# SampleTermName:", cell_type)
)

elements <- read.table(input_file, sep = "\t")[,1:3]
colnames(elements) <- c("ElementChr", "ElementStart", "ElementEnd")
elements$"ElementStart" <- as.integer(elements$"ElementStart" + 1)
elements$"ElementEnd" <- as.integer(elements$"ElementEnd")
elements <- elements %>% 
  mutate(ElementName = paste0(ElementChr, ":", ElementStart, "-", ElementEnd), ElementClass = NA_character_)


writeLines(header, output_file)
write.table(elements, file = output_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE, row.names = FALSE, col.names = TRUE)