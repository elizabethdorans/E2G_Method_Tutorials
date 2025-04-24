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
parser$add_argument("--sample_summary_short",
              help = "Short summary of biosample for header line")
parser$add_argument("--sample_term_id",
              help = "UBERON / CL for cellType/biosample")
parser$add_argument("--metadata",
              help = "IGVF data portal accession")

args <- parser$parse_args()

input_file = args$input_file
output_file = args$output_file
peak_caller_version = args$peak_caller_version
cell_type = args$cell_type
sample_term_id = args$sample_term_id
sample_summary_short = args$sample_summary_short
metadata = args$metadata

# create header lines
header <- c(
  paste("# Source: macs2"),
  paste("# Version:", peak_caller_version),
  "# GenomeBuild: GRCh38",
  "# URL: [add url]",
  "# Assays: 10x multiome",
  "# SampleAgnostic: False",
  paste("# SampleTermName:", cell_type),
  paste("# SampleTermID:", sample_term_id),
  paste("# SampleSummaryShort:", sample_summary_short),
  paste("# Metadata:", metadata)
)

elements <- read.table(input_file, sep = "\t")[,1:3]

colnames(elements) <- c("ElementChr", "ElementStart", "ElementEnd")
elements$"ElementStart" <- as.integer(elements$"ElementStart" + 1)
elements$"ElementEnd" <- as.integer(elements$"ElementEnd")

print(head(elements))

# Restrict to autosomes + chromosome X
included_chromosomes = paste("chr", c(seq(22), "X"), sep = "")
elements = elements[elements$ElementChr %in% included_chromosomes,]

elements <- elements %>% 
  mutate(ElementName = paste0(ElementChr, ":", ElementStart, "-", ElementEnd), ElementClass = NA_character_)

print(head(elements))

# save to output file
message("Writing to output file...")
if (tools::file_ext(output_file) == "gz") {
  
  # save to gzip compressed file
  tmp_file <- tools::file_path_sans_ext(output_file)
  writeLines(header, con = tmp_file)
  fwrite(elements, file = tmp_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  system2("gzip", args = c("-f", tmp_file))
  
} else {
  
  # save to uncompressed file
  writeLines(header, output_file)
write.table(elements, file = output_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE, row.names = FALSE, col.names = TRUE)
  
}