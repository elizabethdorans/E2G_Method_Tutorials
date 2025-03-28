## Reformat pillar project predictors

library(optparse)

# Process input arguments --------------------------------------------------------------------------

# create arguments list
option_list = list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Path to transcripts input file", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "Path to output file", metavar = "character"),
  make_option(c("-g", "--genes_file"), type = "character", default = NULL,
              help = "Path to file containing gene information", metavar = "character"),
  make_option(c("-c", "--cell_type"), type = "character", default = NULL,
              help = "Cell type", metavar = "character"),
  make_option(c("-u", "--sample_summary_short"), type = "character", default = NULL,
              help = "Short summary of biosample for header line", metavar = "character"),
  make_option(c("-d", "--sample_term_id"), type = "character", default = NULL,
              help = "UBERON / CL for cellType/biosample", metavar = "character"),
  make_option(c("-m", "--method"), type = "character", default = NULL,
              help = "E2G method that produced the predictions", metavar = "character"),
  make_option(c("-v", "--version"), type = "character", default = NULL,
              help = "E2G method version", metavar = "character"),
  make_option(c("-e", "--metadata"), type = "character", default = NULL,
              help = "IGVF data portal accession", metavar = "character"),
  make_option(c("-s", "--score_column"), type = "character", default = NULL,
              help = "Column name containing main predictor score", metavar = "character"),  
  make_option(c("-t", "--score_type"), type = "character", default = NULL,
              help = "Type of the used score, e.g. positive_score", metavar = "character"),
  make_option(c("--threshold"), type = "character", default = NULL,
              help = "Used score threshold if applicable", metavar = "character") 
  
)

# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
  }
}

# check that all required parameters are provided
required_args <- c("input_file", "output_file", "method", "version", "score_type")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

# Process file -------------------------------------------------------------------------------------

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# load genes file
message("Loading predictions...")
genes <- fread(opt$genes_file, select = c("GeneSymbol", "GeneEnsemblID", "TSSEnd"))

# load input file
pred <- fread(opt$input_file)

# get all score columns (all columns except EG-pair defining columns)
message("Reformatting predictions...")

# Create chr, start, end columns from peak column
coord_cols = data.frame(str_split_fixed(pred$peak, "-", 3))
colnames(coord_cols) = c("ElementChr", "ElementStart", "ElementEnd")
pred <- cbind(pred, coord_cols)

# Restrict to autosomes + chromosome X
included_chromosomes = paste("chr", c(seq(22), "X"), sep = "")
pred = pred[pred$ElementChr %in% included_chromosomes,]

# Drop peak column
pred <- pred %>% rename(GeneSymbol = gene)
pred <- pred %>% select(-peak)

alt_score_cols <- setdiff(colnames(pred),
                          c("ElementChr", "ElementStart", "ElementEnd", "GeneSymbol", "CellType", opt$score_column))

# set cell type if specified
if (!is.null(opt$cell_type)) {
  pred$CellType <- opt$cell_type
}

# create header lines
header <- c(
  paste("# Source:", opt$method),
  paste("# Version:", opt$version),
  "# GenomeBuild: GRCh38",
  "# URL: [add url]",
  "# Assays: 10x multiome",
  "# SampleAgnostic: False",
  paste("# SampleTermName:", unique(pred$CellType)),
  paste("# SampleTermID:", opt$sample_term_id),
  paste("# SampleSummaryShort:", opt$sample_summary_short),
  paste("# ScoreType:", opt$score_type)
)

# add threshold if applicable
if (!is.null(opt$threshold)) {
  header <- c(header, paste("# ScoreThreshold:", opt$threshold))
}

header <- c(header, paste("# Metadata:", opt$metadata))

# add Ensembl id and gene TSS columns
pred <- left_join(pred, genes, by = "GeneSymbol")

# add additional columns and extract output columns
pred <- pred %>% 
  mutate(ElementName = paste0(ElementChr, ":", ElementStart, "-", ElementEnd),
         ElementClass = NA_character_) %>% 
  select(ElementChr, ElementStart, ElementEnd, ElementName, ElementClass,
         GeneSymbol, GeneEnsemblID, GeneTSS = TSSEnd,
         Score = all_of(opt$score_column), all_of(alt_score_cols))

# save to output file
message("Writing to output file...")
if (tools::file_ext(opt$output_file) == "gz") {
  
  # save to gzip compressed file
  tmp_file <- tools::file_path_sans_ext(opt$output_file)
  writeLines(header, con = tmp_file)
  fwrite(pred, file = tmp_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  system2("gzip", args = c("-f", tmp_file))
  
} else {
  
  # save to uncompressed file
  writeLines(header, con = opt$output_file)
  fwrite(pred, file = opt$output_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  
}

message("Done!")
