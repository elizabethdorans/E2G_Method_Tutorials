#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

options <- list(
  make_option(c("-i", "--input_file"), type = "character", help = "MACS2 peak BED"),
  make_option(c("-o", "--output_file"), type = "character", help = "Portal element TSV or TSV.GZ"),
  make_option(c("-g", "--genes_file"), type = "character", help = "GENCODE gene bounds BED (CollapsedGeneBounds)"),
  make_option("--cell_metadata_file", type = "character", help = "Cell annotation report TSV"),
  make_option("--dataset_cell_type", type = "character", help = "Dataset/cell-type identifier using an underscore, for example igvf1_k562"),
  make_option("--peak_caller_version", type = "character", default = "2.2.7.1", help = "MACS2 version [default %default]"),
  make_option("--url", type = "character", default = "https://github.com/macs3-project/MACS", help = "MACS source URL [default %default]"),
  make_option("--genome_reference", type = "character", default = "IGVFDS0280IQAI", help = "IGVF genome reference [default %default]"),
  make_option("--assays", type = "character", default = "10x Multiome", help = "Assays header value [default %default]")
)
parser <- OptionParser(option_list = options)
opt <- parse_args(parser)

required <- c("input_file", "output_file", "genes_file", "cell_metadata_file",
              "dataset_cell_type", "peak_caller_version")
missing <- required[vapply(required, function(x) {
  is.null(opt[[x]]) || !nzchar(opt[[x]])
}, logical(1))]
if (length(missing)) {
  print_help(parser)
  stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

input_files <- c(opt$input_file, opt$genes_file, opt$cell_metadata_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Input file(s) not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
}
if (!grepl("^[^_]+_.+$", opt$dataset_cell_type)) {
  stop("--dataset_cell_type must use an underscore between dataset and cell type, for example igvf1_k562",
       call. = FALSE)
}

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(startsWith(x, "chr"), x, paste0("chr", x))
}

read_gene_universe <- function(path) {
  genes <- fread(path, header = TRUE, select = 1:8)
  setnames(genes, c("chr", "start", "end", "GeneSymbol", "score", "strand",
                    "GeneEnsemblID", "gene_type"))
  genes[, `:=`(
    chr = normalize_chr(chr),
    start = as.integer(start),
    end = as.integer(end),
    tss = fifelse(strand == "+", as.integer(start), as.integer(end))
  )]
  if (anyNA(genes$start) || anyNA(genes$end) || any(genes$end <= genes$start)) {
    stop("Gene universe contains invalid coordinates", call. = FALSE)
  }
  if (any(!genes$strand %in% c("+", "-"))) {
    stop("Gene universe contains invalid strands", call. = FALSE)
  }
  unique(genes, by = c("chr", "start", "end", "GeneSymbol",
                       "GeneEnsemblID", "strand"))
}

classify_elements <- function(elements, genes) {
  query <- elements[, .(
    chr = ElementChr,
    start = BedStart,
    overlap_end = ElementEnd - 1L,
    element_id
  )]
  gene_bodies <- genes[, .(chr, start, overlap_end = end - 1L)]
  promoters <- genes[, .(
    chr,
    start = pmax(0L, tss - 500L),
    overlap_end = tss + 499L
  )]
  setkey(gene_bodies, chr, start, overlap_end)
  setkey(promoters, chr, start, overlap_end)
  promoter_ids <- unique(foverlaps(
    query, promoters,
    by.x = c("chr", "start", "overlap_end"),
    by.y = c("chr", "start", "overlap_end"),
    type = "any", nomatch = 0L
  )$element_id)
  genic_ids <- unique(foverlaps(
    query, gene_bodies,
    by.x = c("chr", "start", "overlap_end"),
    by.y = c("chr", "start", "overlap_end"),
    type = "any", nomatch = 0L
  )$element_id)
  classes <- data.table(
    element_id = elements$element_id,
    ElementClass = "intergenic"
  )
  classes[element_id %in% genic_ids, ElementClass := "genic"]
  classes[element_id %in% promoter_ids, ElementClass := "promoter"]
  classes
}

write_portal_file <- function(x, header, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("igvf_elements_", tmpdir = dirname(path), fileext = ".tsv")
  on.exit(unlink(temporary), add = TRUE)
  writeLines(header, temporary)
  fwrite(x, temporary, sep = "\t", quote = FALSE, na = "NA",
         append = TRUE, col.names = TRUE)
  if (tools::file_ext(path) == "gz") {
    status <- system2("gzip", c("-c", shQuote(temporary)), stdout = path)
    if (!identical(status, 0L)) stop("gzip failed for ", path, call. = FALSE)
  } else if (!file.rename(temporary, path)) {
    stop("Failed to create ", path, call. = FALSE)
  }
}

message("Loading metadata from ", opt$cell_metadata_file)
metadata <- fread(opt$cell_metadata_file, na.strings = "NA")
metadata_columns <- c("Dataset_Cluster", "SampleTermName", "SampleTermID",
                      "CellAnnotation")
missing_metadata <- setdiff(metadata_columns, names(metadata))
if (length(missing_metadata)) {
  stop("Metadata file is missing: ", paste(missing_metadata, collapse = ", "),
       call. = FALSE)
}
metadata[, DatasetCellType := sub("^([^-]+)-", "\\1_", Dataset_Cluster)]
metadata_row <- metadata[DatasetCellType == opt$dataset_cell_type]
if (nrow(metadata_row) != 1) {
  stop("Expected one metadata row for '", opt$dataset_cell_type,
       "'; found ", nrow(metadata_row), call. = FALSE)
}
for (column in metadata_columns[-1]) {
  value <- metadata_row[[column]][1]
  if (is.na(value) || !nzchar(trimws(value))) {
    stop("Metadata field ", column, " is empty", call. = FALSE)
  }
  if (grepl(" | ", value, fixed = TRUE)) {
    stop("Metadata field ", column, " has unresolved multiple values",
         call. = FALSE)
  }
}
sample_term_name <- metadata_row$SampleTermName[1]
sample_term_id <- metadata_row$SampleTermID[1]
cell_annotation <- metadata_row$CellAnnotation[1]

message("Loading elements from ", opt$input_file)
elements <- fread(opt$input_file, header = FALSE, select = 1:3)
setnames(elements, c("ElementChr", "BedStart", "ElementEnd"))
elements[, `:=`(
  ElementChr = normalize_chr(ElementChr),
  BedStart = suppressWarnings(as.integer(BedStart)),
  ElementEnd = suppressWarnings(as.integer(ElementEnd))
)]
if (anyNA(elements$BedStart) || anyNA(elements$ElementEnd)) {
  stop("Peak coordinates must be integers", call. = FALSE)
}
if (any(elements$BedStart < 0L | elements$ElementEnd <= elements$BedStart)) {
  stop("Peak file contains invalid BED intervals", call. = FALSE)
}
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
excluded <- sum(!elements$ElementChr %in% standard_chromosomes)
if (excluded) {
  message("Removing ", excluded, " element(s) on nonstandard chromosomes")
  elements <- elements[ElementChr %in% standard_chromosomes]
}
elements <- unique(elements, by = c("ElementChr", "BedStart", "ElementEnd"))
if (!nrow(elements)) {
  stop("No elements remain after filtering", call. = FALSE)
}
elements[, element_id := .I]

message("Classifying elements")
genes <- read_gene_universe(opt$genes_file)
classes <- classify_elements(elements, genes)
elements <- merge(elements, classes, by = "element_id", all.x = TRUE, sort = FALSE)

# MACS2 BED starts are zero-based. Portal element coordinates and peak IDs from
# the E2G methods use one-based starts, as shown in the IGVF example file.
elements[, ElementStart := BedStart + 1L]
elements[, ElementName := paste0(ElementChr, ":", ElementStart, "-", ElementEnd)]
output_columns <- c("ElementChr", "ElementStart", "ElementEnd", "ElementName",
                    "ElementClass")
formatted <- elements[, ..output_columns]

if (anyNA(formatted) || any(!nzchar(formatted$ElementName))) {
  stop("Formatted output contains missing required values", call. = FALSE)
}
if (any(formatted$ElementStart < 1L | formatted$ElementEnd < formatted$ElementStart)) {
  stop("Formatted output contains invalid coordinates", call. = FALSE)
}
if (!all(formatted$ElementClass %in% c("promoter", "genic", "intergenic"))) {
  stop("Formatted output contains invalid ElementClass values", call. = FALSE)
}
expected_names <- paste0(formatted$ElementChr, ":", formatted$ElementStart,
                         "-", formatted$ElementEnd)
if (!identical(formatted$ElementName, expected_names)) {
  stop("ElementName does not match element coordinates", call. = FALSE)
}
if (anyDuplicated(formatted$ElementName)) {
  stop("Formatted output contains duplicated elements", call. = FALSE)
}

header <- c(
  "# Source: MACS2",
  paste("# Version:", opt$peak_caller_version),
  paste("# GenomeReference:", opt$genome_reference),
  paste("# URL:", opt$url),
  paste("# Assays:", opt$assays),
  "# SampleAgnostic: False",
  paste("# SampleTermName:", sample_term_name),
  paste("# SampleTermID:", sample_term_id),
  paste("# CellAnnotation:", cell_annotation)
)

message("Writing ", nrow(formatted), " elements to ", opt$output_file)
write_portal_file(formatted, header, opt$output_file)
message("Done")
