#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

options <- list(
  make_option(c("-i", "--input_file"), type = "character", help = "Input peak-gene prediction TSV"),
  make_option(c("-o", "--output_file"), type = "character", help = "Unthresholded output TSV or TSV.GZ"),
  make_option(c("-g", "--genes_file"), type = "character", help = "GENCODE gene bounds BED (CollapsedGeneBounds)"),
  make_option("--element_file", type = "character", help = "Formatted portal element TSV or TSV.GZ"),
  make_option("--cell_metadata_file", type = "character", help = "Cell annotation report TSV"),
  make_option("--dataset_cell_type", type = "character", help = "Dataset/cell-type identifier using an underscore, for example igvf1_k562"),
  make_option(c("-m", "--method"), type = "character", help = "E2G method"),
  make_option(c("-v", "--version"), type = "character", help = "E2G method version"),
  make_option(c("-s", "--score_column"), type = "character", help = "Input column to publish as Score"),
  make_option(c("-t", "--score_type"), type = "character", help = "IGVF ScoreType"),
  make_option("--url", type = "character", default = NULL, help = "Method URL; known methods are filled automatically"),
  make_option("--genome_reference", type = "character", default = "IGVFDS0280IQAI", help = "IGVF genome reference [default %default]"),
  make_option("--assays", type = "character", default = "10x Multiome", help = "Assays header value [default %default]"),
  make_option("--optional_columns", type = "character", default = NULL, help = "Comma-separated input columns to retain"),
  make_option("--threshold_column", type = "character", default = "Score", help = "Column used for thresholding [default %default]"),
  make_option("--threshold_value", type = "double", default = NULL, help = "Optional Score cutoff; comparison is determined by --score_type"),
  make_option("--thresholded_output_file", type = "character", default = NULL, help = "Optional thresholded output TSV or TSV.GZ")
)
parser <- OptionParser(option_list = options)
opt <- parse_args(parser)

required <- c("input_file", "output_file", "genes_file", "element_file", "cell_metadata_file",
              "dataset_cell_type", "method", "version", "score_column", "score_type")
missing <- required[vapply(required, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
if (length(missing)) {
  print_help(parser)
  stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

input_files <- c(opt$input_file, opt$cell_metadata_file, opt$genes_file,
                 opt$element_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop("Input file(s) not found: ", paste(missing_files, collapse = ", "), call. = FALSE)

score_types <- c("positive_score", "negative_score", "p_value", "adj_p_value", "divergent", "boolean")
if (!opt$score_type %in% score_types) stop("Invalid --score_type: ", opt$score_type, call. = FALSE)

method_urls <- c(
  cicero = "https://github.com/cole-trapnell-lab/cicero-release",
  pgboost = "https://github.com/elizabethdorans/pgBoost/tree/IGVF_E2G_Pillar_Project",
  scent = "https://github.com/immunogenomics/SCENT",
  signac = "https://github.com/stuart-lab/signac"
)
if (is.null(opt$url)) opt$url <- unname(method_urls[tolower(opt$method)])
if (length(opt$url) != 1 || is.na(opt$url) || !nzchar(opt$url)) {
  stop("No known URL for method '", opt$method, "'; supply --url", call. = FALSE)
}

threshold_args <- c("threshold_value", "thresholded_output_file")
threshold_supplied <- vapply(threshold_args, function(x) !is.null(opt[[x]]), logical(1))
if (any(threshold_supplied) && !all(threshold_supplied)) {
  stop("Thresholding requires both --threshold_value and --thresholded_output_file", call. = FALSE)
}
make_thresholded <- all(threshold_supplied)
thresholdable_score_types <- c("positive_score", "divergent", "p_value")
if (make_thresholded && !opt$score_type %in% thresholdable_score_types) {
  stop("Thresholding is supported for ScoreType positive_score, divergent, or p_value", call. = FALSE)
}
if (make_thresholded && (!is.finite(opt$threshold_value) || opt$threshold_value < 0)) {
  stop("--threshold_value must be a finite non-negative number", call. = FALSE)
}
if (make_thresholded && normalizePath(opt$output_file, mustWork = FALSE) ==
    normalizePath(opt$thresholded_output_file, mustWork = FALSE)) {
  stop("Unthresholded and thresholded output paths must differ", call. = FALSE)
}

metadata <- fread(opt$cell_metadata_file, na.strings = "NA")
metadata_columns <- c("Dataset_Cluster", "SampleTermName", "SampleTermID", "CellAnnotation")
missing_metadata <- setdiff(metadata_columns, names(metadata))
if (length(missing_metadata)) stop("Metadata file is missing: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
if (!grepl("^[^_]+_.+$", opt$dataset_cell_type)) {
  stop("--dataset_cell_type must use an underscore between dataset and cell type, for example igvf1_k562", call. = FALSE)
}
# Kayla's report calls this field Dataset_Cluster and separates the dataset
# token from the cluster with a hyphen. Normalize it at the input boundary so
# this pipeline consistently uses the dataset_cell_type underscore convention.
metadata[, DatasetCellType := sub("^([^-]+)-", "\\1_", Dataset_Cluster)]
metadata_row <- metadata[DatasetCellType == opt$dataset_cell_type]
if (nrow(metadata_row) != 1) {
  stop("Expected one metadata row for '", opt$dataset_cell_type, "'; found ", nrow(metadata_row), call. = FALSE)
}
for (column in metadata_columns[-1]) {
  value <- metadata_row[[column]][1]
  if (is.na(value) || !nzchar(trimws(value))) stop("Metadata field ", column, " is empty", call. = FALSE)
  if (grepl(" | ", value, fixed = TRUE)) stop("Metadata field ", column, " has unresolved multiple values", call. = FALSE)
}
sample_term_name <- metadata_row$SampleTermName[1]
sample_term_id <- metadata_row$SampleTermID[1]
cell_annotation <- metadata_row$CellAnnotation[1]

read_gene_universe <- function(path) {
  genes <- fread(path, header = TRUE, select = c(2, 3, 4, 6, 7))
  setnames(genes, c("start", "end", "GeneSymbol", "strand", "GeneEnsemblID"))
  genes[, `:=`(
    start = as.integer(start),
    end = as.integer(end),
    GeneTSS = fifelse(strand == "+", as.integer(start) + 1L,
                      as.integer(end) + 1L)
  )]
  if (anyNA(genes$start) || anyNA(genes$end) || any(genes$end <= genes$start)) {
    stop("Gene universe contains invalid coordinates", call. = FALSE)
  }
  if (any(!genes$strand %in% c("+", "-"))) {
    stop("Gene universe contains invalid strands", call. = FALSE)
  }
  genes <- unique(genes, by = c("GeneSymbol", "GeneEnsemblID", "GeneTSS"))
  ambiguous <- genes[, .N, by = GeneSymbol][N > 1]
  if (nrow(ambiguous)) {
    stop("Ambiguous gene symbols in annotation, including: ",
         paste(head(ambiguous$GeneSymbol, 10), collapse = ", "), call. = FALSE)
  }
  genes
}

read_element_header <- function(path) {
  con <- if (tools::file_ext(path) == "gz") gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, n = 50L, warn = FALSE)
  lines[startsWith(lines, "#")]
}

header_value <- function(header, field) {
  prefix <- paste0("# ", field, ":")
  line <- header[startsWith(header, prefix)]
  if (length(line) != 1L) {
    stop("Element file must contain exactly one ", prefix, " header", call. = FALSE)
  }
  trimws(sub(prefix, "", line, fixed = TRUE))
}

read_elements <- function(path) {
  header <- read_element_header(path)
  expected_header <- c(
    GenomeReference = opt$genome_reference,
    SampleTermName = sample_term_name,
    SampleTermID = sample_term_id,
    CellAnnotation = cell_annotation
  )
  for (field in names(expected_header)) {
    observed <- header_value(header, field)
    if (!identical(observed, unname(expected_header[field]))) {
      stop("Element file ", field, " is '", observed, "'; expected '",
           expected_header[field], "'", call. = FALSE)
    }
  }

  elements <- fread(path, skip = "ElementChr", header = TRUE,
                    na.strings = c("NA", ""))
  required_element_columns <- c("ElementChr", "ElementStart", "ElementEnd",
                                "ElementName", "ElementClass")
  if (!identical(names(elements), required_element_columns)) {
    stop("Element file must contain exactly: ",
         paste(required_element_columns, collapse = ", "), call. = FALSE)
  }
  missing_counts <- vapply(elements, function(x) {
    sum(is.na(x) | (is.character(x) & !nzchar(x)))
  }, integer(1))
  if (any(missing_counts)) {
    stop("Element file has missing required values: ",
         paste(names(missing_counts)[missing_counts > 0],
               missing_counts[missing_counts > 0], sep = "=",
               collapse = ", "), call. = FALSE)
  }
  if (any(elements$ElementStart < 1L |
          elements$ElementEnd < elements$ElementStart)) {
    stop("Element file has invalid coordinates", call. = FALSE)
  }
  expected_names <- paste0(elements$ElementChr, ":", elements$ElementStart,
                           "-", elements$ElementEnd)
  if (!identical(elements$ElementName, expected_names)) {
    stop("Element file has ElementName values inconsistent with its coordinates",
         call. = FALSE)
  }
  if (!all(elements$ElementClass %in% c("promoter", "genic", "intergenic"))) {
    stop("Element file has invalid ElementClass values", call. = FALSE)
  }
  if (anyDuplicated(elements$ElementName)) {
    stop("Element file has duplicated ElementName values", call. = FALSE)
  }
  elements
}

required_output <- c("ElementChr", "ElementStart", "ElementEnd", "ElementName",
                     "ElementClass", "GeneSymbol", "GeneEnsemblID", "GeneTSS",
                     "CellAnnotation", "Score")
validate_output <- function(x, label) {
  if (!identical(names(x)[seq_along(required_output)], required_output)) stop(label, " has invalid column order", call. = FALSE)
  missing_counts <- vapply(x[, ..required_output], function(z) sum(is.na(z) | (is.character(z) & !nzchar(z))), integer(1))
  if (any(missing_counts)) stop(label, " has missing required values: ",
    paste(names(missing_counts)[missing_counts > 0], missing_counts[missing_counts > 0], sep = "=", collapse = ", "), call. = FALSE)
  if (any(x$ElementStart < 0L | x$ElementEnd <= x$ElementStart)) stop(label, " has invalid coordinates", call. = FALSE)
  if (!all(x$ElementClass %in% c("promoter", "genic", "intergenic"))) stop(label, " has invalid ElementClass", call. = FALSE)
  if (!all(is.finite(x$Score))) stop(label, " has non-finite Score", call. = FALSE)
  if (!all(x$CellAnnotation == cell_annotation)) stop(label, " has inconsistent CellAnnotation", call. = FALSE)
}

build_header <- function(threshold = NULL) {
  header <- c(
    paste("# Source:", opt$method), paste("# Version:", opt$version),
    paste("# GenomeReference:", opt$genome_reference), paste("# URL:", opt$url),
    paste("# Assays:", opt$assays), "# SampleAgnostic: False",
    paste("# SampleTermName:", sample_term_name), paste("# SampleTermID:", sample_term_id),
    paste("# CellAnnotation:", cell_annotation)
  )
  if (!is.null(threshold)) header <- c(header, paste("# ScoreThreshold:", threshold))
  c(header, paste("# ScoreType:", opt$score_type))
}

write_portal_file <- function(x, header, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("igvf_e2g_", tmpdir = dirname(path), fileext = ".tsv")
  on.exit(unlink(temporary), add = TRUE)
  writeLines(header, temporary)
  fwrite(x, temporary, sep = "\t", quote = FALSE, na = "NA", append = TRUE, col.names = TRUE)
  if (tools::file_ext(path) == "gz") {
    status <- system2("gzip", c("-c", shQuote(temporary)), stdout = path)
    if (!identical(status, 0L)) stop("gzip failed for ", path, call. = FALSE)
  } else if (!file.rename(temporary, path)) {
    stop("Failed to create ", path, call. = FALSE)
  }
}

message("Loading predictions from ", opt$input_file)
pred <- fread(opt$input_file, na.strings = c("NA", ""))
has_peak_gene <- all(c("peak", "gene") %in% names(pred))
has_portal_names <- all(c("ElementName", "GeneSymbol") %in% names(pred))
if (has_peak_gene == has_portal_names) {
  stop("Prediction file must contain either peak/gene or ElementName/GeneSymbol columns",
       call. = FALSE)
}
element_input_column <- if (has_peak_gene) "peak" else "ElementName"
gene_input_column <- if (has_peak_gene) "gene" else "GeneSymbol"
optional_columns <- if (is.null(opt$optional_columns) || !nzchar(opt$optional_columns)) character() else trimws(strsplit(opt$optional_columns, ",", fixed = TRUE)[[1]])
if (any(!nzchar(optional_columns)) || anyDuplicated(optional_columns)) stop("Invalid --optional_columns", call. = FALSE)
if (length(intersect(optional_columns, required_output))) {
  stop("--optional_columns cannot repeat required output columns", call. = FALSE)
}
threshold_input_column <- if (opt$threshold_column == "Score") {
  opt$score_column
} else {
  opt$threshold_column
}
needed <- unique(c(
  element_input_column, gene_input_column, opt$score_column, optional_columns,
  if (make_thresholded) threshold_input_column else character()
))
missing_input <- setdiff(needed, names(pred))
if (length(missing_input)) stop("Prediction file is missing: ", paste(missing_input, collapse = ", "), call. = FALSE)

if (has_peak_gene) {
  parts <- tstrsplit(pred$peak, "-", fixed = TRUE)
  if (length(parts) != 3) stop("Peaks must have format chrN-start-end", call. = FALSE)
  peak_start <- suppressWarnings(as.integer(parts[[2]]))
  peak_end <- suppressWarnings(as.integer(parts[[3]]))
  if (anyNA(peak_start) || anyNA(peak_end) ||
      any(peak_start < 1L | peak_end < peak_start)) {
    stop("Peak coordinates must be valid one-based integers", call. = FALSE)
  }
  pred[, ElementName := paste0(parts[[1]], ":", peak_start, "-", peak_end)]
} else {
  pred[, ElementName := as.character(ElementName)]
  if (anyNA(pred$ElementName) || any(!nzchar(pred$ElementName))) {
    stop("ElementName contains missing or empty values", call. = FALSE)
  }
}

message("Matching predictions to ", opt$element_file)
elements <- read_elements(opt$element_file)
missing_elements <- unique(pred[!ElementName %chin% elements$ElementName,
                                ElementName])
if (length(missing_elements)) {
  stop(length(missing_elements),
       " linked element(s) are absent from --element_file, including: ",
       paste(head(missing_elements, 10), collapse = ", "), call. = FALSE)
}
input_rows <- nrow(pred)
pred <- merge(pred, elements, by = "ElementName", all.x = TRUE, sort = FALSE)
if (nrow(pred) != input_rows) stop("Element join changed row count", call. = FALSE)

pred[, (opt$score_column) := suppressWarnings(as.numeric(get(opt$score_column)))]
if (anyNA(pred[[opt$score_column]]) || any(!is.finite(pred[[opt$score_column]]))) {
  stop("Main score column contains missing or non-numeric values", call. = FALSE)
}
message("Annotating genes")
genes <- read_gene_universe(opt$genes_file)
if (has_peak_gene) {
  pred[, GeneSymbol := gene]
} else {
  pred[, GeneSymbol := as.character(GeneSymbol)]
}
outside_gene_universe <- unique(pred[!GeneSymbol %chin% genes$GeneSymbol, GeneSymbol])
if (length(outside_gene_universe)) {
  remove_rows <- !pred$GeneSymbol %chin% genes$GeneSymbol
  message(
    "Removing ", sum(remove_rows), " prediction(s) for ",
    length(outside_gene_universe),
    " gene symbol(s) outside the GENCODE v43 gene universe, including: ",
    paste(head(outside_gene_universe, 10), collapse = ", ")
  )
  pred <- pred[!remove_rows]
}
input_rows <- nrow(pred)
if (input_rows == 0) {
  stop("No predictions remain after restricting to the gene universe", call. = FALSE)
}
pred <- merge(pred, genes[, .(GeneSymbol, GeneEnsemblID, GeneTSS)],
              by = "GeneSymbol", all.x = TRUE, sort = FALSE)
if (nrow(pred) != input_rows) stop("Gene join changed row count", call. = FALSE)

pred[, `:=`(CellAnnotation = cell_annotation, Score = get(opt$score_column))]
output_columns <- c(required_output, optional_columns)
formatted <- pred[, ..output_columns]
validate_output(formatted, "Unthresholded output")
message("Writing unthresholded predictions to ", opt$output_file)
write_portal_file(formatted, build_header(), opt$output_file)

if (make_thresholded) {
  threshold_values <- if (opt$threshold_column == "Score") {
    formatted$Score
  } else {
    suppressWarnings(as.numeric(pred[[threshold_input_column]]))
  }
  if (anyNA(threshold_values) || any(!is.finite(threshold_values))) {
    stop("Threshold column '", opt$threshold_column,
         "' contains missing or non-numeric values", call. = FALSE)
  }
  keep <- switch(opt$score_type,
    positive_score = threshold_values >= opt$threshold_value,
    divergent = abs(threshold_values) >= opt$threshold_value,
    p_value = threshold_values < opt$threshold_value)
  thresholded <- formatted[keep]
  validate_output(thresholded, "Thresholded output")
  threshold_expression <- switch(opt$score_type,
    positive_score = paste(opt$threshold_column, ">="),
    divergent = paste0("abs(", opt$threshold_column, ") >="),
    p_value = paste(opt$threshold_column, "<"))
  threshold_text <- paste(threshold_expression,
                          format(opt$threshold_value, scientific = FALSE, trim = TRUE))
  message("Writing ", nrow(thresholded), " thresholded predictions to ", opt$thresholded_output_file)
  write_portal_file(thresholded, build_header(threshold_text), opt$thresholded_output_file)
}

message("Done: ", nrow(formatted), " unthresholded predictions")
