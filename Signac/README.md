### Signac Tutorial

Required inputs:
1) Seurat object: a Seurat object containing a matrix of raw RNA counts (cell x gene), ATAC peak counts (cell x peak)

Optional inputs:
1) Cell barcodes file: a text file containing one cell barcode per line (cell barcodes must be a subset of barcodes in the Seurat object

parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object with RNA, ATAC, and peak assays [.rds]")
parser$add_argument("--signac_output_dir", default = ".",
                    help = "Path to directory for output files")

args <- parser$parse_args()
