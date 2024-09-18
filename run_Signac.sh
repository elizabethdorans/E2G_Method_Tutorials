seurat_obj=$1
cell_type_taskfile=$2
cell_types_by_barcode=$3

result_dir=$(echo $seurat_obj | rev | cut -d/ -f2- | rev)

if ! [ -d $result_dir/signac/pgl ]
then
    mkdir -p $result_dir/signac/pgl
fi

for chrom in {1..22}
do
    outfile=$result_dir/signac/pgl/chr${chrom}.tsv
    if ! [ -f $outfile ]
    then
        cmd="Rscript run_Signac_single_chromosome.R $chrom --seurat_obj $seurat_obj $result_dir/seurat_object_processed.rds"
        echo $cmd
        sbatch --time=8:00:00 --mem=20G --output=outfiles/run_signac.out --error=outfiles/run_signac.err -p short -c 1 --wrap="$cmd"
    fi