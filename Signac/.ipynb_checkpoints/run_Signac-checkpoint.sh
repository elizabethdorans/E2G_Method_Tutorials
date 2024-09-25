seurat_object=$1
signac_output_dir=$2

# Autosomes
for chrom in {1..22}
do
    outfile=$signac_output_dir/chr${chrom}.tsv
    if ! [ -f $outfile ]
    then
        cmd="Rscript run_Signac_single_chromosome.R $chrom --seurat_object $seurat_object --signac_output_dir $signac_output_dir"
        echo $cmd
        sbatch --time=8:00:00 --mem=20G --output=outfiles/run_signac.out --error=outfiles/run_signac.err -p short -c 1 --wrap="$cmd"
    fi
done

# X chromosome
outfile=$signac_output_dir/chrX.tsv
if ! [ -f $outfile ]
then
    cmd="Rscript run_Signac_single_chromosome.R X --seurat_object $seurat_object --signac_output_dir $signac_output_dir"
    echo $cmd
    sbatch --time=8:00:00 --mem=20G --output=outfiles/run_signac.out --error=outfiles/run_signac.err -p short -c 1 --wrap="$cmd"
fi