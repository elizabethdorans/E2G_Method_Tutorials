seurat_object=$1
signac_output_dir=$2
max_peak_TSS_distance=$3

echo $max_peak_TSS_distance

# Autosomes
for chrom in {1..22}
do
    outfile=$signac_output_dir/chr${chrom}.tsv
    if ! [ -f $outfile ]
    then
        if [[ -z "$max_peak_TSS_distance" ]]
        then
            cmd="Rscript run_Signac_single_chromosome.R $chrom --seurat_object $seurat_object --signac_output_dir $signac_output_dir"
        else
            cmd="Rscript run_Signac_single_chromosome.R $chrom --seurat_object $seurat_object --signac_output_dir $signac_output_dir --max_peak_TSS_distance $max_peak_TSS_distance"
        fi
        echo $cmd
        sbatch --time=12:00:00 --mem=40G -p short -c 1 --wrap="$cmd" --output $signac_output_dir/signac_${chrom}.out
    fi
done

# X chromosome
outfile=$signac_output_dir/chrX.tsv
if ! [ -f $outfile ]
then
    if [[ -z "$max_peak_TSS_distance" ]]
    then
        cmd="Rscript run_Signac_single_chromosome.R X --seurat_object $seurat_object --signac_output_dir $signac_output_dir"
    else
        cmd="Rscript run_Signac_single_chromosome.R X --seurat_object $seurat_object --signac_output_dir $signac_output_dir --max_peak_TSS_distance $max_peak_TSS_distance"
    fi
    echo $cmd
    sbatch --time=12:00:00 --mem=40G -p short -c 1 --wrap="$cmd" --output $signac_output_dir/signac_X.out
fi
