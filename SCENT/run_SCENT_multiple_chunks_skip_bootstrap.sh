seurat_object=$1
candidate_link_folder=$2
scent_output_folder=$3
script=$4

# Save name of R script from script directory if no alternate path provided
if [ -z "${script}" ]
then
    script=run_SCENT.R
fi

for candidate_link_file in $candidate_link_folder/*
do
    # Name output file
    basename=$(echo $candidate_link_file | rev | cut -d/ -f1 | cut -d. -f2- | rev)
    scent_output_file=$scent_output_folder/$basename.tsv
    
    if ! [ -f $scent_output_file ]
    then
        # Run SCENT
        cmd="Rscript $script --seurat_object $seurat_object --candidate_link_file $candidate_link_file --scent_output_file $scent_output_file --skip_bootstrap"
        echo $cmd
        sbatch --time=1-0:00:00 --mem=10G -p medium -c 1 --wrap="$cmd"
    fi
done
