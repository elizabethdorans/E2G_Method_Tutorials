input_file=$1
peaks_bedfile=$2

promoter_bedfile=../TSS_plusminus500bp.bed

peaks_bedfile_dir=$(echo $peaks_bedfile | rev | cut -d/ -f 2- | rev)

bedtools intersect -a <(cat $peaks_bedfile | cut -dr -f2- | cut -f1-3) -b $promoter_bedfile -loj | cut -f1,2,3,7 > $peaks_bedfile_dir/macs2_peaks_X_promoters.bed

Rscript E2G_Method_Tutorials/Cicero/Cicero_postprocessing.R --input_file $input_file --promoter_peaks_bedfile $peaks_bedfile_dir/macs2_peaks_X_promoters.bed