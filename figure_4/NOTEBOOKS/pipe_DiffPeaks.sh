#!/usr/bin/env bash

bams=$(tail -n +2 $1 | cut -f 4 | tr "\t" " ")
mark=$(tail -n +2 $1 | cut -f 1 | uniq | grep "^H3K\|^ATAC")
bed_peaks=$2
out_folder=$3

mkdir -p $out_folder/output_homer
echo "bams: "$bams
echo "mark: "$mark
echo "outdir: "$out_folder

multiBamSummary BED-file -b $bams -o $out_folder/$mark"_counts.mat.gz" --BED $bed_peaks -bl /home/ferrari/ferrari/my_repository/blacklist_ChIP-Seq/GRCm38_General_readAttractingRegions.UseThisOne_DKFZ.bed -p 20 --outRawCounts $out_folder/$mark"_counts.counts" -e --minMappingQuality 3

Rscript ./DESeq2_pipe.R $out_folder/$mark"_counts.counts" $1 $out_folder 

#Annotate Peaks Homer
source activate Homer

annotatePeaks.pl $out_folder/peaks_Annotation.bed mm10 -gtf /home/ferrari/ferrari/my_repository/annotations_gencode/mouse/M18/annotation_snakePipes/gencode.vM18.annotation.sorted.gtf > $out_folder/output_homer/$mark"_Homer_PeaksAnnotation.annot"

conda deactivate

#Annotate Peaks Uropa
source activate uropa

uropa -b $out_folder/peaks_Annotation.bed -g /home/ferrari/ferrari/my_repository/annotations_gencode/mouse/M18/annotation_snakePipes/gencode.vM18.annotation.sorted.gtf --feature transcript --distance 1500 --internals 1 -o $out_folder/output_uropa -t 10 --show_attributes gene_id transcript_id gene_name gene_type transcript_type #--filter_attribute gene_type --attribute_values protein_coding

conda deactivate

extract_from_gtf.py -o $out_folder -f TSS -w transcript \
/home/ferrari/ferrari/my_repository/annotations_gencode/mouse/M18/annotation_snakePipes/gencode.vM18.annotation.sorted.gtf

create_final_table.py $out_folder $mark


