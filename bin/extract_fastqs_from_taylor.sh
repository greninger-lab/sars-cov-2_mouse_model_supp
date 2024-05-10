#!/bin/bash

# get the current directory so we know where to find remove_soft_clips.R
bin_dir=$(pwd)

# define the directory where we'll find the taylor output bams with soft clips
# e.g. /Volumes/lizso_backup_drive/GSU_CL_VOCBeta_Swift/April2024/swift_out/
parent_dir=$1
# provide the full path to the Wuhan-Hu-1 reference fasta from the TAYLOR pipeline
# required to interpret the bam files...
#e.g. ~/Documents/covid_swift_pipeline/NC_045512.2.fasta
ref_fasta=$2

# go into the parent directory
cd $parent_dir

# make a directory to put the hard-clipped fastqs
mkdir taylor_output_fastqs/

for i in *.clipped.bam
do
base=$(basename $i .clipped.bam)
echo $base
samtools sort -o ${base}_sorted.bam ${base}.clipped.bam
samtools view -h -o to_clip.sam ${base}_sorted.bam
Rscript $bin_dir/remove_soft_clips.R $parent_dir
echo "Converting to Bam and adding header"
samtools view -T $ref_fasta -S -b hard_clipped.sam > ${base}_hard_clipped.bam 
echo "Extracting fastqs"
bedtools bamtofastq -i ${base}_hard_clipped.bam -fq ./taylor_output_fastqs/${base}.fastq 
gzip ./taylor_output_fastqs/*.fastq
rm *.sam
done


