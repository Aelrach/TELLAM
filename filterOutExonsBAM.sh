#!/bin/bash
sample_folder=$1
exons_bed=$2
directory=$3
THREADS=$4

SampleList=$(find "$sample_folder" -type f -name "*.bam")

L=${#SampleList[@]}

cd directory

mkdir -p filtered_bam

for ((i=0;i<=L-1;i++))
do
    #Filter out exonic reads using bedtools intersect
    nohup bedtools intersect -v -abam $sample_folder/sample${i}_fwd.bam -b Database/chrPREFIX_exons_h38.bed > filtered_bam/filtered_fwd_sample${i}.bam &
    nohup bedtools intersect -v -abam $sample_folder/sample${i}_rev.bam -b exons_bed > filtered_bam/filtered_rev_sample${i}.bam &
done

for file in filtered_bam/*.bam; 
do
    nohup samtools -@ $THREADS index "$file" &
done
