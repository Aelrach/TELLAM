#!/bin/bash
script_dir=$1
directory=$2
filtered_folder=$3
THREADS=$4

sample_list=$(find "$filtered_folder" -type f -name "*.bam")

L=${#sample_list[@]}

cd filtered_folder

for ((i=0;i<=L-1;i++))
do
        nohup bash $script_dir/Fast_splitR1_R2.sh ${filtered_folder}/${sample_list[i]} sample${i} $THREADS > sample${i}_fastSplit.out &
done