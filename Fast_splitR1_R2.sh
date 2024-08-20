#!/bin/bash
set -ue

# Get the bam file from the command line
DATA=$1
NAME=${2=sample}
THREADS=$3  # Set the number of threads to use

echo "Starting processing for forward strand..."

# Forward strand
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse strand
echo "Extracting forward strand alignments (first part)..."
samtools view -@ $THREADS -b -f 128 -F 16 $DATA > ${NAME}_fwd1.bam &
pid1=$!

echo "Extracting forward strand alignments (second part)..."
samtools view -@ $THREADS -b -f 80 $DATA > ${NAME}_fwd2.bam &
pid2=$!
wait $pid1 $pid2

echo "Merging forward strand alignments..."
samtools merge -@ $THREADS -f ${NAME}_fwd_unsorted.bam ${NAME}_fwd1.bam ${NAME}_fwd2.bam

echo "Sorting forward strand alignments..."
samtools sort -@ $THREADS -o ${NAME}_fwd.bam ${NAME}_fwd_unsorted.bam
samtools index -@ $THREADS ${NAME}_fwd.bam

echo "Forward strand processing complete."

echo "Starting processing for reverse strand..."

# Reverse strand
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
echo "Extracting reverse strand alignments (first part)..."
samtools view -@ $THREADS -b -f 144 $DATA > ${NAME}_rev1.bam &
pid3=$!

echo "Extracting reverse strand alignments (second part)..."
samtools view -@ $THREADS -b -f 64 -F 16 $DATA > ${NAME}_rev2.bam &
pid4=$!
wait $pid3 $pid4

echo "Merging reverse strand alignments..."
samtools merge -@ $THREADS -f ${NAME}_rev_unsorted.bam ${NAME}_rev1.bam ${NAME}_rev2.bam

echo "Sorting forward strand alignments..."
samtools sort -@ $THREADS -o ${NAME}_rev.bam ${NAME}_rev_unsorted.bam
samtools index -@ $THREADS ${NAME}_rev.bam

echo "Reverse strand processing complete."

# Cleanup intermediate files
echo "Cleaning up intermediate files..."
rm ${NAME}_fwd1.bam ${NAME}_fwd2.bam ${NAME}_rev1.bam ${NAME}_rev2.bam

echo "All processing complete."

