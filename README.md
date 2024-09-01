# TELLAM (Transposable Elements Locus Level Analysis Metric)
Pipeline for Transposable Elements Locus Level Analysis. Created to analyze results derived from RNAseq experiments. It uses DESeQ2 tables, exon-les and strand separated bam files to compute a metric, TELLAM, on each annotated loci given in input. It then determines which are activated or not using a small Random Forest Classifier.

# Requirements
- bedtools
- samtools
- python
- The following python libraries :
  - numpy
  - pandas
  - pickle
  - biopython
  - pysam

# Usage
Usage:: 
```bash
bash run_TELLAM.sh -d <deseq2_table> -state <bam_state> -chr <chr_prefix> -name <condition_name> -control <control_name> -elements <element_pattern> -annotation <annotation> [options]
```

Required arguments:
      -d          Path to the DESeq2 table (e.g., results.csv)
      -state      State of bam files, 'r' for raw bam files or 'f' for strand separated AND exon-less filtered bam
      -chr        Prefix used to identify chromosomes in bam files, 'chr' or 'None'
      -name       Name to refer to the condition (e.g., 'AZA')
      -control    Name to refer to the control (e.g., 'DMSO')
      -elements   List of patterns found in the IDs of the annotation file used to filter rows of deseq2 table (e.g., 'LTR')
      -annotation Path to the annotation file (default is 'TELLAM/annotation.gtf')
Optional arguments:
      -raw        Must be set if -state is 'r', Path to the folder containing all raw bam files
      -exons      Must be set if -state is 'r', Path to the bed file containing all exons of genome. VERIFY that chromosome naming of this file and that of your bam files match
      -threads    Number of threads to use if -state is set to 'r' (default is 8)
      -fb         Must be set if -state is 'f', Path to the folder containing all forward and reverse bam files
      -directory  Project directory to save files (default is \${name}_TELLAM)
      -consensus  Path to the consensus fasta file (default is 'TELLAM/consensus.fasta')
*Analysis optional parameters*
      -window     Window size in base pairs used to compute 3' and 5' context (default is 3000bp)
      -context    between (0,1), context ratio threshold below which exponential decrease of the context effect starts. (default is 0.95)
      -size       between (0,1), size ratio threshold below which exponential decrease of the size effect starts. (default is 0.1)
      -coverage   between (0,1), Should be very low values, Reads per base coverage threshold below which exponential decrease of the coverage effect starts. (default is 0.06)
      -full       between (0,1), Proportion of consensus size above which the loci is considered full_length (default is 0.9)
Example:
```bash
      bash run_TELLAM.sh -d DESEQ.txt -state f -fb path_to_filtered_bam -chr chr -name AZA -control DMSO -elements 'L1:LINE',LTR -annotation path_to_TE_locus_annotation"
```
In this command the user is specifying that exon-less, strand separated bam files (both forward and reverse in the same folder) are located in path_to_filtered_bam. Specifies that the bam files were generated using chr1, chr2 etc.. chromosome naming and to only analyze loci which are L1:LINE or LTR. Beware that the rows of your deseq table MUST have as ID the name of the loci. To this end, we recommand users to use the loci annotation file provided by TElocal (https://www.mghlab.org/software/telocal) when making their DESEQ2 tables. If you have a specially curated set of loci, make sure that every locus' ID is giving the Type, sub-Type and Family.

# Output
The pipeline will produce a processed version of your DESeQ2 table as an intermediary file but the main output is a Table in bed format containing the following columns : 
"#chr", "start", "end", "TE", "family", "position", "strand", "score", "3v5_effect", "MeanCoverage", "size_ratio", "full_length", "size_effect", "Metric", "Activated"

This table contains all informations needed to analyze the activation of your loci, providing positions, names and TELLAM features for each loci in bed format. A subset of this table containing only activated loci will be produced. All files are saved in either the directory specified by -directory or in an automatically generated folder in the TELLAM folder.




