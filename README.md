# TELLAM (Transposable Elements Locus Level Analysis Metric)
Pipeline for Transposable Elements (TEs) Locus Level Analysis. This pipeline was created to analyze results derived from RNAseq experiments. It leverages differential expression analysis and genomic context to analyze TE activity by using DESeQ2 tables, exon-less and strand separated bam files to compute a metric, TELLAM, on each annotated loci given in input. It then determines which are activated or not using a small Random Forest Classifier.

# Requirements
- bedtools
- samtools
- python 3.10.14 tested
  - biopython 1.81 tested
  - numpy 1.25.0 tested
  - pandas 2.0.3 tested
  - pysam 0.22.0 tested
  - swifter 1.3.4 tested
  - scikit-learn 1.5.0 tested
 
# Installation

### 1. Clone the Repository

```bash
git clone https://github.com/Aelrach/TELLAM.git
cd TELLAM
```
### 2. Create a Virtual Environment

It's recommended to use a virtual environment to manage dependencies:

```bash
conda create --name TELLAM python=3.10
conda activate TELLAM
pip install -r requirements.txt 
```
Nice step-by-step tutorial on how to get miniconda up and running from the command line : https://medium.com/@rajiitmandi21/installing-miniconda-and-setting-up-conda-environments-on-linux-0c114e76a324

### 3. Verify that you have these auxilliary files
Before using TELLAM, make sure you have these 2 things : a **locus-level annotation file** of transposable elements insertions in the genome of interest and a **FASTA file** containing **consensus sequences** for each element you wish to analyze.
- If you are analyzing human or mouse samples, we recommend that you download the locus-level annotation file provided by the **TElocal team at (https://www.mghlab.org/software/telocal)**
- The **FASTA file with consensus sequences is provided within TELLAM**. It was downloaded from UCSC's Database (https://hgwdev.gi.ucsc.edu/~max/kznf/hg38reps/hg38/seqs/). **You can provide your own.**

# Usage
Usage: 
```bash
bash run_TELLAM.sh -d <deseq2_table> -state <bam_state> -chr <chr_prefix> -name <condition_name> -control <control_name> -elements <element_pattern> [options]
```

### Required arguments:
- `-d`          Path to the DESeq2 table (e.g., `results.csv`)
- `-state`      State of BAM files: `'r'` for raw BAM files or `'f'` for strand-separated AND exon-less filtered BAM files
- `-chr`        Prefix used to identify chromosomes in BAM files: `'chr'` or `'None'`
- `-name`       Name to refer to the condition (e.g., `'AZA'`)
- `-control`    Name to refer to the control (e.g., `'DMSO'`)
- `-elements`   List of patterns found in the IDs of the annotation file used to filter rows of DESeq2 table (e.g., `'LTR'`)
- `-raw`        IF AND ONLY IF `-state` is `'r'`; Path to the folder containing all TREATED CONDITION ONLY raw BAM files
- `-exons`      IF AND ONLY IF `-state` is `'r'`; Path to the BED file containing all exons of the genome. Verify that the chromosome naming of this file matches that of your BAM files
- `-fb`         IF AND ONLY IF `-state` is `'f'`; Path to the folder containing all TREATED CONDITION ONLY forward and reverse BAM files
  
### Optional arguments:
- `-threads`    Number of threads to use (default is 8)
- `-directory`  Project directory to save files (default is `${name}_TELLAM`)
- `-consensus`  Path to the consensus FASTA file (default is `'TELLAM/consensus.fasta'`)

### Analysis optional parameters:
- `-window`     Window size in base pairs used to compute 3' and 5' context (default is 3000 bp)
- `-context`    Between (0,1), context ratio threshold below which exponential decrease of the context effect starts (default is 0.95)
- `-size`       Between (0,1), size ratio threshold below which exponential decrease of the size effect starts (default is 0.1)
- `-coverage`   Between (0,1); should be very low values; reads per base coverage threshold below which exponential decrease of the coverage effect starts (default is 0.06)
- `-full`       Between (0,1), proportion of consensus size above which the locus is considered full_length (default is 0.9)

Example:
```bash
      bash run_TELLAM.sh -d DESEQ_AZAvsDMSO.txt -state f -fb path_to_AZA_filtered_bam -chr chr -name AZA -control DMSO -elements 'L1:LINE',LTR
```

In this command the user is specifying that exon-less, strand separated bam files (both forward and reverse in the same folder) for the **TREATED CONDITION**, are **all** located in path_to_filtered_bam. 
The command also specifies that the bam files were generated using chr1, chr2 etc.. chromosome naming and to only analyze loci which are L1:LINE or LTR.  

Note that to specify several patterns you must give each pattern in a comma separated format (No spaces allowed between patterns). 
**e.g VALID : L1:LINE,LTR,AluJo / NOT VALID : L1:LINE, LTR, AluJo**

Beware that the rows of your deseq table MUST have as ID the name of the loci. To this end, we recommend users to use the loci annotation file provided by TElocal (https://www.mghlab.org/software/telocal) when making their DESEQ2 tables. If you have a specially curated set of loci, make sure that every locus' ID follows the convention : 

"family"_dup"copy_name or number":"family":"Type":"Class"  

**e.g AluJo_dup71176:AluJo:Alu:SINE**

# Using TELLAM with raw BAM files (-state set to 'r')
If you plan to use TELLAM with raw BAM files, you **MUST** specify a **BED file containing all exons** of your genome of interest in the **-elements argument**. Before running the actual TELLAM pipeline, the script will produce exon-less, strand separated bam files. This will take significantely more time to run.

# Output
The pipeline will produce a processed version of your DESeQ2 table as an intermediary file but the main output is a Table in bed format containing the following columns : 

"#chr", "start", "end", "TE", "family", "position", "strand", "score", "3v5_effect", "MeanCoverage", "size_ratio", "full_length", "size_effect", "Metric", "Activated"

This table contains all informations needed to analyze the activation of your loci, providing positions, names and TELLAM features for each loci in bed format. A subset of this table containing only loci considered as reactivated by TELLAM will be produced. All files are saved in either the directory specified by  -directory or in an automatically generated folder in the TELLAM folder.
