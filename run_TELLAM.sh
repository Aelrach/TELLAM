#!/bin/bash

# Get the directory where the script is located
script_dir=$(dirname "$(pwd "${BASH_SOURCE[0]}")")

# Function to display the help message
function show_help {
    echo "Usage: $(basename $0) -d <deseq2_table> -state <bam_state> -chr <chr_prefix> -name <condition_name> -control <control_name> -elements <element_pattern> [options]"
    echo
    echo "Required arguments:"
    echo "  -d          Path to the DESeq2 table (e.g., results.csv)"
    echo "  -state      State of bam files, 'r' for raw bam files or 'f' for strand separated AND exon-less filtered bam"
    echo "  -chr        Prefix used to identify chromosomes in bam files, 'chr' or 'None'"
    echo "  -name       Name to refer to the condition (e.g., 'AZA')"
    echo "  -control    Name to refer to the control (e.g., 'DMSO')"
    echo "  -elements   List of patterns found in the IDs of the annotation file used to filter rows of deseq2 table (e.g., 'LTR')"
    echo
    echo "Optional arguments:"
    echo "  -raw        Must be set if -state is 'r', Path to the folder containing all raw bam files"
    echo "  -exons      Must be set if -state is 'r', Path to the bed file containing all exons of genome. VERIFY that chromosome naming of this file and that of your bam files match"
    echo "  -threads    Number of threads to use if -state is set to 'r' (default is 8)"
    echo "  -fb         Must be set if -state is 'f', Path to the folder containing all forward and reverse bam files"
    echo "  -directory  Project directory to save files (default is \${name}_TELLAM)"
    echo "  -annotation Path to the annotation file (default is 'TELLAM/GRCh38_Ensembl_rmsk_TE.gtf.locInd.locations')"
    echo "  -consensus  Path to the consensus fasta file (default is 'TELLAM/UCSC_TE_consensus.fa.txt')"
    echo "  -window     Window size in base pairs used to compute 3' and 5' context (default is 3000bp)"
    echo "  -context    [0,1], context ratio threshold below which exponential decrease of the context effect starts. (default is 0.95)"
    echo "  -size       [0,1], size ratio threshold below which exponential decrease of the size effect starts. (default is 0.1)"
    echo "  -coverage   [0,1], Should be very low values, Reads per base coverage threshold below which exponential decrease of the coverage effect starts. (default is 0.06)"
    echo "  -full       [0,1], Proportion of consensus size above which the loci is considered full_length (default is 0.9)"
    echo
    echo "Example:"
    echo "  $(basename $0) -d DESEQ.txt -state f -fb path_to_filtered_bam -chr chr -name AZA -control DMSO -elements 'L1:LINE',LTR "
    exit 1
}

# Parse command-line arguments
# Default values
threads="${threads:-8}"
directory="${directory:-$script_dir/${condition_name}vs${control_name}_TELLAM}"
annotation="${annotation:-$script_dir/GRCh38_Ensembl_rmsk_TE.gtf.locInd.locations}"
consensus="${consensus:-$script_dir/UCSC_TE_consensus.fa.txt}"
window="${window:-3000}"
context="${context:-0.95}"
size="${size:-0.1}"
coverage="${coverage:-0.06}"
full="${full:-0.9}"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d) deseq2_table="$2"; shift ;;
        -state) bam_state="$2"; shift ;;
        -chr) chr_prefix="$2"; shift ;;
        -name) condition_name="$2"; shift ;;
        -control) control_name="$2"; shift ;;
        -elements) element_pattern="$2"; shift ;;
        -raw) raw_bam_folder="$2"; shift ;;
        -exons) exons_bed="$2"; shift ;;
        -threads) threads="$2"; shift ;;
        -fb) filtered_bam_folder="$2"; shift ;;
        -directory) directory="$2"; shift ;;
        -annotation) annotation="$2"; shift ;;
        -consensus) consensus="$2"; shift ;;
        -window) window="$2"; shift ;;
        -context) context="$2"; shift ;;
        -size) size="$2"; shift ;;
        -coverage) coverage="$2"; shift ;;
        -h|--help) show_help ;;
        *) echo "Unknown parameter passed: $1"; show_help ;;
    esac
    shift
done

# Check for required arguments
if [ -z "$deseq2_table" ] || [ -z "$bam_state" ] || [ -z "$chr_prefix" ] || [ -z "$condition_name" ] || [ -z "$control_name" ] || [ -z "$element_pattern" ] || [ -z "$annotation" ]; then
    echo "Error: Missing required arguments."
    show_help
fi

if [ "$bam_state" = "r" ] && [ -z "$raw_bam_folder" ]; then
     echo "Error: -state is set to 'r' but -raw not defined. Pass to -raw argument the path to raw bam files folder"
     show_help
fi

if [ "$bam_state" = "r" ] && [ -z "$exons_bed" ]; then
     echo "Error: -state is set to 'r' but -exons not defined. Pass to -exons argument the path to exons bed file"
     show_help
fi

if [ "$bam_state" = "f" ] && [ -z "$filtered_bam_folder" ]; then
     echo "Error: -state is set to 'f' but -fb not defined. Pass to -fb argument the path to all filtered bam files folder (forward and reverse bam files must be in the same directory)"
     show_help
fi

# Create the project directory if it doesn't exist
mkdir -p "$directory"

if [ "$bam_state" = "r" ]; then
    # Make Exon-less bam files
    bash "$script_dir/filterOutExonsBAM.sh" "$raw_bam_folder" "$exons_bed" "$directory" "$threads"
    wait $!
    
    filtered_bam_folder="$directory/filtered_bam"
    
    # Make Forward and Reverse bam files
    bash "$script_dir/launchFastSplit.sh" "$script_dir" "$directory" "$filtered_bam_folder" "$threads"
    wait $!
fi
# Create the configuration table for the Python script
config_file="$directory/config.txt"
cat > "$config_file" <<EOL
TELLAM_DIR=$script_dir
DESEQ2_TABLE=$deseq2_table
BAM_STATE=$bam_state
CHR_PREFIX=$chr_prefix
CONDITION_NAME=$condition_name
CONTROL_NAME=$control_name
ELEMENT_PATTERN=$element_pattern
ANNOTATION=$annotation
RAW_FOLDER=$raw_bam_folder
EXONS=$exons_bed
THREADS=$threads
FILTERED_FOLDER=$filtered_bam_folder
DIRECTORY=$directory
CONSENSUS=$consensus
WINDOW=$window
CONTEXT=$context
SIZE=$size
COVERAGE=$coverage
FULL=$full
EOL

echo "Configuration file created at $config_file"
echo "Launching the Python pipeline..."

# Launch the Python script with the configuration file
python=$(command -v python)
$python "$script_dir/TELLAM.py" "$config_file"

echo "Pipeline completed."
