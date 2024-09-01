#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 09:48:49 2024

@author: Raphaël RUBRICE

Locus-Level analysis Pipeline
"""
import pandas as pd
import numpy as np
import pickle as pkl
from Bio import SeqIO
import pysam
import time
import re
import os

def load_config(config_file):
    """
    Loads configuration settings from a file into a dictionary.

    Parameters:
    config_file: str
        Path to the configuration file.

    Returns:
    dict
        A dictionary with configuration key-value pairs.
    """
    config = {}
    try:
        with open(config_file, 'r') as file:
            for line in file:
                # Ignore comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                # Split the line into key and value
                key, value = line.strip().split('=', 1)
                config[key.strip()] = value.strip()
    except FileNotFoundError:
        print(f"Error: Configuration file {config_file} not found.")
        exit(1)
    except Exception as e:
        print(f"Error while reading the configuration file: {e}")
        exit(1)
    
    return config

# Helper functions 

def clean_vector(vector, query, n=1, keep="before"):
    """
    Cleans elements of a pandas Series by removing a specified portion around a regex match.

    Parameters:
    vector: pd.Series
        The input Series to be cleaned.
    query: str
        The regex pattern to search for in each element of the Series.
    n: int
        The nth occurrence of the pattern to clean around (default is 1).
    keep: str
        Whether to keep the portion "before" or "after" the match (default is "before").

    Returns:
    pd.Series
        The cleaned Series.
    """
    def clean_element(element):
        matches = list(re.finditer(query, element))
        if len(matches) < n:
            return element
        match = matches[n - 1]
        if keep == "before":
            return element[:match.start()]
        else:
            return element[match.end():]
    return vector.apply(clean_element)

def remove_pos(vector, n):
    """
    Removes the last n characters from each string in a pandas Series.

    Parameters:
    vector: pd.Series
        The input Series.
    n: int
        Number of characters to remove from the end of each string.

    Returns:
    pd.Series
        The modified Series with truncated strings.
    """
    #Removes the last n characters from each string in a pandas Series.
    return vector.str.slice(stop=-n)

def get_orientations(positions, query):
    """
    Cleans each element of a pandas Series twice by removing specified portions around a regex match.

    Parameters:
    positions: pd.Series
        Series of genomic positions.
    query: str
        Regex pattern to search for.

    Returns:
    pd.Series
        Cleaned Series with orientation information.
    """
    #Cleans each element of a pandas Series twice by removing specified portions around a regex match.
    return clean_vector(clean_vector(positions, query, keep="after"), query, keep="after")

def browser_position(positions, prefix="chr"):
    """
    Adds a prefix to each element of a pandas Series if it does not already start with the prefix.

    Parameters:
    positions: pd.Series
        Series of genomic positions.
    prefix: str
        Prefix to add (default is "chr").

    Returns:
    pd.Series
        Series with positions formatted for browser visualization.
    """
    #Adds a prefix to each element of a pandas Series if it does not already start with the prefix.
    cleaned_positions = remove_pos(positions, 2)
    return np.where(cleaned_positions.str.startswith(prefix), cleaned_positions, prefix + cleaned_positions)

def get_family(IDs, query):
    """
    Extracts a portion from each element of a pandas Series based on a regex match.

    Parameters:
    IDs: pd.Series
        Series of element IDs.
    query: str
        Regex pattern to extract family information.

    Returns:
    pd.Series
        Series containing family names.
    """
    #Extracts a portion from each element of a pandas Series based on a regex match.
    return clean_vector(clean_vector(IDs, query, keep="after"), query, keep="before")

def get_regulation(values):
    """
    Determines the regulation status (UP or DN) for each value in a pandas Series.

    Parameters:
    values: pd.Series
        Series of numerical values representing log2FoldChange.

    Returns:
    pd.Series
        Series containing "UP" or "DN" labels based on the values.
    """
    #Determines the regulation status (UP or DN) for each value in a pandas Series.
    return np.where(values > 0, "UP", "DN")

def get_start_end(positions) -> pd.DataFrame:
    """
    Splits a pandas Series of genomic positions into separate columns for chromosome, start, and end positions.

    Parameters:
    positions: pd.Series
        Series containing genomic positions in the format "chr:start-end".

    Returns:
    pd.DataFrame
        DataFrame with separate columns for chromosome, start, and end positions.
    """
    #Splits a pandas Series of genomic positions into separate columns for chromosome, start, and end positions.
    chr_col = clean_vector(positions, ":", keep="before")
    cleaned_positions = clean_vector(positions, ":", keep="after")
    start_end = cleaned_positions.str.split('-', expand=True)
    start_end.columns = ['start', 'end']
    start_end.insert(0, 'chr', chr_col)
    return start_end

def makeNiceDESEQ(deseq2_table, annotation, element_pattern, directory, condition, control):
    """
    Processes a DESeq2 results table to extract and annotate TE loci, and prepares it for downstream analysis.

    Parameters:
    deseq2_table: str
        Path to the DESeq2 results file.
    annotation: str
        Path to the TE locus level annotation file.
    element_pattern: list
        List of patterns for filtering elements.
    directory: str
        Output directory for saving results.
    condition: str
        Name of the experimental condition.
    control: str
        Name of the control condition.

    Returns:
    tuple
        A DataFrame with annotated TE loci and a focus label.
    """
    #Load TE locus level annotation file 
    position_table = pd.read_table(annotation, sep='\t')
    
    data = pd.read_table(deseq2_table, sep='\t')
    
    data["ID"] = np.array(data.index)
    pattern = '|'.join(element_pattern)
    
    #Keep only specified elements
    data = data[data["ID"].str.contains(pattern, regex=True)]

    #Keep only Loci
    data = data[data["ID"].str.contains("dup")]

    regulation = data["log2FoldChange"]

    # Rank Metric
    baseMean_scaled = (data["baseMean"] - np.mean(data["baseMean"]))/np.std(data["baseMean"])
    baseMean_scaled = np.abs(np.min(baseMean_scaled)) + 1 + baseMean_scaled

    raw_p_scaled = (data["pvalue"] - np.mean(data["pvalue"]))/np.std(data["pvalue"])
    raw_p_scaled = np.abs(np.min(raw_p_scaled)) + 1 + raw_p_scaled

    data["RankMetric"] = np.log((1/raw_p_scaled)*baseMean_scaled + 1)*regulation

    og_IDs = data['ID'].astype(str).tolist()

    # Preprocess IDs
    IDs = clean_vector(pd.Series(og_IDs), ":")

    # Keep conserved loci
    niceDESEQ = position_table[position_table['TE'].isin(IDs)]

    # Reorder to keep original data order
    niceDESEQ = niceDESEQ.set_index('TE').loc[IDs].reset_index()

    # Rename columns
    niceDESEQ.columns = ['TE', 'position']

    # Add family
    family = clean_vector(niceDESEQ['TE'], "_")
    niceDESEQ['family'] = family

    # Retrieve locus strand
    niceDESEQ['strand'] = get_orientations(niceDESEQ['position'], ":")

    # Add regulation
    niceDESEQ['regulation'] = get_regulation(regulation)

    # Process positions to be browser ready
    rdy = browser_position(niceDESEQ['position'])
    niceDESEQ['position'] = rdy

    # Add start and end columns for easy check-up
    start_end = get_start_end(niceDESEQ['position'])
    niceDESEQ['chr'] = start_end['chr']
    niceDESEQ['start'] = start_end['start'].astype(int)
    niceDESEQ['end'] = start_end['end'].astype(int)

    # Add ranking metric
    # Copy data
    new_data = data.copy()

    # Get short version of row indices in original order
    data_TE_copy = clean_vector(new_data['ID'], ":")

    # Get indices of new_data where the id matches the id in niceDESEQ
    index = data_TE_copy.isin(niceDESEQ['TE'])

    # Reorder new_data based on the index
    new_data = new_data[index]

    # Extract RankMetric column values in the correct order
    score = new_data['RankMetric']

    niceDESEQ['score'] = score.values

    Focus = element_pattern[0] if len(element_pattern) == 1 else str(element_pattern[0]) + '_' + str(element_pattern[-1])
    if directory[-1] == '/':
        directory = directory[:-1]
    niceDESEQ.to_csv(f"{directory}/{Focus}_{condition}vs{control}.txt", sep='\t', index=False)
    return niceDESEQ, Focus

#Helper functions
def getBAM(filtered_bam_folder):
    """
    Retrieves forward and reverse BAM files for sample processing.

    Parameters:
    filtered_bam_folder: str
        Path to the folder containing filtered BAM files.

    Returns:
    list
        A list of dictionaries with keys: 'sample', 'fwd', 'rev'.
    """
    listdir_filtered = os.listdir(filtered_bam_folder)
    listdir_fwd = [file for file in listdir_filtered if 'fwd' in file and 'bai' not in file]
    listdir_rev = [file for file in listdir_filtered if 'rev' in file and 'bai' not in file]
    
    assert len(listdir_fwd) == len(listdir_rev), "Number of Forward and Reverse bam files do not match"
    dico_L = []
    for i in range(len(listdir_fwd)):
        fwd_bam = pysam.AlignmentFile(filtered_bam_folder + '/' + listdir_fwd[i], 'rb')
        rev_bam = pysam.AlignmentFile(filtered_bam_folder + '/' + listdir_rev[i], 'rb')
        dico_L.append({'sample':i,'fwd':fwd_bam, 'rev':rev_bam})
    return dico_L

def count_Coverage(chrom, start, end, bamfile):
    # Count reads in the region and print their flags
    read_count = 0
    for read in bamfile.fetch(chrom, start, end):
        read_count += 1
    return read_count

def getCoverage(region, fwd, rev, prefix=None):
    """
    Counts reads in a specified genomic region from a BAM file.

    Parameters:
    chrom: str
        Chromosome name.
    start: int
        Start position.
    end: int
        End position.
    bamfile: pysam.AlignmentFile
        BAM file to fetch reads from.

    Returns:
    int
        Number of reads in the specified region.
    """
    chr_n, start, end, strand = region

    if prefix:
        if prefix == 'remove_chr':
            chr_n = chr_n[3:]
        else:
            chr_n = prefix + chr_n

    if end - start == 0:
        end = start + 2

    bam_file = fwd if strand == '+' else rev
        
    try:
        coverage = count_Coverage(chr_n, start, end, bam_file)
    except:
        if 'KI' in chr_n or 'GL' in chr_n:
            chr_n = chr_n[:3] + 'Un_' + chr_n[3:-2] + "v" + chr_n[-1]
            try:
                coverage = count_Coverage(chr_n, start, end, bam_file)
            except:
                print(f"Chromosome not found ({chr_n})")
                return np.full(end - start, np.nan)
        else:
            print(f"Chromosome not found ({chr_n})")
            return np.full(end - start, np.nan)

    tot_cov = np.sum(coverage, axis=0)/(end-start)

    return tot_cov


def computeMetric(row, sample_List, window, decrease_indicator, prefix):
    """
    Computes 5' and 3' metrics for a given locus.

    Parameters:
    row: pd.Series
        A row from a DataFrame representing a genomic locus.
    sample_List: list
        List of sample dictionaries containing forward and reverse BAM files.
    window: int
        Window size used to compute 5' and 3' metrics.
    decrease_indicator: float
        Threshold for adjusting the 3' vs 5' ratio.
    prefix: str
        Prefix modification for chromosome names.

    Returns:
    tuple
        A tuple containing the computed metrics: (score5, score3, effect, coverage).
    """
    chr_n = row["chr"]
    start = row["start"]
    end = row["end"]
    strand = row["strand"]
    locus = (chr_n, start, end, strand)
            
    if strand == '+':
        # Get upstream region (compared to TE orientation)
        woi_5 = (chr_n,max(start - window, 0),start,strand)
        
        # Get downstream region (compared to TE orientation)
        woi_3 = (chr_n,end,end + window,strand)
    else:
        # Get upstream region (compared to TE orientation)
        woi_5 = (chr_n,end,end + window,strand)
        
        # Get downstream region (compared to TE orientation)
        woi_3 = (chr_n,max(start - window,0),start,strand)
    
    L =[]
    for i in range(len(sample_List)):
        # Get coverage in each sample
        sample = sample_List[i]
        
        L.append(getCoverage(woi_5, sample['fwd'], sample['rev'], prefix=prefix))
    
    score5 = np.mean(np.array(L) + 1)
        
    L =[]
    for i in range(len(sample_List)):
        # Get coverage in each sample
        sample = sample_List[i]
        L.append(getCoverage(woi_3, sample['fwd'], sample['rev'], prefix=prefix))

    score3 = np.mean(np.array(L) + 1)
    
    L =[]
    for i in range(len(sample_List)):
        # Get coverage in each sample
        sample = sample_List[i]
        L.append(getCoverage(locus, sample['fwd'], sample['rev'], prefix=prefix))
    
    coverage = np.mean(L)
    ratio = score3/score5
    coeff = 1/(1 - decrease_indicator)
    effect = np.power(ratio, 1 + coeff*ratio) if ratio <= 1 else np.tanh(ratio) + 1 - np.tanh(decrease_indicator)
    return (score5, score3, effect, coverage)
    

def getTupleInfo(row, INDEX = 0):
    """
    Extracts specific information from a tuple stored in a DataFrame row.

    Parameters:
    row: pd.Series
        A row containing a tuple.
    INDEX: int
        Index of the tuple element to retrieve.

    Returns:
    Any
        The value of the specified tuple element.
    """
    return row[0][INDEX]

def makeGenomicContextTable(loci_table, sample_List, window, prefix, decrease_3v5):
    """
    Generates a table containing genomic context scores for a list of loci.

    Parameters:
    loci_table: pd.DataFrame
        DataFrame containing loci information.
    sample_List: list
        List of sample dictionaries containing forward and reverse BAM files.
    window: int
        Window size used to compute metrics.
    prefix: str
        Prefix modification for chromosome names.
    decrease_3v5: float
        Threshold for adjusting the 3' vs 5' effect ratio.

    Returns:
    pd.DataFrame
        DataFrame containing computed genomic context scores.
    """
    MetricVector = pd.DataFrame(loci_table.apply(computeMetric, axis=1, 
                                      sample_List=sample_List,
                                      window=window,
                                      prefix=prefix,
                                      decrease_indicator=decrease_3v5))
    
    scores5 = MetricVector.apply(getTupleInfo, axis=1, INDEX=0)
    scores3 = MetricVector.apply(getTupleInfo, axis=1, INDEX=1)
    ratio = MetricVector.apply(getTupleInfo, axis=1, INDEX=2)
    coverage = MetricVector.apply(getTupleInfo, axis=1, INDEX=3)
    
    dico = {"5prime":scores5, "3prime":scores3, "3v5_effect":ratio, "MeanCoverage":coverage}
    return pd.DataFrame(data=dico)

# Get consensus File
def getConsensusSize(row):
    """
    Calculates the length of the consensus sequence for a given TE element.

    Parameters:
    row: pd.Series
        A row containing reverse and forward sequences of a TE.

    Returns:
    int
        The length of the consensus sequence.
    """
    return len(row["reverse"])
    
def read_fasta_to_df(file_path):
    """
    Reads a FASTA file and creates a DataFrame with element names, consensus sequences, and sequence lengths.

    Parameters:
    file_path: str
        Path to the FASTA file.

    Returns:
    tuple
        A DataFrame containing TE information and a dictionary mapping TE names to their lengths.
    """
    data = []
    
    for record in SeqIO.parse(file_path, "fasta"):
        element_name = record.id
        consensus_seq = str(record.seq)
        length = len(consensus_seq)
        
        data.append([element_name, consensus_seq, length])
    
    df = pd.DataFrame(data, columns=['TE', 'Consensus_Sequence', 'length'])
    dico = {family:df["length"][np.where(df["TE"] == family, True, False)].iloc[0] for family in df["TE"]}
    return df, dico

# make final Table
def IsFull(row, ref=None, coeff=0.9):
    """
    Determines if a locus is considered full-length based on its size relative to a reference consensus sequence.

    Parameters:
    row: pd.Series
        A row containing TE information.
    ref: dict
        A dictionary mapping TE families to their consensus sequence lengths.
    coeff: float
        The size coefficient threshold (default is 0.9).

    Returns:
    int
        1 if the locus is full-length, 0 otherwise.
    """
    try:
        if row["end"] - row["start"] < coeff*ref[row["family"]]: 
            return 0
        else:
            return 1
    except:
        print(f"\nNo consensus found for {row['family']}")
        return np.nan
    
def getSizeRatio(row, ref=None):
    """
    Calculates the size ratio of a locus relative to its consensus sequence.

    Parameters:
    row: pd.Series
        A row containing TE information.
    ref: dict
        A dictionary mapping TE families to their consensus sequence lengths.

    Returns:
    float
        The size ratio of the locus.
    """
    try:
        return (row["end"] - row["start"])/ref[row["family"]]
    except:
        print(f"\nNo consensus found for {row['family']}")
        return np.nan

def vec_LeakyReLU(array, slope):
    """
    Applies a vectorized Leaky ReLU activation function to an array.

    Parameters:
    array: np.array
        Input array.
    slope: float
        Slope of the negative side of the activation function.

    Returns:
    np.array
        Array after applying the Leaky ReLU function.
    """
    return np.where(array > 0, array, slope*array)
    
def vec_transformCoverage(array, start_decrease):
    """
    Transforms coverage values based on a threshold for decreased effect.

    Parameters:
    array: np.array
        Array of coverage values.
    start_decrease: float
        Threshold below which the effect decreases.

    Returns:
    np.array
        Array with transformed coverage values.
    """
    # Calculate the constant 'c'
    c = (1/start_decrease)*np.log(1/start_decrease) if start_decrease != 0 else np.nan
    
    return np.where(array > start_decrease, np.tanh(array) + (1 - np.tanh(start_decrease)), array*np.exp(c*array))

def vec_transformSize(array, start_decrease):
    """
    Transforms size ratios using a logarithmic decay function based on a threshold.

    Parameters:
    array: np.array
        Array containing size ratios to be transformed.
    start_decrease: float
        The threshold below which the transformation applies a decay function.

    Returns:
    np.array:
        Transformed size ratios.
    """
    # Calculate the constant 'c'
    c = (1/start_decrease)*np.log(1/start_decrease) if start_decrease != 0 else np.nan
    
    return np.where(array > start_decrease, 1, array*np.exp(c*array))

def PredictActivation(Table, path_to_model):
    """
    Predicts activation status of loci using a pre-trained model.

    Parameters:
    Table: pd.DataFrame
        DataFrame containing the feature columns needed for the prediction.
    path_to_model: str
        Path to the pickle file containing the pre-trained model.

    Returns:
    np.array:
        Predicted activation status for each locus.
    """
    
    features = ["score", "3v5_effect", 'MeanCoverage', 'size_ratio', "Metric"]
    df = Table[features]
    
    with open(path_to_model, 'rb') as f:
        model = pkl.load(f)
    activation = model.predict(df[features])
    
    return activation

def makeTable(deseq, sample_dico, window, prefix, decrease_3v5, consensus_dico, size_threshold, coverage_threshold, full_length_threshold, path_to_model):
    """
    Constructs the final table by computing genomic context scores, size metrics, and predicting activation.

    Parameters:
    deseq: pd.DataFrame
        DataFrame containing the DESeq2 output with TE loci information.
    sample_dico: list of dict
        List containing dictionaries with BAM file information for each sample.
    window: int
        Window size for 5' and 3' coverage calculations.
    prefix: str
        Chromosome prefix or modification for BAM file chromosome names.
    decrease_3v5: float
        Threshold for the decrease in 3' vs 5' metric.
    consensus_dico: dict
        Dictionary mapping TE family names to their consensus sequence length.
    size_threshold: float
        Threshold for size ratio transformation.
    coverage_threshold: float
        Threshold for coverage transformation.
    full_length_threshold: float
        Threshold for defining a locus as full-length.
    path_to_model: str
        Path to the pre-trained activation prediction model.

    Returns:
    tuple:
        (Table, summary_dict) where `Table` is the final DataFrame and `summary_dict` contains summary statistics.
    """
    Table = deseq
        
    t0 = time.time()
    print("\nStarted computing genomic context scores..")
    debut = time.time()
    context = makeGenomicContextTable(Table, sample_dico, window=window, prefix=prefix, decrease_3v5=decrease_3v5)
    
    Table["3v5_effect"] = context["3v5_effect"]
    
    Table["MeanCoverage"] = context["MeanCoverage"] 
    Table["Coverage_effect"] = vec_transformCoverage(Table["MeanCoverage"], coverage_threshold)
    print(f"Finished. Computation time is {time.time()-debut} seconds.\n")
    
    
    
    print("Processing length informations of loci..")
    debut = time.time()
    Table["size_ratio"] = Table.apply(getSizeRatio, axis=1, ref=consensus_dico)
    Table["full_length"] = Table.apply(IsFull, axis=1, ref=consensus_dico, coeff=full_length_threshold)
    Table["size_effect"] = vec_transformSize(Table["size_ratio"], size_threshold)
    print(f"Finished. Computation time is {time.time()-debut} seconds.\n")
    
    print("Computing Metric and identifying activated loci..")
    debut = time.time()
    Table["Metric"] = Table["score"]*Table["3v5_effect"]*Table["Coverage_effect"]*Table["size_effect"]
    
    noNA_table = Table[~pd.isna(Table["Metric"])]
    NA_table = Table[pd.isna(Table["Metric"])]
    
    active = PredictActivation(noNA_table, path_to_model)
    
    noNA_table["Activated"] = active
    
    NA_table["Activated"] = [np.nan for i in range(Table[pd.isna(Table["Metric"])].shape[0])]
    Table = pd.concat([noNA_table, NA_table], ignore_index=True)
    
    UP_mask = (Table["regulation"] == "UP") & (~pd.isna(Table["Metric"]))
    
    print(f"Finished. Computation time is {time.time()-debut} seconds.\n")
    print(f"Table ready. Total computation time was {time.time()-t0} seconds. Happy analysis :D\n")
    DN_mask = (Table["regulation"] == "DN") & (~pd.isna(Table["Metric"]))
    summary_dict = {"window":window,
                    "decrease_3v5":decrease_3v5,
                    "size_threshold":size_threshold,
                    "coverage_threshold":coverage_threshold,
                    "full_length_threshold":full_length_threshold,
                    "dataset_size":Table.shape[0],
                    "mean_metric_UP":np.mean(Table[UP_mask]["Metric"]),
                    "std_metric_UP":np.std(Table[UP_mask]["Metric"]),
                    "mean_metric_DN":np.mean(Table[DN_mask]["Metric"]),
                    "std_metric_DN":np.std(Table[DN_mask]["Metric"]),
                    "max_coverage_full_length_UP":np.max(Table[UP_mask]["MeanCoverage"]),
                    "max_metric":np.max(Table[UP_mask]["Metric"]),
                    "UP_proportion":Table[UP_mask].shape[0]/Table.shape[0],
                    "DN_proportion":Table[DN_mask].shape[0]/Table.shape[0]}
    return Table, summary_dict

# Main function to run your pipeline
def run_pipeline(config):
    """
    Main function to run the locus-level analysis pipeline using the provided configuration.

    Parameters:
    config: dict
        Configuration settings loaded from the configuration file.

    Returns:
    None
    """
    print("Configuration Loaded:")
    for key, value in config.items():
        print(f"{key}: {value}")
    
    # Access configuration values using config.get('<KEY>')
    tellam_dir = config.get('TELLAM_DIR')
    deseq2_table = config.get('DESEQ2_TABLE')
    bam_state = config.get('BAM_STATE')
    chr_prefix = config.get('CHR_PREFIX')
    if chr_prefix == 'None':
        chr_prefix = None
    condition_name = config.get('CONDITION_NAME')
    control_name = config.get('CONTROL_NAME')
    element_pattern = config.get('ELEMENT_PATTERN').split(',')
    filtered_bam_folder = config.get('FILTERED_FOLDER')
    directory = config.get('DIRECTORY')
    annotation = config.get('ANNOTATION')
    consensus = config.get('CONSENSUS')
    window = int(config.get('WINDOW'))
    context = float(config.get('CONTEXT'))
    size = float(config.get("SIZE"))
    coverage = float(config.get('COVERAGE'))
    full = float(config.get('FULL'))
    
    print(f"\nRunning pipeline for {condition_name} vs {control_name}...")
    
    # Load filtered BAM files 
    sample_List = getBAM(filtered_bam_folder)
    print("\nLoaded filtered BAM files.")
    
    # Load consensus file
    table, consensus_dico = read_fasta_to_df(consensus)
    print("\nLoaded consensus file.")
    
    print("\nProcessing DESeQ2 table...")
    # Process deseq table
    niceDESEQ, Focus = makeNiceDESEQ(deseq2_table, 
                              annotation, 
                              element_pattern, 
                              directory, 
                              condition_name, 
                              control_name)
    print("\nProcessed DESeQ2 table.")
    
    print("\nComputing TELLAM...")
    # Making table
    Table, summary_dict = makeTable(niceDESEQ, sample_List, window, chr_prefix,
                            context,
                            consensus_dico, 
                            size, # High decrease in size effect below 10% of consensus size
                            coverage, # High decrease in coverage effect below 0.06 reads per base
                            full, 
                            f"{tellam_dir}/ActivationModel.pkl") # proportion of consensus length above which loci are considered full_length
    print("\nTELLAM computed...")
    # Formating and Saving tables of interest
    Table.rename(columns={'chr': '#chr'}, inplace=True)

    Activated = Table[Table["Activated"] == 1]

    bed_format = ["#chr", "start", "end", "TE", "family", "position", "strand", "score", "3v5_effect", "MeanCoverage", "size_ratio", "full_length", "size_effect", "Metric", "Activated"]

    Table[bed_format].to_csv(f"{directory}/{Focus}_{condition_name}vs{control_name}.bed", sep='\t', index=False)
    Activated[bed_format].to_csv(f"{directory}/active_{Focus}_{condition_name}vs{control_name}.bed", sep='\t', index=False)

    with open(f"{directory}/{Focus}_{condition_name}vs{control_name}_summary.pkl", 'wb') as f:
        pkl.dump(summary_dict, f)
    
    print("\nTable saved. TELLAM Pipeline completed.")

if __name__ == "__main__":
    # Load the configuration file path from command-line arguments
    import sys
    if len(sys.argv) != 2:
        print("Usage: python3 your_pipeline_script.py <config_file>")
        exit(1)
    
    config_file = sys.argv[1]
    config = load_config(config_file)
    run_pipeline(config)
