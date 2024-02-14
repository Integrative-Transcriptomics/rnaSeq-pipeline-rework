import os
import pandas as pd

# File opens counts file and saves data to data structure
def prepare_counts_file(counts_file):
    
    # each entry contains one gene
    counts_data = [] 
    #[[Geneid, Chr, Start, End, Strand, Length, gene_name, saamples ...], [SCO0001, ...], ...]
    
    # extract all lines from file, except top line
    with open(counts_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                line_list = line.strip().split('\t')
                counts_data.append(line_list)
    
    return counts_data

# File opens specific bedgraph file and returns dictionary of the genomic references and their
# positions and counts how many reads map to each position
def prepare_bedgraph_file(bedgraph_file):

    bedgraph_data = {} # {'NC_03888.3' : {1 : 0; 2 : 0, 3 : 0, ...}, 'NC_03888.2' : {1 : 0; 2 : 0; ....}}
    genome_reference = ""
    
    with open(bedgraph_file, "r") as file:
        for line in file:
            line_list = line.strip().split('\t')
            if genome_reference != line_list[0]:
                genome_reference = line_list[0]
                bedgraph_data[genome_reference] = {}
            else:
                key = int(line_list[1])
                value = int(line_list[2])
                bedgraph_data[genome_reference][key] = value
    
    return bedgraph_data

# Extracts start and end position of specific gene
def get_start_and_end_position(counts_data, gene, genome_reference):
    
    # iterate over each entry of counts data
    for gene_entry in counts_data:
        # if gene id of entry in counts data corresponds to searched gene and genome reference
        # return start and end position 
        if gene_entry[0] == gene and gene_entry[1] == genome_reference:
            return([gene_entry[2], gene_entry[3]]) # [start, end]
     
    
    
# Takes start and end postion of the gene and calculates all postion inbetween
# returns data set, that can be used for a line chart
def prepare_dataset_for_graph(start_and_end_position, bedgraph_data, genome_reference):
    
    # Data 
    x = []
    y = []
    start = 0
    end = 0
    
    # If there are multiple genome references, choose current one
    current_bedgraph_data = bedgraph_data[genome_reference]
    
    # Start and end position
    if len(start_and_end_position) == 2:
        start = int(start_and_end_position[0])
        end = int(start_and_end_position[1])
    else:
        raise ValueError('Gene could not be found and has no start and end position.')
    
    
    # Iterate over all positions between start and end
    for i in range(start, end):
        y.append(current_bedgraph_data[i])
        x.append(i)
    
    # Prepate dataset for plotly graph
    data = pd.DataFrame({'x' : x, 'y' : y})

    return data 

