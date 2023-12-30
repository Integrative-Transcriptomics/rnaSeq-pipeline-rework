import plotly.express as px
import pandas as pd
import argparse
import plotly.graph_objects as go
from numpy import log as ln
import numpy as np
import sys, getopt
from glob import glob

# calculates the empirical variance
def calculate_empirical_variance(rRNA_genes_type, counts_txt):
    
    genes_data = [] #[[gene_id, type], [gene_id, type], ...]
    counts_data = [] # [[Geneid, Chr, Start, End, Strand, Length, gene_name, sample1, sample2, ...], [..., ..., ..., ...]];field 7 -> end 
    not_rRNA_genes = []
    
    with open(rRNA_genes_type) as genes_file:
        for line in genes_file:
            gene, type = line.strip().split(",")
            genes_data.append(gene)
    
    # Add all non rRNA genes to counts_data     
    with open(counts_txt) as counts_file:
        for line in counts_file:
            current_line = line.strip().split('\t')
            if current_line[0] not in genes_data[1:]:
                if not current_line[0].startswith("#"):
                    not_rRNA_genes.append(current_line[0])
                    counts_data.append(current_line)

    not_rRNA_genes = not_rRNA_genes[1:]
    mean_values = []

    # 1) Caclulate mean value
    for gene in counts_data[1:]:
        sum = 0
        number_of_sanples = len(gene[7:])
        for value in gene[7:]:
           sum += float(value) 
        mean_values.append(sum/number_of_sanples)
    
    empirical_variance = []
    count = 0
    # calculate variance
    for gene in counts_data[1:]:
        sum = 0
        number_of_sanples = len(gene[7:])
        for value in gene[7:]:
            sum += (float(value) - mean_values[0]) ** 2
            
        empirical_variance.append(sum/(number_of_sanples-1))    
        count += 1
    
    
    difference_between_varaiance_mean = 0   
    
    results = [[x, y, z] for x,y,z in zip(mean_values, empirical_variance, not_rRNA_genes)]
    
    for result in results:
        difference_between_varaiance_mean = abs(result[1] - result[0])
        result.append(difference_between_varaiance_mean)
    
    sorted_results = sorted(results, key=lambda x: x[-1])
    print(sorted_results[-10:])


def main():
        # Create the argument parser
    parser = argparse.ArgumentParser(description="An example script with command-line arguments.")

    # Add command-line arguments
    parser.add_argument("-genesf", "--rRNA_genes", help="Path to the file containing all rRNA genes", required=True)
    parser.add_argument("-fcf", "--feature_counts_file", help="Output file counts.txt of featureCounts", required=True)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    rRNA_genes = args.rRNA_genes
    feature_counts = args.feature_counts_file
    
    calculate_empirical_variance(rRNA_genes, feature_counts)
    
main()

""" 
command for execution:

"""