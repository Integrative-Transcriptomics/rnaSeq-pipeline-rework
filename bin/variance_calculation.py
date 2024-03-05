#!/usr/bin/env python3

import plotly.express as px
import pandas as pd
import argparse
import plotly.graph_objects as go
from numpy import log as ln
import numpy as np
import math
import sys, getopt
from glob import glob
import gff_parser
import tpm_calculation as tpm

# calculates the empirical variance
def calculate_empirical_variance(rRNA_genes_type, counts_txt):
    
    # Datasets 
    genes_data = [] #[[gene_id, type], [gene_id, type], ...]
    counts_data = [] # [[Geneid, Chr, Start, End, Strand, Length, gene_name, sample1, sample2, ...], [..., ..., ..., ...]];field 7 -> end 
    not_rRNA_genes = []
    
    # Values in TPM format
    tpm_values = tpm.tpm_calculation(counts_txt)

    with open(rRNA_genes_type) as genes_file:
        for line in genes_file:
            gene, type = line.strip().split(",")
            genes_data.append(gene)
            
    counter = -2
            
    with open(counts_txt) as counts_file:
        for line in counts_file:
            current_line = line.strip().split('\t')
            if current_line[0] not in genes_data[1:]:
                if not current_line[0].startswith('#'):
                    if counter < 0 :
                        counts_data.append(current_line)
                    else:
                        counts_data.append([current_line[0]] + tpm_values[counter])
            counter += 1
            
    number_of_genes = len(counts_data)
    number_of_samples = len(counts_data[0][1:])
    
    mean_values = []    
    variance_calculated_genes = []
    
    
    # 1) Calculate mean value
    for i in range(1, number_of_genes):
        sum = 0
        for value in counts_data[i][1:]:
            sum += float(value)
            
        mean = round(sum/number_of_samples, 2)
        
        # Limit for minimum expression rate: 10
        if mean > 10:
            mean_values.append(mean)
            variance_calculated_genes.append(counts_data[i][0])
            
    empirical_variance = []
    standard_deviation = []
    count = 0
    
    # 2) Calculate variance and standard deviation
    for i in range(1, number_of_genes):
        sum = 0
        variance = 0
        
        if counts_data[i][0] in variance_calculated_genes:
            for value in counts_data[i][1:]:
                sum += (float(value) - mean_values[count]) ** 2
                
            variance = round(sum/(number_of_samples - 1), 2)
            empirical_variance.append(variance)
            
            standard_deviation.append(round(math.sqrt(variance), 2))
            count += 1
            
    coefficient_of_variation = []
    # 3) Calculate coefficient of variation
    for std, mean in zip(standard_deviation, mean_values):
        coefficient_of_variation.append(round(std/mean, 2))
        
    results = [[a,b,c,d,e] for a,b,c,d,e in zip(mean_values, empirical_variance, standard_deviation, coefficient_of_variation, variance_calculated_genes)]
        
    # Sort results in two different ways to find subset
    # 1) Expression rate from highest do lowest (mean value, with TPM calculation)
    sorted_mean = sorted(results, key=lambda x: (x[0]), reverse=True)
    
    first_100_genes = sorted_mean[:100]
    
    #2) Variance from lowest to highest (empirical variance)
    sorted_variance_first_100 = sorted(first_100_genes, key=lambda x: (x[1]))
    
    return(sorted_variance_first_100[:10])
            
def create_table(rRNA_genes_type, counts_txt, gff_file):
    #["mean_values" , "empirical_variance", "standard_deviation", "coefficient_of_variation", "gene"]
    header_values = ['Mean', 'Standard deviation', 'Coefficient of variation', 'Gene ID', 'Product']
    genes = []
    gff_content = gff_parser.gff_file_parser(gff_file)
    
    top_10_results = calculate_empirical_variance(rRNA_genes_type, counts_txt)
    
    # Extract all genes from the top results
    for result in top_10_results:
        genes.append(result[-1])
        
    for line in gff_content:
        product = None
        if "locus_tag" in line[-1]:
            gene_id = line[-1]["locus_tag"]
            if gene_id in genes:
                if "product" in line[-1]:
                    product = line[-1]["product"]
                    for gene in top_10_results:
                        if gene[-1] == gene_id:
                            gene.append(product)
        
    # Swap columns with rows
    cells_values = [[row[i] for row in top_10_results] for i in range(len(top_10_results[0]))]
    
    # Create csv file for dash table with tooltip
    with open("table_data.csv", "w") as file:
        # Add header values to file
        file.write(",".join(header_values) + "\n")
        number_of_columns = len(header_values)
        number_of_rows = 10
        
        for i in range(0, number_of_rows):
            new_line = ""
            for j in range(0, number_of_columns + 1 ):
                if j != 1:#
                    new_line = new_line + str(cells_values[j][i]) + ","
            # Write each line to file    
            file.write(new_line[:-1] + "\n")

    

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="An example script with command-line arguments.")

    # Add command-line arguments
    parser.add_argument("-genesf", "--rRNA_genes", help="Path to the file containing all rRNA genes", required=True)
    parser.add_argument("-fcf", "--feature_counts_file", help="Output file counts.txt of featureCounts", required=True)
    parser.add_argument("-gff", "--gff_file", help="Gff file", required=True)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    rRNA_genes = args.rRNA_genes
    feature_counts = args.feature_counts_file
    gff = args.gff_file
    
    create_table(rRNA_genes, feature_counts, gff)
    
main()