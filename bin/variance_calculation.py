import plotly.express as px
import pandas as pd
import argparse
import plotly.graph_objects as go
from numpy import log as ln
import numpy as np
import math
import sys, getopt
from glob import glob

""" 
command for execution:
python3 variance_calculation.py -genesf /Users/sarina/Documents/Bachelorarbeit/rnaSeq-pipeline-rework/rRNA_genes_type.txt 
-fcf /Users/sarina/Documents/Bachelorarbeit/rnaSeq-pipeline-rework/rRNADepletion/counts.txt
"""

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

    number_of_genes = len(counts_data)
    number_of_samples = len(counts_data[0][7:])
    # print(number_of_genes)
    # print(number_of_samples)
    
    variance_calculated_genes = []
    
    # 1) Caclulate mean value
    for i in range(1, number_of_genes):
        sum = 0
        for value in counts_data[i][7:]:
            sum += float(value)
            
        mean = sum/number_of_samples
        
        # Limit for minimum expression rate: 10
        if mean > 10:
            mean_values.append(mean)
            variance_calculated_genes.append(counts_data[i][0])     # When the mean expression rate is above 10 save it to the mean values
                
    empirical_variance = []
    standard_deviation = []
    count = 0
    
    # 2) Calculate variance and standard deviation
    for i in range(1, number_of_genes):
        sum = 0
        variance = 0
        if counts_data[i][0] in variance_calculated_genes:
            for value in counts_data[i][7:]:
                sum += (float(value) - mean_values[count]) ** 2
            
            variance = sum/(number_of_samples-1)
            empirical_variance.append(variance)
            
            standard_deviation.append(math.sqrt(variance))
            count += 1 
    
    coefficient_of_variation = []        
    # 3) Calculate coefficient of variation
    for std, mean in zip(standard_deviation, mean_values):
        coefficient_of_variation.append(std/mean)
    
    # print(len(variance_calculated_genes))
   
    results = [[a,b,c,d,e] for a,b,c,d,e in zip(mean_values, empirical_variance, standard_deviation, coefficient_of_variation, variance_calculated_genes)]
    
    # Sort resulst in two different ways and find subset
    # 1) Expression rate from highest to lowest (mean value)
    sorted_results_mean_expression = sorted(results, key=lambda x: (x[0]), reverse=True)
    first_100_genes = sorted_results_mean_expression[:100]
    
    # 2) Variance from lowest to highest (empirical variance)
    # sorted_results_variance = sorted(results, key=lambda x: (x[1]))
    sorted_first_100_genes_after_variance = sorted(first_100_genes, key=lambda x: (x[1]))
    
    # Limit to find intersection: 2200 (in this case only two elements are the same) 
    #intersection = [element for element in sorted_results_mean_expression[:2200] if element in sorted_results_variance[:2200]]
    
    # Ouput of the genes with the smallest variance but high expression
    print(sorted_first_100_genes_after_variance[:10])
    return(sorted_first_100_genes_after_variance[:10])
    
    # print(intersection)
    
    

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
