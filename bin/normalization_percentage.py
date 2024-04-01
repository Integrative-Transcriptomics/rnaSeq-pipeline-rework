import argparse

def normalization_percentage(feature_counts_file):

    counts_data = []
    read_counts_sum = []
    # Read File input
    with open(feature_counts_file, 'r') as file:
        for line in file:
            if (not line.startswith('#')) and (not line.startswith("Geneid")):
                line_list = line.strip().split('\t')
                counts_data.append(line_list)
                
    # Add Pseudocount to data
    number_of_genes = len(counts_data)
    number_of_samples = len(counts_data[0])
    for i in range(0, number_of_genes):
        for j in range(7, number_of_samples):
            counts_data[i][j] = float(counts_data[i][j]) + 1

    samples = []
    for gene in counts_data:
        samples.append(gene[7:])
    
    number_of_samples = len(samples[0])
    
    # 1) Sum up all readcounts
    for i in range(0, number_of_samples):
        sum = 0
        for gene in samples:
            sum += float(gene[i])
            
        read_counts_sum.append(sum)
        
    results = []
    # 2) Divide all read counts by the sum (percentage)
    for gene in samples:
        percentage_values = []
        value = 0
        for i in range(0, number_of_samples):
            value = float(gene[i]) / read_counts_sum[i]
            percentage_values.append(value)
            
        results.append(percentage_values)

    return results
