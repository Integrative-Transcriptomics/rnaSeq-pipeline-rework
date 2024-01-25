import argparse

"""
    command for execution:
    python3 normalization_percentage.py -fcf /Users/sarina/Documents/Bachelorarbeit/rnaSeq-pipeline-rework/rRNADepletion/counts.txt
"""

def normalization_percentage(feature_counts_file):

    counts_data = []
    read_counts_sum = []
    # Read File input
    with open(feature_counts_file, 'r') as file:
        for line in file:
            if (not line.startswith('#')) and (not line.startswith("Geneid")):
                line_list = line.strip().split('\t')
                counts_data.append(line_list)

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
    
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="An example script with command-line arguments.")

    # Add command-line arguments
    parser.add_argument("-fcf", "--feature_counts_file", help="Output file counts.txt of featureCounts", required=True)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    feature_counts_file = args.feature_counts_file
    
    normalization_percentage(feature_counts_file)

#main()