# Calculates TPM Value for each gene
"""
1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.
"""
def tpm_calculation(feature_counts_file):
    
    counts_data = []
    
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
    
    rpk_values = []
    per_million_values = []
    tpm_values = []
    
    # 1) Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    for gene in counts_data:
        
        length = float(gene[5])
        length_in_kbp = length/1000
        
        samples = gene[7:]
        rpk_sample = [float(element)/length_in_kbp for element in samples]
        
        rpk_values.append(rpk_sample)

    number_of_samples = len(rpk_values[0])
    per_million_values = []
    
    # 2) Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
    for i in range(0, number_of_samples):
        sum = 0
        for gene in rpk_values:
            sum = sum + gene[i]
        
        per_million_factor = sum/1000000
        per_million_values.append(per_million_factor)
            
    # 3) Divide the RPK values by the “per million” scaling factor. This gives you TPM.
    for i in range(0, number_of_samples):
        tpm = 0
        tpm_genes = []
        
        for gene in rpk_values:
            tpm = gene[i] / per_million_values[i]
            tpm_genes.append(tpm)
            
        tpm_values.append(tpm_genes)
        
    tpm_values_transposed = [[row[i] for row in tpm_values] for i in range(len(tpm_values[0]))] #[[sample1, sample2, ....], [sample1, sample2, ...], ...]
  
    return tpm_values_transposed
