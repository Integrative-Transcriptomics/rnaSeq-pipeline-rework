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

# Takes counts file and extracts gene names from file
def extract_gene_names(counts_file):
    
    # Call function that extracts data from counts file
    counts_data = prepare_counts_file(counts_file)
    gene_names = []
    
    # add each gene to gene_names list
    for gene in counts_data[1:]:
        gene_names.append(gene[0])
    
    return gene_names

# takes counts file and extracts sample names from file
def extract_sample_names(counts_file):
    
    # call function that extracts data from file
    counts_data = prepare_counts_file(counts_file)
    sample_names = []
    
    # add each sample name to sample_names list
    for sample in counts_data[0][7:]:
        sample_names.append(sample)
        
    return sample_names
    
    
    