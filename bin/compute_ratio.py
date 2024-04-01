#!/usr/bin/env python3

import pandas as pd
import sys, getopt
#from glob import glob



def main(argv):
    countsfile = ''
    rRNA_genes = ''
    
    try:
        opts, args = getopt.getopt(argv,"hc:r:",["counts=","rRNAgenes="])
    except getopt.GetoptError:
        print('test.py -c <counts> -r <rRNAgenes>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-c", "--counts"):
            countsfile = arg
        elif opt in ("-r", "--rRNAgenes"):
            rRNA_genes_path = arg
    # read input file
    with open(rRNA_genes_path) as file:
        contents = file.read()
    # split lines
    rRNA_genes_and_type = contents.split('\n')
    # change data structure of rRNA_genes
    rRNA_genes = []
    # split gene from type and save only the gene 
    for gene_and_txpe in rRNA_genes_and_type[1:]:
        rRNA_genes.append(gene_and_txpe.split(",")[0])

    count_cols = "file total_count rRNA_count ratio".split()
    counts = []
    counts_df = pd.read_table(countsfile, skiprows=1)
    for col_name in counts_df.columns[7:]:
        rRNA_df = counts_df.loc[counts_df["Geneid"].isin(rRNA_genes)]
        total_count = counts_df[col_name].sum()
        rRNA_count = rRNA_df[col_name].sum()
        counts.append([col_name, total_count, rRNA_count, rRNA_count/total_count])
    pd.DataFrame(counts, columns = count_cols).to_csv("rRNA_remaining.csv", index=False)
    print("Done.")

if __name__ == "__main__":
   main(sys.argv[1:])
