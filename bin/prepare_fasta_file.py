#!/usr/bin/env python3

import argparse
import pickle

# Extracts sequence from fasta file
def parse_fasta_file(fasta_file, genome_reference):
    
    sequence = {} # {'genome_reference' : 'sequence'}
    genomic_reference = ""
    found_genome_reference = False
    
    # Open fasta file and write into data structure
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">") and genome_reference in line:
                line_list = line[1:].strip().split(" ")
                genomic_reference = line_list[0]
                sequence[genomic_reference] = ''
                found_genome_reference = True
            elif line.startswith(">") and genome_reference not in line:
                found_genome_reference = False
            elif found_genome_reference:
                current_sequence = sequence[genomic_reference]
                sequence_part = line.strip()
                sequence[genomic_reference] = current_sequence + sequence_part
        
    with open('fasta_data', 'wb') as file:
        pickle.dump(sequence, file)


def main():
    
    # Create the argument parser
    parser = argparse.ArgumentParser(description="An example script with command-line arguments.")

    # Add command-line arguments
    parser.add_argument("-f", "--fasta_file", help="Path to the fasta file", required=True)
    parser.add_argument("-gf", "--genome_reference", help="Genome reference")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    
    fasta_file = args.fasta_file
    genome_reference = args.genome_reference
    
    parse_fasta_file(fasta_file, genome_reference)

main()