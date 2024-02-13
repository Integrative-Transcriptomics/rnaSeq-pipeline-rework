import argparse
import pickle

# Extracts sequence from fasta file
def parse_fasta_file(fasta_file):
    
    sequence = {} # {'genome_reference' : 'sequence'}
    genomic_reference = ""
    
    # Open fasta file and write into data structure
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                line_list = line[1:].strip().split(" ")
                genomic_reference = line_list[0]
                sequence[genomic_reference] = ''
            else:
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
    
    # Parse the command-line arguments
    args = parser.parse_args()
    
    fasta_file = args.fasta_file
    
    parse_fasta_file(fasta_file)

main()