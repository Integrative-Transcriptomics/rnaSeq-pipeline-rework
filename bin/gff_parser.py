import argparse

"""
    command for execution:
    python3 gff_parser.py -gff /Users/sarina/Bachelorarbeit/rnaSeq-pipeline-rework/GCF_000203835.1_ASM20383v1_genomic.gff
"""

# This function parses a gff file, and saves it's content into a list
def gff_file_parser(gff_file):
    
    gff_content = [] #[[seqname, source, feature, start, end, score, strand, frame, attribute], [], ...]
    attribute_dict = {}
    
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                line_list = line.strip().split('\t')
                gff_content.append(line_list)
                
    for element in gff_content:
        attribute = element[-1]
        attribute_list = attribute.split(';')
        for pair in attribute_list:
            key, value = pair.split('=')
            attribute_dict[key] = value
       # print(attribute_dict)
        element[-1] = attribute_dict # e.g. element  {"Name": "GeneID"}
        attribute_dict = {}

    #print(gff_content[:5])
    return gff_content
        
    
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="An example script with command-line arguments.")

    # Add command-line arguments
    parser.add_argument("-gff", "--gff_file", help="Path to the gff file.", required=True)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    gff_file = args.gff_file
    
    gff_content = gff_file_parser(gff_file)
    
#main()