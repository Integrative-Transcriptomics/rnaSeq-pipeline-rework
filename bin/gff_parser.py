
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

        element[-1] = attribute_dict # e.g. element  {"Name": "GeneID"}
        attribute_dict = {}

    return gff_content
