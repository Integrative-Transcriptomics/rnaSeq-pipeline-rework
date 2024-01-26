import plotly.express as px
import pandas as pd
import argparse
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import tpm_calculation as tpm 
import normalization_percentage as norm_percent
import variance_calculation as vc
import gff_parser as gff_parser


# creates a bar chart, with all samples and the given rRNA percentage
def bar_chart_procent(rRNA_remaining):
    
    file_data = []
    samples = []
    ratios =[]
    
    with open(rRNA_remaining) as file:
        for line in file:
            file_data.append(line.strip().split(','))
    
    for sample in file_data[1:]:
        samples.append(sample[0])
        ratios.append(float(sample[-1]))
        
    data_dict = {
        'Samples' : samples,
        'Ratio' : ratios,
    }
    
    data = pd.DataFrame(data_dict)
        
    fig = px.bar(data, x='Ratio', y='Samples')
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        
    fig.write_html('ratio_plot.html', auto_open=True)
    

# Creates a bar chart for each sample, indicating how much 5S, 16S and 23S rRNA is still present.
# tpm normalization and normalization in percent
def bar_chart_different_rRNA_tpm(rRNA_genes_type, counts_txt):
    
    genes_data = [] #[[gene_id, type], [gene_id, type], ...]
    rRNA_genes = [] #[gene_id, gene_id, ...]
    counts_data = [] # [[Geneid, Chr, Start, End, Strand, Length, gene_name, sample1, sample2, ...], [..., ..., ..., ...]];field 7 -> end 
    rRNA_data_tpm = [] # [[Geneid, Chr, Start, End, Strand, Length, gene_name, sample1, sample2, ...], [..., ..., ..., ...]];field 7 -> end only rRNA genes
    rRNA_data_percent = [] # [[Geneid, Chr, Start, End, Strand, Length, gene_name, sample1, sample2, ...], [..., ..., ..., ...]];field 7 -> end only rRNA genes
    genes = []
    type = []
    number_of_arrays = 0
    data_dict_tpm = {}
    data_dict_percent = {}
    tpm_values = tpm.tpm_calculation(counts_txt)
    norm_percent_values = norm_percent.normalization_percentage(counts_txt)
    
    # Extract rRNA genes and their tyoe from file
    with open(rRNA_genes_type) as genes_file:
        for line in genes_file:
            line_list = line.strip().split(",")
            genes_data.append(line_list)
            rRNA_genes.append(line_list[0])

    # Read File input
    with open(counts_txt, 'r') as file:
        for line in file:
            if line.startswith("Geneid"):
                top_line = line.strip().split('\t')
                rRNA_data_tpm.append(top_line)
                rRNA_data_percent.append(top_line)
            if (not line.startswith('#')) and (not line.startswith("Geneid")):
                line_list = line.strip().split('\t')
                counts_data.append(line_list)
    
    # Replace read counts with calculated tpm values.
    number_of_genes = len(counts_data)
    for i in range(0, number_of_genes):
        entery = None
        if counts_data[i][0] in rRNA_genes:
            entery_tpm = counts_data[i][:7] + tpm_values[i]
            entery_percent = counts_data[i][:7] + norm_percent_values[i]
            rRNA_data_tpm.append(entery_tpm)
            rRNA_data_percent.append(entery_percent)
    
    # Add genes and their type to the corresponding lists 
    for gene in genes_data[1:]:
        if len(gene) == 2:
            genes.append(gene[0])
            type.append(gene[1])
        else:
            raise NameError("Input file that contains gene id and type is not correct.")
      
    # Add genes to dictionary      
    data_dict_tpm["genes"] = genes 
    data_dict_percent["genes"] = genes
    
    # create correct number of arrays
    number_of_arrays = len(rRNA_data_tpm[0]) - 7
    print("Feature counts has information about " + str(number_of_arrays) + " samples.")
    
    for i in range(0, number_of_arrays):
        data_dict_tpm[f"{rRNA_data_tpm[0][i + 7]}"] = []
        data_dict_percent[f"{rRNA_data_tpm[0][i + 7]}"] = []
    
    for count in rRNA_data_tpm[1:]:
        for i in range(0, number_of_arrays):
            #number_tpm = float(count[i + 7])
            number_tpm = np.log10(float(count[i + 7]))
            data_dict_tpm[f"{rRNA_data_tpm[0][i + 7]}"].append(number_tpm)
    
    for count in rRNA_data_percent[1:]:
        for i in range(0, number_of_arrays):
            #number_percent = np.log10(float(count[i+7]))
            number_percent = float(count[i+7])
            data_dict_percent[f"{rRNA_data_tpm[0][i + 7]}"].append(number_percent)
    
    # Add type to data dict
    data_dict_tpm["type"] = type  
    data_dict_percent["type"] = type
    
    # Prepare for bar chart
    data_tpm = pd.DataFrame(data_dict_tpm)
    data_percent = pd.DataFrame(data_dict_percent)
    
    #create dict with unique color values
    color_dict={}
    x=0
    for t in type:
        if t not in color_dict:
            color_dict[t]=x
            x += 1
            
    #create color list
    color_list=[]
    for t in type:
        color_list.append(color_dict[t])
    
    # Create Figure
    fig = make_subplots(rows=1, cols=2)
    
    for column in data_tpm.columns[1:-1]:
        fig.add_trace(go.Bar(
            x=data_tpm["genes"],
            y= data_tpm[column],
            name=column,
            marker_color=color_list,
            legendgrouptitle_text = "Normalization with TPM:",
            legendgroup = 'a',
        ),
        row=1,
        col=1)
    
    for column in data_percent.columns[1:-1]:
        fig.add_trace(go.Bar(
            x=data_percent["genes"],
            y= data_percent[column],
            name=column,
            marker_color=color_list,
            legendgrouptitle_text = "Normalization in percent:",
            legendgroup = 'b',
        ),
        row=1,
        col=2)
        
        
    fig.update_layout(
        barmode='group',
        title='Overview over genes and their type',
        xaxis1_title="Genes",
        yaxis1_title="Expression number (logarithmic scaling)",
        xaxis2_title="Genes",
        yaxis2_title="Expression number",
        autotypenumbers='convert types',
        legend=dict(groupclick="toggleitem"),
        updatemenus=[
            dict(type="buttons",
                buttons=[
                    dict(label="Select all",
                        method="restyle",
                        args=[{"visible": True}]),
                    dict(label="Select none",
                         method="restyle",
                         args=[{"visible": "legendonly"}])
                ])],    
    )
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        
    fig.write_html('genes_barchart.html', auto_open=True)


# Histogramm for readcounts
def readcounts_histogram(counts_txt):
    
    counts_data = []
    counts_data_transformed = []
    dict_data = {}
    
    with open(counts_txt) as file:
        for line in file:
            if not(line.startswith("#")):
                counts_data.append(line.strip().split("\t")[7:])
    
    number_of_samples = len(counts_data[0])
    number_of_genes = len(counts_data)
    
    for i in range(0, number_of_samples):
        intermediate_step = []
        
        for j in range(1, number_of_genes):
            counts = float(counts_data[j][i])
            if counts == 0: # Pseudocount
                intermediate_step.append(counts + 1)
            else:
                intermediate_step.append(counts)

        counts_data_transformed.append(intermediate_step)
        
    fig = go.Figure()
    
    number_of_sample = 0
    
    for sample in counts_data_transformed:
        for count in sample:
            if str(count) in dict_data:
                dict_data[str(count)] += 1
            else:
                dict_data[str(count)] = 1
            
        x_values = []
        y_values = []
    
        for key, value in dict_data.items():
            x_values.append(np.log10(float(key)))
            y_values.append(value)
        
        #fig.add_trace(go.Histogram(x=x_values, name=counts_data[0][number_of_sample], marker=dict(color='rgba(0,0,0,0)', line=dict(color='black', width=2))))
        fig.add_trace(go.Histogram(x=x_values, name=counts_data[0][number_of_sample]))
        number_of_sample += 1
    
    fig.update_layout(
        barmode='overlay',
        title='Histogram of the distribution of read counts per sample (output of featureCounts)',
        xaxis_title="Read frequency per Gene (logarithmic scaling)",
        yaxis_title="Frequency",
        bargap=0,
        showlegend=True,
        updatemenus=[
            dict(type="buttons",
                buttons=[
                    dict(label="Select all",
                        method="restyle",
                        args=[{"visible": True}]),
                    dict(label="Select none",
                         method="restyle",
                         args=[{"visible": "legendonly"}])
        ]
    )])
    fig.update_traces(opacity=0.25)
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        
    fig.write_html('read_counts_histogram.html', auto_open=False)


def create_table(rRNA_genes_type, counts_txt, gff_file):
    #["mean_values" , "empirical_variance", "standard_deviation", "coefficient_of_variation", "gene"]
    header_values = ['Mean', 'Empirical variance', 'Standard deviation', 'Coefficient of variation', 'Gene ID', 'Product']
    genes = []
    gff_content = gff_parser.gff_file_parser(gff_file)
    
    top_10_results = vc.calculate_empirical_variance(rRNA_genes_type, counts_txt)
    
    # Extract all genes from the top results
    for result in top_10_results:
        genes.append(result[-1])
        
    for line in gff_content:
        product = None
        if "locus_tag" in line[-1]:
            gene_id = line[-1]["locus_tag"]
            if gene_id in genes:
                if "product" in line[-1]:
                    product = line[-1]["product"]
                    for gene in top_10_results:
                        if gene[-1] == gene_id:
                            gene.append(product)
        
    # Swap columns with rows
    cells_values = [[row[i] for row in top_10_results] for i in range(len(top_10_results[0]))]
        
    
    fig = go.Figure(data=[go.Table(header=dict(values=header_values, align='left'),
                 cells=dict(values=cells_values, align='left'))
                     ])
        
    fig.update_layout(
        title = "Top 10 genes with highest expression and low variance"
    )
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
    
    fig.write_html('genes_barchart.html', auto_open=True)

# Calls plot functions 
def main():
    
    # Create the argument parser
    parser = argparse.ArgumentParser(description="An example script with command-line arguments.")

    # Add command-line arguments
    parser.add_argument("-if", "--rRNA_remaining", help="Path to the remaining rRNA file", required=True)
    parser.add_argument("-genesf", "--rRNA_genes", help="Path to the file containing all rRNA genes", required=True)
    parser.add_argument("-fcf", "--feature_counts_file", help="Output file counts.txt of featureCounts", required=True)
    parser.add_argument("-gff", "--gff_file", help="Gff file input.", required=True)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    rRNA_remaining = args.rRNA_remaining
    rRNA_genes = args.rRNA_genes
    feature_counts = args.feature_counts_file
    gff_file = args.gff_file
    
    open("all_plots.html", "w").close()
    
    bar_chart_procent(rRNA_remaining)
    bar_chart_different_rRNA_tpm(rRNA_genes, feature_counts)
    readcounts_histogram(feature_counts)
    create_table(rRNA_genes, feature_counts, gff_file)
    
main()