import plotly.express as px
import pandas as pd
import argparse
import plotly.graph_objects as go
from numpy import log as ln
import numpy as np

# creates a bar chart, with all samples and the given rRNA percentage
def bar_chart_procent(rRNA_remaining):
    
    file_data = []
    samples = []
    ratios =[]
    
    with open(rRNA_remaining) as file:
        for line in file:
            file_data.append(line.strip().split(','))
    
    for sample in file_data[2:]:
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
def bar_chart_different_rRNA(rRNA_genes_type, counts_txt):
    
    genes_data = [] #[[gene_id, type], [gene_id, type], ...]
    counts_data = [] # [[Geneid, Chr, Start, End, Strand, Length, gene_name, sample1, sample2, ...], [..., ..., ..., ...]];field 7 -> end 
    top_line_added = False
    genes = []
    type = []
    number_of_arrays = 0
    data_dict = {
        "genes" : [],
        #"type" : []
    }
    
    with open(rRNA_genes_type) as genes_file:
        for line in genes_file:
            genes_data.append(line.strip().split(","))
    
    with open(counts_txt) as counts_file:
        for line in counts_file:
            for gene in genes_data[1:]:
                if line.startswith("Geneid") and not(top_line_added):
                    counts_data.append(line.strip().split('\t'))
                    top_line_added = True
                if line.startswith(gene[0]):
                    counts_data.append(line.strip().split('\t'))    # weiß nicht ob das so effizient ist??
    
    for gene in genes_data[1:]:
        if len(gene) == 2:
            genes.append(gene[0])
            type.append(gene[1])
        else:
            raise NameError("Input file that contains gene id and type is not correct.")
        
        
    data_dict["genes"] = genes
    
    
    # create correct number of arrays
    number_of_arrays = len(counts_data[0]) - 7
    print("Feature counts has information about" + str(number_of_arrays) + " samples.")
    
    for i in range(0, number_of_arrays):
        data_dict[f"{counts_data[0][i + 7]}"] = []
    
    for count in counts_data[1:]:
        for i in range(0, number_of_arrays):
            number = ln(float(count[i + 7]))
            data_dict[f"{counts_data[0][i + 7]}"].append(number)
    
    
    data_dict["type"] = type  
    
    # Prepare for bar chart
    data = pd.DataFrame(data_dict)
    
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
    
    fig = go.Figure()

    x = 0
    for column in data.columns[1:-1]:
        fig.add_trace(go.Bar(
            x=data["genes"],
            y= data[column],
            name=column,
            marker_color=color_list,
        ))
        

        
    fig.update_layout(
        barmode='group',
        title='Overview over genes and their type',
        xaxis_title="Genes and Type",
        yaxis_title="Expression number (logarithmic scaling)",
        autotypenumbers='convert types',
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
            x_values.append(ln(float(key)))
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
        
    fig.write_html('read_counts_histogram.html', auto_open=True)




# Calls plot functions 
def main():
    
    # Create the argument parser
    parser = argparse.ArgumentParser(description="An example script with command-line arguments.")

    # Add command-line arguments
    parser.add_argument("-if", "--rRNA_remaining", help="Path to the remaining rRNA file", required=True)
    parser.add_argument("-genesf", "--rRNA_genes", help="Path to the file containing all rRNA genes", required=True)
    parser.add_argument("-fcf", "--feature_counts_file", help="Output file counts.txt of featureCounts", required=True)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    rRNA_remaining = args.rRNA_remaining
    rRNA_genes = args.rRNA_genes
    feature_counts = args.feature_counts_file
    
    open("all_plots.html", "w").close()
    
    bar_chart_procent(rRNA_remaining)
    bar_chart_different_rRNA(rRNA_genes, feature_counts)
    readcounts_histogram(feature_counts)
     
    
main()