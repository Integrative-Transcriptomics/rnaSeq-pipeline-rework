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
import random


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
        
    fig = px.bar(data, x='Ratio', y='Samples', color_discrete_sequence=['darkgrey']*len(data))
    
    fig.update_layout(
        title='Percentage of the number of rRNA reads in the individual samples',
        xaxis_title="Number of reads in percent",
        yaxis_title="Individual samples",
        yaxis=dict(tickfont=dict(size=12),),
        height= 600,
        width = 1400
    )
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        
    fig.write_html('ratio_plot.html', auto_open=False)
    

# Creates a bar chart for each sample, indicating how much 5S, 16S and 23S rRNA is still present.
# tpm normalization and normalization in percent
def bar_chart_different_rRNA(rRNA_genes_type, counts_txt):
    
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
            number_percent = float(count[i+7])
            data_dict_percent[f"{rRNA_data_tpm[0][i + 7]}"].append(number_percent)
            
    # Add type to data dict
    data_dict_tpm["type"] = type  
    data_dict_percent["type"] = type
    
    # Prepare for bar chart
    data_tpm = pd.DataFrame(data_dict_tpm)
    data_percent = pd.DataFrame(data_dict_percent)
    
    data_sum_rRNA_genes = calculate_bar_plot_total_percent_rRNA(data_dict_percent, rRNA_data_percent[0][7:])
    
    #create dict with unique color values
    color_dict={
        "5S": "darkblue",
        "16S": "yellow",
        "23S": "hotpink"
    }

    #create color list
    color_list=[]
    for t in type:
        color_list.append(color_dict[t])
    
    # Create Figure
    fig = make_subplots(rows=2, cols=2,
                        specs=[[{}, {}], [{'colspan': 2}, None]],
                        subplot_titles=("Normalization with TPM", "Normalization in percent", "Sum of all percent values, per rRNA type"))
    
    for column in data_tpm.columns[1:-1]:
        fig.add_trace(go.Bar(
            x =data_tpm["genes"],
            y = data_tpm[column],
            name = column,
            marker_color = color_list,
            legendgrouptitle_text = "Normalization with TPM:",
            legendgroup = 'a',
        ),
        row=1,
        col=1)
    
    for column in data_percent.columns[1:-1]:
        fig.add_trace(go.Bar(
            x = data_percent["genes"],
            y = data_percent[column],
            name = column,
            marker_color = color_list,
            legendgrouptitle_text = "Normalization in percent:",
            legendgroup = 'b',
        ),
        row=1,
        col=2)
    
    color_sum_plot = ['darkblue', 'yellow', 'hotpink']
    
    for column in data_sum_rRNA_genes.columns[:-1]:
        fig.add_trace(go.Bar(
            x = data_sum_rRNA_genes["type"],
            y = data_sum_rRNA_genes[column],
            name = column,
            marker_color = color_sum_plot,
            legendgrouptitle_text = "Sum of all percent values",
            legendgroup = 'c',
            ),
            row=2, col=1) 
        
    # Add color legend
    for type, color_type in color_dict.items():
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(color=color_dict[type], opacity=1),
            name=type,
            legendgrouptitle_text = "Color Overview",
            hoverinfo='none',
            legendgroup = 'd',
            legendrank=1,
        ),
)
        
        
        
    fig.update_layout(
        barmode='group',
        #title='5S, 16S and 23S Expression in the different samples',
        xaxis1_title="Genes",
        yaxis1_title="Expression number (logarithmic scaling)",
        xaxis2_title="Genes",
        yaxis2_title="Expression number",
        autotypenumbers='convert types',
        legend=dict(groupclick='toggleitem'),
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
        width=1400,  
    )
    
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        
    fig.write_html('genes_barchart.html', auto_open=False)
    
    return fig
    
    
# Bar chart showing the percentage of reads that are 5S, 16S or 23S rRNA
def calculate_bar_plot_total_percent_rRNA(data_dict_percent, sample_names):

    type_list = data_dict_percent['type']    
    percent_data = [] # [[sample_name, [(type, percent), (type, percent), ...]], [sample_name, [(type, percent), (type, percent), ...]], ]
    
    for sample_name in sample_names:
        sample_to_add = [sample_name]
        sample_to_add.append(list(zip(type_list, data_dict_percent[sample_name])))
        percent_data.append(sample_to_add)
    
    percent_sum_data_dict = {} # [[sample_name, sum 5S genes, sum 16S genes, sum 23S genes], [sample_name, sum 5S genes, sum 16S genes, sum 23S genes]]
    
    for sample in percent_data:
        fiveS_sum = 0
        sixteen_sum = 0
        twentythree_sum = 0
        sample_name = sample[0]
        sample_data = []
        for data in sample[1]:
            type = data[0]
            percent_value = data[1]
            if type == "5S":
                fiveS_sum += percent_value
            if type == "16S":
                sixteen_sum += percent_value
            if type == "23S":
                twentythree_sum += percent_value
                
        sample_data.append(fiveS_sum)
        sample_data.append(sixteen_sum)
        sample_data.append(twentythree_sum)

        percent_sum_data_dict[sample_name] = sample_data

    percent_sum_data_dict['type'] = ["5S", "16S", "23S"]
    
    return pd.DataFrame(percent_sum_data_dict)    
    

# Histogramm for readcounts
def readcounts_histogram(counts_txt, rRNA_genes_type):
    
    counts_data = []
    counts_data_without_rRNA = []
    rRNA_genes = []
    counts_data_transformed = []
    counts_data_transformed_without_rRNA = []
    dict_data1 = {}
    dict_data2 = {}
    
    # Add all rRNA genes from file to rRNA_genes list
    with open(rRNA_genes_type) as file:
        for line in file:
            if not line.startswith('gene_id'):
                line_list = line.strip().split(",")
                rRNA_genes.append(line_list[0])
    
    # Create data with and without rRNA
    with open(counts_txt) as file:
        for line in file:
            if not(line.startswith("#")):
                line_list = line.strip().split("\t")
                new_list = [line_list[0]] + line_list[7:]
                if line.startswith("Geneid"):
                    counts_data.append(new_list)
                    counts_data_without_rRNA.append(new_list)
                elif line_list[0] in rRNA_genes:
                    counts_data_without_rRNA.append(new_list)
                else:
                    counts_data.append(new_list)
    
    number_of_samples = len(counts_data[0][1:])
    number_of_genes = len(counts_data)
    
    # Add Pseudocount for each entry
    for i in range(1, number_of_samples + 1):
        intermediate_step = []
        intermediate_step_without_rRNA = []
        for j in range(1, number_of_genes): 
            counts = float(counts_data[j][i])
            intermediate_step.append(counts+1)
            
        for j in range(1, len(counts_data_without_rRNA)):
            counts = float(counts_data_without_rRNA[j][i])
            intermediate_step_without_rRNA.append(counts+1)
            
        counts_data_transformed.append(intermediate_step)
        counts_data_transformed_without_rRNA.append(intermediate_step_without_rRNA)
            
    fig = make_subplots(rows = 2,
                        cols = 1,
                        shared_xaxes=True,
                        shared_yaxes=True,
                        vertical_spacing = 0,
                        subplot_titles=("Upper plot without rRNA genes, lower plot only with rRNA genes", None)
                        )
    
    all_colors = ["#{:06x}".format(i) for i in range(0xFFFFFF + 1)]  # Liste aller möglichen Farben
    random_colors = random.sample(all_colors, number_of_samples) 
    
    max_value_for_y_axis = 0
    # Iterate over each sample
    for i in range(0, number_of_samples):
        # count how often read count appears for all genes
        for count in counts_data_transformed[i]:
            if str(count) in dict_data1:
                dict_data1[str(count)] += 1
            else:
                dict_data1[str(count)] = 1
                
        for count in counts_data_transformed_without_rRNA[i]:
            # count how ofen read count appears only for rRNA genes
            if str(count) in dict_data2:
                dict_data2[str(count)] += 1
            else:
                dict_data2[str(count)] = 1

        x_values_rRNA = []
        y_values_rRNA = []
    
        x_values_normal = []
        y_values_normal = []
    
        for key, value in dict_data1.items():
            x_values_normal.append(np.log10(float(key)))
            y_values_normal.append(value)
    
        for key, value in dict_data2.items():
            x_values_rRNA.append(np.log10(float(key)))
            y_values_rRNA.append(value)
        
        curent_max_value = max(y_values_normal)
        if curent_max_value > max_value_for_y_axis:
            max_value_for_y_axis = curent_max_value

        fig.add_trace(go.Histogram(x=x_values_normal, name=counts_data[0][i+1], 
                                   legendgroup ='a', 
                                   legendgrouptitle_text = "Samples without rRNA genes",
                                   marker_color=random_colors[i],
                                   xbins=dict(size=0.09)
                                   ),
                      row=1, col=1)
        fig.add_trace(go.Histogram(x=x_values_rRNA, y=y_values_rRNA, name=counts_data_without_rRNA[0][i+1],
                                   legendgroup='b',
                                   legendgrouptitle_text = "Samples with only rRNA genes",
                                   marker_color=random_colors[i],
                                   xbins=dict(size=0.09)
                                   ),
                      row=2, col=1)
    
    #fig.update_yaxes1(range=[0,90])
    fig.update_layout(
        barmode='overlay',
        #title='Histogram of the distribution of read counts per sample',
        yaxis1_title="Frequency",
        yaxis2_title="Frequency",
        #yaxis2=dict(autorange='reversed', range=[0,90]),
        yaxis1_range=[0,350],
        yaxis2_range=[0,350],
        #yaxis2_autorange='reversed',
        bargap=0,
        showlegend=True,
        legend=dict(groupclick='toggleitem'),
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
        width=1400)
    fig.update_traces(opacity=0.25)
    fig.update_xaxes(title_text="Read frequency per Gene (logarithmic scaling)", row=2, col=1)

    
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        
    fig.write_html('read_counts_histogram.html', auto_open=False)
    
    return fig


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
    
    # Create csv file for dash table with tooltip
    with open("table_data.csv", "w") as file:
        # Add header values to file
        file.write(",".join(header_values) + "\n")
        number_of_columns = len(header_values)
        print(number_of_columns)
        number_of_rows = 10
        
        for i in range(0, number_of_rows):
            new_line = ""
            for j in range(0, number_of_columns):
                new_line = new_line + str(cells_values[j][i]) + ","
            # Write each line to file    
            file.write(new_line[:-1] + "\n")
    
    fig = go.Figure(data=[go.Table(header=dict(values=header_values, align='left'),
                 cells=dict(values=cells_values, align='left'))
                     ])
        
    fig.update_layout(
        title = "Top 10 genes with highest expression and low variance",
        width=1400
    )
    
    with open("all_plots.html", "a") as plot_file:
        plot_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
    
    fig.write_html('table_highest_expression_low_variance.html', auto_open=False)

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
    bar_chart_different_rRNA(rRNA_genes, feature_counts)
    readcounts_histogram(feature_counts, rRNA_genes)
    create_table(rRNA_genes, feature_counts, gff_file)
    
#main()