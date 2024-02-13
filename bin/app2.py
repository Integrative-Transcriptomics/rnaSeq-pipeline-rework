from dash import Dash, html, dash_table, dcc, callback, Output, Input
import pandas as pd
import plotly.express as px
import extract_names
import sys
import dash_bootstrap_components as dbc
import os
import prepare_data_line_plot as pdlp
import pickle

with open('fasta_data', "rb") as file:
    sequence_data = pickle.load(file)

# Initialize the app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = 'rRNA depletion quality control'
# Prepare dropdown

gene_names = extract_names.extract_gene_names('/Users/sarina/Bachelorarbeit/rnaSeq-pipeline-rework/rRNADepletion/counts.txt')
sample_names = extract_names.extract_sample_names('/Users/sarina/Bachelorarbeit/rnaSeq-pipeline-rework/rRNADepletion/counts.txt')

# App layout
app.layout = html.Div([
    html.Header(className='header', children=['Quality control for rRNA depletion in prokaryotes'],
                style=dict(fontSize = '25px')),
    html.Br(),
    html.Div(className='row',
             children = [
                 html.Div(className = 'dropdown',
                          children = [
                              dcc.Dropdown(
                                  id = 'dropdown_samples',
                                  options = sample_names,
                                  value = str(sample_names[0]),
                                  placeholder = 'Select sample'
                              )
                          ],
                          style = dict(width='50%')),
                 html.Div(className = 'dropdown',
                          children = [
                              dcc.Dropdown(
                                  id = 'dropdown_genes',
                                  options = gene_names,
                                  value = str(gene_names[0]),
                                  placeholder = 'Select gene'
                              )
                          ],
                          style = dict(width='50%')),
             ],
             style=dict(display='flex')),
    html.Div(className = 'Line Plot',
             children = [
                 dcc.Graph(id = 'line_plot')
             ]),
    html.Div(className = 'text', 
             children=['Select an area in the figure above to display the corresponding sequence.']),
    html.Div(className = 'Sequence output', id='sequence_output',style={'white-space': 'pre-wrap', 'word-wrap': 'break-word'}),
    ], style={'margin' : '40px'})

@callback(
    Output("line_plot", "figure"),
    [Input('dropdown_samples', "value"),
     Input('dropdown_genes', 'value')]
)

def update_line_plot(sample, gene):
    
    current_dict = os.getcwd()
    parrent_dict = os.path.dirname(current_dict)
    sample_name, extension = os.path.splitext(sample)
    
    if '.sorted' in sample_name:
        sample_name = sample_name.replace('.sorted', '')

    path_counts_txt = parrent_dict + '/rRNADepletion/counts.txt'
    path_bedgraph = parrent_dict + '/rRNADepletion/' + sample_name + '.bedgraph'
    
    genome_reference = 'NC_003888.3'
    
    counts_data = pdlp.prepare_counts_file(path_counts_txt)
    bedgraph_data = pdlp.prepare_bedgraph_file(path_bedgraph)
    
    start_and_end_position = pdlp.get_start_and_end_position(counts_data, gene, genome_reference)
    data = pdlp.prepare_dataset_for_graph(start_and_end_position, bedgraph_data, genome_reference)
    
    fig = px.line(data, x = 'x', y = 'y', color_discrete_sequence=['darkgrey'])
    
    fig.update_layout(
        title = f'Overview over position coverage for the {gene} gene in {sample}',
        xaxis_title='Position',
        yaxis_title='Coverage per position',
    )
    
    return fig

@callback(
    Output("sequence_output", "children"),
    [Input("line_plot", 'relayoutData'),
    #Input("fasta_file_input", "value")
    ]
)

def get_sequence(relayout_data):
    
    start = '0'
    end = '0'
    genome_reference = 'NC_003888.3'
    sequence = ''
    
    if relayout_data:
        if 'xaxis.range[0]' in relayout_data and 'xaxis.range[1]' in relayout_data:
            start = round(float(relayout_data["xaxis.range[0]"]))
            end = round(float(relayout_data["xaxis.range[1]"]))
        else:
            sequence = ""

    if start != '0' and end != '0':
        # Berechnung für Position 
        sequence = sequence_data[genome_reference][start:end]
        
    return f'Sequence from position {start} to {end}:\t {sequence}'
        

#Run the app
if __name__ == '__main__':
    app.run(debug=True)
