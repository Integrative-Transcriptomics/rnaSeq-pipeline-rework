from dash import Dash, html, dash_table, dcc, callback, Output, Input
import pandas as pd
import plotly.express as px
import extract_names
import sys
import dash_bootstrap_components as dbc
import os
import prepare_data_line_plot as pdlp
import pickle
import plot
import time

try:
    with open('fasta_data', "rb") as file:
        sequence_data = pickle.load(file)

    # Initialize the app
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.title = 'rRNA depletion quality control'

    # Table Data
    current_directory = os.getcwd()
    parrent_directory = os.path.dirname(current_directory)
    data_rRNA_ratio = pd.read_csv(parrent_directory + '/Results/rRNADepletion/rRNA_remaining.csv')
    counts_file = parrent_directory + '/Results/rRNADepletion/counts.txt'
    rRNA_type_file= parrent_directory + '/rRNA_genes_type_ecoli.txt'
    data_expression = pd.read_csv(parrent_directory + '/Results/rRNADepletion/table_data.csv')

    # Prepare dropdown
    gene_names = extract_names.extract_gene_names(counts_file)
    sample_names = extract_names.extract_sample_names(counts_file)

    # Genome reference names
    genome_reference = list(sequence_data.keys())
    if not genome_reference:
        raise ValueError("No Reference was defined")
        
    # App layout
    app.layout = html.Div([
        html.Header(className='header', children=['Quality control for rRNA depletion in prokaryotes'],
                    style=dict(fontSize = '25px')),
        html.Br(),
        html.Div(className = 'RadioItems', style = {'color': 'black', 'fontSize': 15, 'fontFamily': 'Arial'},
                children  = [
                    html.H5(children = ['Percentage of the number of rRNA read in the individual samples']),
                    dcc.RadioItems(options = ['total_count', 'rRNA_count', 'ratio'], value='ratio', id='controls-and-radio-item'),]),
        html.Div(className = 'table',
                style=dict(font = 'arial'),
                children = [
                    dash_table.DataTable(data=data_rRNA_ratio.to_dict('records'),
                                        columns = [{"name" :i, "id": i } for i in data_rRNA_ratio.columns],
                                        style_cell={'textAlign': 'left'},
                                        page_size=10,
                                        sort_action='native',)
                ]),
        html.Div(className = 'barplot',
                children = [
                    dcc.Graph(figure = {}, id = 'controls-and-graph')]),
        html.Div(className='row',
                children = [
                    html.H5(children = ['Wiggle plot for all annotated genes']),
                    html.H5(children = ['Given Reference Genome:\t' + str(genome_reference[0])]),
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
                            style = dict(width='50%'))
                ],
                style=dict(display='flex')),
        html.Br(),
        html.Div(className = 'Line Plot',
                children = [
                    dcc.Graph(id = 'line_plot'),
                ]),
        html.Div(className = 'text', 
                children=['Select an area in the figure above to display the corresponding sequence.']),
        html.Div(className = 'Sequence output', id='sequence_output',style={'white-space': 'pre-wrap', 'word-wrap': 'break-word'}),
        html.Br(),
        html.Div(className = 'table', id='table',
                children = [html.H5(children = ['Overview over the top 10 genes with a high expression rate and low variance, normalized with TPM']),
                            dash_table.DataTable(
                                data = data_expression.to_dict('records'),
                                columns = [{"name" :i, "id": i } for i in data_expression.columns],
                                tooltip_header={
                                    'Mean': 'The mean value represents the average expression rate of the genes. Calculated from the different samples. ',
                                    #'Empirical variance': 'The empirical variance indicates how much the individual expression rates of the samples fluctuate around the mean value.',
                                    'Standard deviation': 'The standard deviation is a number that represents the dispersion or variation of the different expression rates around the mean value. The value is attributed to the same scale as the initial values.',
                                    'Coefficient of variation': 'The coefficient of variation is a relative measure of the dispersion or variation from the data to the mean value of the expression rates.',
                                    'Gene ID' : 'Is a unique designation that is assigned to a gene.',
                                    'Product' : 'The product attribute can contain various information, such as the name of the protein encoded by the gene or a description of its function or property.'
                                    },
                                tooltip_delay = 0,
                                tooltip_duration = None,
                                sort_action='native',
                                style_cell={'textAlign': 'left', 'font-family': 'arial'},
                                css=[{
                                    'selector': '.dash-table-tooltip',
                                    'rule': 'background-color: grey; font-family: monospace; color: white'
                                    }],
                                )],
                style=dict(font='arial')
                            ),
        html.Br(),
        html.Div(className =  'bar_chart_different_rRNA_types',
                children = [
                    html.H5(children = ['5S, 16S and 23S expression in the different samples']),
                    dcc.Graph(figure = plot.bar_chart_different_rRNA(rRNA_type_file, counts_file), style = {'height' : '800px'}),
                ]),
        html.Br(),
        html.Div(className = 'readcounts_histogram',
                children = [
                    html.H5(children = ['Histogram of the distribution of read counts per sample']),
                    dcc.Graph(figure = plot.readcounts_histogram(counts_file, rRNA_type_file), style = {'height' : '800px'})
                ])
        ], style={'margin' : '40px', 'height' : 'auto'})

    # Bar Plot with rRNA ratio
    @callback(
        Output(component_id='controls-and-graph', component_property='figure'),
        Input(component_id='controls-and-radio-item', component_property='value')
    )

    def update_graph(col_chosen):
        
        fig = px.bar(data_rRNA_ratio, x=col_chosen, y='file', color_discrete_sequence =['darkgrey']*len(data_rRNA_ratio))
        
        fig.update_layout(
            xaxis_title= col_chosen,
            yaxis_title="Individual samples",
            yaxis=dict(tickfont=dict(size=12),)
        )
        return fig


    # Line plot with 
    @callback(
        Output("line_plot", "figure"),

        [Input('dropdown_samples', "value"),
        Input('dropdown_genes', 'value'),]
    )

    def update_line_plot(sample, gene):
        time.sleep(5)
        current_dict = os.getcwd()
        try:
            parrent_dict = os.path.dirname(current_dict)
            sample_name, extension = os.path.splitext(sample)
            
            if '.sorted' in sample_name:
                sample_name = sample_name.replace('.sorted', '')

            path_counts_txt = parrent_dict + '/Results/rRNADepletion/counts.txt'
            path_bedgraph = parrent_dict + '/Results/rRNADepletion/' + sample_name + '.bedgraph'
            
            counts_data = pdlp.prepare_counts_file(path_counts_txt)
            bedgraph_data = pdlp.prepare_bedgraph_file(path_bedgraph, genome_reference[0])
            
            start_and_end_position = pdlp.get_start_and_end_position(counts_data, gene, genome_reference[0])
            data = pdlp.prepare_dataset_for_graph(start_and_end_position, bedgraph_data, genome_reference[0])
            
            fig = px.line(data, x = 'x', y = 'y', color_discrete_sequence=['darkgrey'])
            
            fig.update_layout(
                title = f'Overview over position coverage for the {gene} gene in {sample} ',
                xaxis_title='Position',
                yaxis_title='Coverage per position',
            )
            return fig
        except:
            data = pd.DataFrame({'x': [], 'y': []})
            
            fig = px.line(data, x = 'x', y = 'y', color_discrete_sequence=['darkgrey'])
            
            fig.update_layout(
                title = 'Please choose a sample and a gene',
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
        
        
except:
    raise NameError('Please start Pipeline before executing the app')