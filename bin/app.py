from dash import Dash, html, dash_table, dcc, callback, Output, Input
import pandas as pd
import plotly.express as px

# Initialize the app
app = Dash(__name__)


df = pd.read_csv('/Users/sarina/Bachelorarbeit/rnaSeq-pipeline-rework/rRNADepletion/rRNA_remaining.csv')

# App layout
app.layout = html.Div([
    html.Div(className='row', children='My First App with Data',
             style={'textAlign': 'center', 'color': 'black', 'fontSize': 30, 'fontFamily': 'Arial'}),
    html.Hr(),
    html.Div(className ='row', style =
             {'color': 'black', 'fontSize': 15, 'fontFamily': 'Arial'},
             children = 
             [dcc.RadioItems(options=['ratio', 'rRNA_count'], value='ratio', inline=True, id='controls-and-radio-item')]),
    html.Div(className='row', 
             style = {'color': 'black', 'fontSize': 10, 'fontFamily': 'Arial'},
             children=[dash_table.DataTable(data=df.to_dict('records'), page_size=10)]),
    dcc.Graph(figure={}, id='controls-and-graph'),
])

# Add controld to build interaction
@callback(
    Output(component_id='controls-and-graph', component_property='figure'),
    Input(component_id='controls-and-radio-item', component_property='value')
)

def update_graph(col_chosen):
    fig = px.histogram(df, x=col_chosen, y='file')
    return fig

#Run the app
if __name__ == '__main__':
    app.run(debug=True)
    