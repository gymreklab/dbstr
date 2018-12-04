# -*- coding: utf-8 -*-

import dash
from dash.dependencies import Output, Input, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
import plotly
import plotly.graph_objs as go
from collections import deque
import pandas as pd
import sqlite3
import numpy as np


def run_query_withparms(sql):
    conn = sqlite3.connect("/storage/resources/source/dbSTR/dbSTR.db")
    df   = pd.read_sql_query( sql  , conn)
    return df

app = dash.Dash()


app.layout = html.Div([
    html.H4('DataTable'),
    html.Label('Report type:', style={'font-weight': 'bold'}),
    dcc.Input(
        id='field-dropdown',
        type='text',
        value="MMP9"
    ),
    html.Button('Submit',id='Button-1'),
    html.Div(id='selected-indexes'),
    dcc.Graph(
        id='Main-graphic'
    ),
    dt.DataTable(
        # Initialise the rows
        rows=[{}],
        row_selectable=True,
        filterable=True,
        sortable=True,
        selected_row_indices=[],
        id='table2'
    ),
    html.Div(id='selected-indexes')
], className='container')

@app.callback(Output('table2', 'rows'), [Input('Button-1','n_clicks')],
                                      state=[State('field-dropdown', 'value')])
def update_table(numc,user_selection):
    """
    For user selections, return the relevant table
    """

    the_attrib = "gene_name"

    if user_selection[:4] == "ENSG":
        the_attrib = "gene_id"
        print(user_selection)
    else:
        the_attrib = "gene_name"

    sql2 = ("select fe.seqid Chrom,fe.featuretype Type,minmaxstart.start Begin,minmaxstart.end Ending,str.motif,str.start,str.end,avg(str.period) peri,avg(str.length) length" 
           " from" 
           " strlocmotif str," 
           " features fe," 
           " attribs at," 
           " (select seqid,min(start)-10000 start, max(end)+10000 end" 
           " from features fe," 
           " attribs at" 
           " where"  
           " at.id = fe.id" 
           " and at.value = '{}' "  
           " and at.attrib = '{}' "  
           " group by seqid) minmaxstart" 
           " where"  
           " str.chrom =  minmaxstart.seqid" 
           " and str.start >= minmaxstart.start" 
           " and str.end   <= minmaxstart.end" 
           " and at.id = fe.id" 
           " and at.value = '{}' "  
           " and at.attrib = '{}' "  
           " and fe.seqid = minmaxstart.seqid" 
           " group by fe.seqid,fe.featuretype,minmaxstart.start,minmaxstart.end,str.start,str.end,str.motif;").format(user_selection,the_attrib,user_selection,the_attrib)

    df = run_query_withparms(sql2)
    return df.to_dict('records')



#@app.callback(Output('Main-graphic','figure'),
#             [Input('table','rows')])
#def update_figure(rows):
#    dff = pd.DataFrame(rows)
#    tr1 = go.Scatter(
#    x = np.linspace(min(dff['begin']),max(dff['ending']),num=max(dff['ending'])-min(dff['begin'])),
#    y = dff['peri']
#    )
#    
#    data=[tr1]
#    
#    layout = go.Layout(
#     xaxis = dict(range=(min(dff['begin']),max(dff['ending'])))
#     )
#    fig = go.Figure(data=data, layout = layout)
#    return fig

@app.callback(Output('Main-graphic','figure'),
              [Input('table2','rows')])
def update_figure2(rows):
    dff = pd.DataFrame(rows)
    tr1 = go.Scatter(
    x = np.linspace(min(dff['start']),max(dff['end']),num=max(dff['end'])-min(dff['start'])),
    y = dff['peri']
    )

    data=[tr1]

    layout = go.Layout(
     xaxis = dict(range=(min(dff['start']),max(dff['end'])))
     )
    fig = go.Figure(data=data, layout = layout)
    return fig

   
if __name__ == '__main__':
    app.run_server(debug=True)
