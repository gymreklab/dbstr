import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
from dash.exceptions import PreventUpdate
import pandas as pd

from dbutils import *

def SetupDashApp(app):
    app.css.append_css({'external_url': 'https://codepen.io/amyoshino/pen/jzXypZ.css'})

    dcc._css_dist[0]['relative_package_path'].append('datatable.css')
    
    app.layout =  html.Div([
        html.Link(
            rel='stylesheet',
            href='/static/stylesheet.css'),
        html.Div(className='row',children=[
            html.Div([
                html.Nav(className = "navbar", children=[
                    html.A('Main', className="navbar", href='/dbSTR'),
                    html.A('FAQ',  className="navbar", href='/faq') 
                ]),
                html.Label('Enter Gene:', style={'font-weight': 'bold'}),
                dcc.Location(id='url', refresh=False),
                dcc.Input(
                    id='field-dropdown',
                    type='text',
                    value=[{}]
                )
            ]),
            
            html.Div([
                html.Div(id='selected-indexes'),
                dcc.Graph(
                    id='Main-graphic',
                    style={'height': 400, 'width': 600}),
                dt.DataTable(
                    # Initialise the rows
                    rows=[{}],
                    row_selectable=True,
                    filterable=True,
                    sortable=True,
                    selected_row_indices=[],
                    id='STRtable'
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
            ],style={'height': '25%','display': 'block'})
        ])
    ])
    

def display_page(href):
    # this is called every page load and every URL change
    # you can update your components with this URL in herae
    global glochrom
    if href is None:
       raise PreventUpdate
    findq1,findq2,findq3 = href.split("?")
    findq1e,findq2e = findq2.split("=")
    findq1q,findq2q = findq3.split("=")
    glochrom = findq2q
    #t1 = getdata(findq2q,findq2e)
    return findq2e

def update_table(user_selection, BasePath):
    """
    For user selections, return the relevant table
    """

    #the_attrib = "gene_name"

    #if user_selection[:4] == "ENSG":
    #    the_attrib = "gene_id"
    #    print(user_selection)
    #else:
    #    the_attrib = "gene_name"

    sql2 = ("select * from Sample_szes" 
            " where str_id = '{}' "
            " ; ").format(user_selection)

    df = run_query_withparms(sql2, BasePath, glochrom) # TODO remove glochrom
    return df.to_dict('records')

def getdata(user_selection, BasePath):
    global glochrom
    connt = sqlite3.connect(BasePath + "dbSTR" + str(glochrom) + ".db")

    sql3 = ("SELECT lent, "
           " Case mult when 2 then max(lent) else 0 end, "
           " CASE mult when 1 then max(lent) else 0 end, "
           " count(sample) strcnt from ( "
           " select base.*,hzyg.sample,' ' alt,length(base.ref) lent,2 mult from"
           " vcfBase as base,"
           " vcfhomozyg as hzyg"
           " where"
           "     base.str_id=hzyg.str_id"
           " and hzyg.str_id = '{}'"
           " UNION ALL"
           " select base.*,GT.sample,alt.alt,length(alt.alt) lent,1 mult from"
           " vcfBase as base,"
           " vcfAlt as alt,"
           " altGT as gt"
           " where "
           "     base.pos=alt.pos"
           " and base.str_id=alt.str_id"
           " and base.str_id = gt.str_id"
           " and base.str_id = '{}'"
           " and alt.str_id = gt.str_id"
           " and gt.altref_gt = alt.altorder)"
           " group by lent;").format(user_selection,user_selection)

    df   = pd.read_sql_query( sql3  , connt)
    return df.to_dict('records')

def update_figure(rows):
    dff = pd.DataFrame(rows)

    genes = dff
    traceg = go.Histogram(
             x = genes['lent'],
             text = genes['alt'])

    layout = go.Layout(
         xaxis=dict(
              range= [min(genes['lent']),max(genes['lent'])],
              dtick=1
         )
    )


    data =[traceg]
    fig = go.Figure(data=data, layout = layout)
    return fig
