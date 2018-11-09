import dash
import dash_core_components as dcc
import dash_html_components as html
from flask import Flask, render_template, request, session 
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
import dash_table_experiments as dt
import plotly
import plotly.graph_objs as go
from collections import deque
import pandas as pd
import sqlite3
import numpy as npa
import json
from textwrap import dedent as d
import os

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}


server = Flask(__name__)
server.secret_key = 'dbSTR' 
global BasePath
glochrom = 1
BasePath = "/storage/resources/dbase/dbSTR/SS1/"

app = dash.Dash(__name__, server=server, url_base_pathname='/dashapp')

def run_query_withparms(sql):
    print(glochrom)
    conn = sqlite3.connect(BasePath + "dbSTR" + glochrom + ".db")
    df   = pd.read_sql_query( sql  , conn)
    return df

def run_query_withparms2(sql):
    conn = sqlite3.connect(BasePath + "dbSTR.db")
    df   = pd.read_sql_query( sql  , conn)
    return df

def run_query_withparmsi(sql,value):
    conn = sqlite3.connect(BasePath + "dbSTR.db")
    cur = conn.cursor()
    cur.execute(sql,value)
    conn.close()
    return cur


def run_query_withparmstr(sql):
    conn = sqlite3.connect(BasePath + "dbSTR.db")
    c = conn.cursor()
    chr = c.execute(sql).fetchall()
    return chr[0][0]

def getdata(chrom,sqlip):
    connt = sqlite3.connect(BasePath + "dbSTR" + chrom + ".db")
    ct = connt.cursor()

    sqltemp = ("DROP TABLE IF EXISTS SMTBL;")
    chrd = ct.execute(sqltemp).fetchall()

    sql3 = (" CREATE TEMPORARY TABLE SMTBL AS"
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
            " and gt.altref_gt = alt.altorder;").format(sqlip,sqlip)
    sqlstbl = ct.execute(sql3).fetchall()

    sqlt = ("SELECT "
           " Case mult when 2 then max(lent) else 0 end, "
           " CASE mult when 1 then max(lent) else 0 end, "
           " count(sample) strcnt, mult, totcnt from SMTBL ss,"
           " (select pos,count(*) totcnt from SMTBL group by pos) as tcnt"
           " where tcnt.pos = ss.pos "
           " group by lent ")
    sumtbl = ct.execute(sqlt).fetchall()
    dfsum = pd.DataFrame.from_records(sumtbl,columns=['lenRef','lenAlt','Samp Cnt','Mult','total cnt'])
    print(dfsum)


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


@app.callback(dash.dependencies.Output('field-dropdown','value'),
              [dash.dependencies.Input('url', 'href')])
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

@app.callback(Output('table2', 'rows'), [Input('field-dropdown', 'value')])
def update_table(user_selection):
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

    df = run_query_withparms(sql2)

    return df.to_dict('records')

@app.callback(Output('STRtable', 'rows'), [Input('field-dropdown', 'value')])
def getdata2(user_selection):
    global glochrom
    connt = sqlite3.connect("dbSTR" + glochrom + ".db")

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



@app.callback(Output('Main-graphic','figure'),
              [Input('table2','rows')])
def update_figure2(rows):
    dff = pd.DataFrame(rows)

    genes = dff
    traceg = go.Histogram(
             x = genes['length'],
             text = genes['alt'])

    layout = go.Layout(
         xaxis=dict(
              range= [min(genes['length']),max(genes['length'])],
              dtick=1
         )
    )


    data =[traceg]
    fig = go.Figure(data=data, layout = layout)
    return fig

@server.route('/')
def index():
    return '''
<html>
<div>
    <h1>Flask App</h1>
</div>
</html>
'''

@server.route('/dbSTR')
def test1():
    db = get_db()
    sqlip = 'select * from COHORTS where active = "Y";'
    print(sqlip)
    ct = db.cursor()
    df = ct.execute(sqlip).fetchall()
    print(df)
    return render_template('homepage.html', option_list = df)

@server.route('/faq')
def faq1():
    return '''
<html>
<div>
    <h1>Flask App faqs... test 1</h1>
</div>
</html>
'''

def connect_db():
    """
    Connects to the specific database.
    """
    tfile = BasePath + "dbSTR.db"
    print(tfile)
    conn = sqlite3.connect(tfile)
    return conn

def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    db_conn = connect_db()
    return db_conn


@server.route('/awesome')
def awesome():
    db = get_db()
    sqlip = request.args.get('query')
    print(sqlip)

    connt = sqlite3.connect(BasePath + "dbSTR.db")
    ct = connt.cursor()

    the_attrib = "gene_name"

    if sqlip[:3] == "ENS":
        the_attrib = "gene_id"
        print(sqlip)
    else:
        the_attrib = "gene_name"

    sql2 = ("select fe.seqid,fe.featuretype,minmaxstart.start begin,minmaxstart.end ending,str.strid,str.motif,str.start,str.end,substr(str.chrom,4,length(str.chrom)),avg(str.period) peri,avg(str.length)"
    " from"
    " strlocmotif str,"
    " features fe,"
    " newattrib at,"
    " (select seqid,min(start)-10000 start, max(end)+10000 end"
    " from features fe,"
    " newattrib at"
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
    " group by fe.seqid,fe.featuretype,minmaxstart.start,minmaxstart.end,str.start,str.end,str.motif,str.strid,str.chrom;").format(sqlip,the_attrib,sqlip,the_attrib)
    
    df = ct.execute(sql2).fetchall()
    
    df_df = pd.DataFrame.from_records(df)

    trace1 = go.Scatter(
       x = npa.linspace(1,len(df),len(df)),
       y = df_df[7]
    )

    data = [trace1]
    #grappJSON = json.dumps(result,
    #                       default=lambda df: json.loads(df.to_json()))
    #mytest = json.loads(grappJSON)
    ids = range(1,len(df),1)
    mytest = json.dumps(data,cls=plotly.utils.PlotlyJSONEncoder)
    return render_template('view2.html',table=df,
                           graphJSON=mytest )

if __name__ == '__main__':
    server.run(debug=True)
