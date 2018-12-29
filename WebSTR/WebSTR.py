#!/usr/bin/env python3
"""
WebSTR database application

"""

import argparse
#import dash
from flask import Flask, redirect, render_template, request, session, url_for
#from dash.dependencies import Output, Input, State
from collections import deque
import pandas as pd
import numpy as npa
import json
from textwrap import dedent as d
import sys
import os

#from locus_view_dash import *
from locus_view import *
from region_view import *

#################### Database paths ###############
PLATFORM = "snorlax" # or AWS
if PLATFORM == "snorlax":
    BasePath = "/storage/resources/dbase/dbSTR/SS1/"
    DbSTRPath = "/storage/resources/dbase/dbSTR/"
    RefFaPath = "/storage/resources/dbase/human/hg19/hg19.fa"
elif PLATFORM == "AWS":
    BasePath = ""
    DbSTRPath = ""
    RefFaPath = "" # TODO
else:
    sys.stderr.write("Could not locate database files\n")
    sys.exit(1)

#################### Set up flask server ###############
server = Flask(__name__)
server.secret_key = 'dbSTR' 

#################### Render locus page ###############
#app = dash.Dash(__name__, server=server, url_base_pathname='/dashapp')
#app.config['suppress_callback_exceptions']=True
#SetupDashApp(app)

#@app.callback(dash.dependencies.Output('field-dropdown','value'),
#              [dash.dependencies.Input('url', 'href')])
#def main_display_page(href): return display_page(href)

#@app.callback(Output('table2', 'rows'), [Input('field-dropdown', 'value')])
#def main_update_table(user_selection): return update_table(user_selection, BasePath)

#@app.callback(Output('STRtable', 'rows'), [Input('field-dropdown', 'value')])
#def main_getdata(user_selection): return getdata(user_selection, BasePath)

#@app.callback(Output('Main-graphic','figure'),
#              [Input('table2','rows')])
#def main_update_figure(rows): return update_figure(rows)

#################### Render region page ###############

@server.route('/awesome')
def awesome():
    region_query = request.args.get('query')
    region_data = GetRegionData(region_query, DbSTRPath)
    if region_data.shape[0] > 0:
        plotly_plot_json, plotly_layout_json = GetGenePlotlyJSON(region_data, region_query, DbSTRPath)
        return render_template('view2.html',table=region_data.to_records(index=False),
                               graphJSON=plotly_plot_json, layoutJSON=plotly_layout_json,
                               chrom=region_data["chrom"].values[0].replace("chr",""),
                               strids=list(region_data["strid"]))
    else:
        return render_template('view2_nolocus.html')

reffa = pyfaidx.Fasta(RefFaPath)
@server.route('/locus')
def locusview():
    str_query = request.args.get('STRID')
    chrom, start, end, seq = GetSTRInfo(str_query, DbSTRPath, reffa)
    gtex_data = GetGTExInfo(str_query, DbSTRPath)
    mut_data = GetMutInfo(str_query, DbSTRPath)
    imp_data = GetImputationInfo(str_query, DbSTRPath)
    imp_allele_data = GetImputationAlleleInfo(str_query, DbSTRPath)
    if len(mut_data) != 1: mut_data = None
    else:
        mut_data = list(mut_data[0])
        mut_data[0] = 10**mut_data[0]
    if len(imp_data) != 1: imp_data = None
    else:
        imp_data = list(imp_data[0])

    if len(gtex_data) == 0: gtex_data = None
    if len(imp_allele_data) == 0: imp_allele_data = None
    return render_template('locus.html', strid=str_query,
                           chrom=chrom.replace("chr",""), start=start, end=end, strseq=seq,
                           estr=gtex_data, mut_data=mut_data,
                           imp_data=imp_data, imp_allele_data=imp_allele_data)

#################### Render HTML pages ###############
@server.route('/')
@server.route('/dbSTR')
def dbSTRHome():
    return render_template('homepage.html')

@server.route('/faq')
def dbSTRFAQ():
    return render_template("faq.html")

@server.route('/contact')
def dbSTRContact():
    return render_template("contact.html")

@server.route('/about')
def dbSTRAbout():
    return render_template("about.html")

@server.route('/downloads')
def dbSTRDownloads():
    return render_template("downloads.html")

@server.route('/terms')
def dbSTRTerms():
    return render_template("terms.html")

#################### Set up and run the server ###############
def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--host", help="Host to run app", type=str, default="0.0.0.0")
    parser.add_argument("--port", help="Port to run app", type=int, default=5000)
    args = parser.parse_args()
    server.run(debug=True, host=args.host, port=args.port)

if __name__ == '__main__':
    main()
