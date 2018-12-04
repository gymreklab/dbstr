#!/usr/bin/env python3
"""
WebSTR database application

"""

import argparse
import dash
from flask import Flask, redirect, render_template, request, session, url_for
from dash.dependencies import Output, Input, State
from collections import deque
import pandas as pd
import numpy as npa
import json
from textwrap import dedent as d
import os

from locus_view import *
from region_view import *

#################### Make these as command line arguments ###############
global BasePath
glochrom = 1
BasePath = "/storage/resources/dbase/dbSTR/SS1/"
BasePathM = "/storage/resources/dbase/dbSTR/"

#################### Set up flask server ###############
server = Flask(__name__)
server.secret_key = 'dbSTR' 

#################### Render locus page ###############
app = dash.Dash(__name__, server=server, url_base_pathname='/dashapp')
app.config['suppress_callback_exceptions']=True
SetupDashApp(app)

@app.callback(dash.dependencies.Output('field-dropdown','value'),
              [dash.dependencies.Input('url', 'href')])
def main_display_page(href): return display_page(href)

@app.callback(Output('table2', 'rows'), [Input('field-dropdown', 'value')])
def main_update_table(user_selection): return update_table(user_selection, BasePath)

@app.callback(Output('STRtable', 'rows'), [Input('field-dropdown', 'value')])
def main_getdata(user_selection): return getdata(user_selection, BasePath)

@app.callback(Output('Main-graphic','figure'),
              [Input('table2','rows')])
def main_update_figure(rows): return update_figure(rows)
#################### Render region page ###############

@server.route('/awesome')
def awesome():
    region_query = request.args.get('query')
    region_data = GetRegionData(region_query, BasePathM)
    if region_data.shape[0] > 0:
        chrom = region_data["chrom"].values[0].replace("chr","")
        plotly_plot_json, plotly_layout_json = GetGenePlotlyJSON(region_data, region_query, chrom)
        return render_template('view2.html',table=region_data.to_records(index=False),
                               graphJSON=plotly_plot_json, layoutJSON=plotly_layout_json,
                               chrom=chrom, strids=list(region_data["strid"]))
    else:
        return render_template('view2_nolocus.html')

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
