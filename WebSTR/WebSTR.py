#!/usr/bin/env python3
"""
WebSTR database application
"""

import argparse
from flask import Flask, redirect, render_template, request, session, url_for
from collections import deque
import pandas as pd
import numpy as np
import json
from textwrap import dedent as d
import sys
import os

from locus_view import *
from region_view import *

#################### Database paths ###############
PLATFORM = "snorlax" # or AWS
if PLATFORM == "snorlax":
    BasePath = "/storage/resources/dbase/dbSTR/SS1/" # TODO this is allele freq. not used now
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

#################### Render region page ###############
@server.route('/awesome')
def awesome():
    region_queryOrg = request.args.get('query')
    region_query = region_queryOrg.upper()
    region_data = GetRegionData(region_query, DbSTRPath)
    if region_data.shape[0] > 0:
        strs_id = region_data.strid.unique()
        H_data = GetHCalc(strs_id,DbSTRPath)
        estr_data = GetestrCalc(strs_id,DbSTRPath)
        Regions_data = pd.merge(region_data, H_data, left_on='strid', right_on = 'str_id')
        Regions_data = pd.merge(Regions_data, estr_data, left_on='strid', right_on = 'str_id', how='left')
        Regions_data = Regions_data.replace(np.nan, '', regex=True)
        plotly_plot_json, plotly_layout_json = GetGenePlotlyJSON(Regions_data, region_query, DbSTRPath)
        return render_template('view2.html',table=Regions_data.to_records(index=False),
                               graphJSON=plotly_plot_json, layoutJSON=plotly_layout_json,
                               chrom=region_data["chrom"].values[0].replace("chr",""),
                               strids=list(Regions_data["strid"]))
    else:
        return render_template('view2_nolocus.html')

#################### Render locus page ###############

reffa = pyfaidx.Fasta(RefFaPath)
@server.route('/locus')
def locusview():
    str_query = request.args.get('STRID')
    chrom, start, end, seq = GetSTRInfo(str_query, DbSTRPath, reffa)
    gtex_data = GetGTExInfo(str_query, DbSTRPath)
    mut_data = GetMutInfo(str_query, DbSTRPath)
    imp_data = GetImputationInfo(str_query, DbSTRPath)
    imp_allele_data = GetImputationAlleleInfo(str_query, DbSTRPath)
    freq_dist = GetFreqSTRInfo(str_query, DbSTRPath)
    if len(mut_data) != 1: mut_data = None
    else:
        mut_data = list(mut_data[0])
        mut_data[0] = 10**mut_data[0]
    if len(imp_data) != 1: imp_data = None
    else:
        imp_data = list(imp_data[0])

    if len(gtex_data) == 0: gtex_data = None
    if len(imp_allele_data) == 0: imp_allele_data = None
    if len(freq_dist) > 0:
        plotly_plot_json_datab, plotly_plot_json_layoutb = GetFreqPlotlyJSON2(freq_dist)
        return render_template('locus.html', strid=str_query,
                           graphJSONx=plotly_plot_json_datab,graphlayoutx=plotly_plot_json_layoutb, 
                           chrom=chrom.replace("chr",""), start=start, end=end, strseq=seq,
                           estr=gtex_data, mut_data=mut_data,
                           imp_data=imp_data, imp_allele_data=imp_allele_data,freq_dist=freq_dist)

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

@server.route('/url')
def my_method():
    try:
        call_method_that_raises_exception()
    except Exception as e:
            render_template("500.html", error= str(e))

@server.errorhandler(404)
def internal_server_error(error):
    server.logger.error('Server Error: %s', (error))
    return render_template('500.htm', emsg = error), 404

@server.errorhandler(Exception)
def unhandled_exception(e):
    server.logger.error('Unhandled Exception: %s', (e))
    return render_template('500.html', emsg = e), 500

#################### Set up and run the server ###############
def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--host", help="Host to run app", type=str, default="0.0.0.0")
    parser.add_argument("--port", help="Port to run app", type=int, default=5000)
    args = parser.parse_args()
    server.run(debug=False, host=args.host, port=args.port)

if __name__ == '__main__':
    main()
