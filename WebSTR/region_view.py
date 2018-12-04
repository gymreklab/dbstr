import json
import pandas as pd
import plotly
import plotly.graph_objs as go
from dbutils import *

def GetRegionData(region_query, BasePathM):
    ct = connect_db(BasePathM).cursor()
    colpos = region_query.find(":")
    genebuf = 5000
    if colpos < 0: # search is by gene
        the_attrib = "gene_name"
        if region_query[:3] == "ENS":
            the_attrib = "gene_id"
        else: 
            the_attrib = "gene_name"
        gene_query = ("select fe.seqid,min(fe.start)-{},min(fe.end)+{}"
                      " from features fe, newattrib at "
                      " where at.value='{}' and at.attrib='{}' and fe.id=at.id").format(genebuf, genebuf, region_query, the_attrib)
        gene_df = ct.execute(gene_query).fetchall()
        if len(gene_df) == 0 or None in gene_df[0]:
            chrom = None
            start = None
            end = None
        else:
            chrom = "chr"+gene_df[0][0].replace("chr","")
            start = int(gene_df[0][1])
            end = int(gene_df[0][2])
    else:
        chrom = "chr"+region_query.split(":")[0].replace("chr","")
        start = int(region_query.split(":")[1].split("-")[0])
        end = int(region_query.split(":")[1].split("-")[1])

    if chrom is not None:
        region_query = ("select str.chrom,str.strid,str.motif,str.start,str.end,str.period,str.length"
                        " from"
                        " strlocmotif str"
                        " where str.chrom = '{}' and str.start >= {} and str.end <= {}").format(chrom, start, end)
        df = ct.execute(region_query).fetchall()
        if len(df) == 0: return pd.DataFrame({})
        df_df = pd.DataFrame.from_records(df)
        df_df.columns = ["chrom","strid", "motif", "str.start","str.end","period","str.length"] 
        df_df["featuretype"] = "NA" # TODO set
        df_df["chrom"] = df_df["chrom"].apply(lambda x: x.replace("chr",""))
        df_df["str.length"] = df_df["str.length"].round(2)
        df_df = df_df[["chrom","str.start","str.end","motif","period","str.length","strid","featuretype"]].sort_values("str.start")
        df_df.drop_duplicates(inplace=True)
    else: df_df = pd.DataFrame({})
    return df_df

def GetColor(period):
    colors = ["gray","red","gold","blue","purple","green"]
    return colors[int(period)-1]

def GetGenePlotlyJSON(region_data, region_query, chrom):
    trace1 = go.Scatter(
        x = (region_data["str.start"]+region_data["str.end"])/2,
        y = [0]*region_data.shape[0],
        mode="markers",
        marker=dict(size=10, color=region_data["period"].apply(lambda x: GetColor(x)), line=dict(width=2)),
        text=region_data.apply(lambda x: x["chrom"]+":"+str(x["str.start"]) + " ("+x["motif"]+")", 1),
        hoverinfo='text'
    )

    plotly_data = [trace1]
    plotly_layout= go.Layout(
        title= region_query,
        hovermode= 'closest',
        showlegend= False,
        xaxis=dict(
            title="Position (chr%s)"%chrom,
            autorange=True,
            showgrid=True,
            zeroline=False,
            showline=True,
            showticklabels=True,
            tickformat = '.0f'
        ),
        yaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            showticklabels=False
        )
    )
    plotly_plot_json = json.dumps(plotly_data, cls=plotly.utils.PlotlyJSONEncoder) 
    plotly_layout_json = json.dumps(plotly_layout, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_plot_json, plotly_layout_json

