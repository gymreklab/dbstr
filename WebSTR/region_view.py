import json
import pandas as pd
import plotly
import plotly.graph_objs as go
from dbutils import *

def GetRegionData(region_query, BasePathM):
    ct = connect_db(BasePathM).cursor()
    the_attrib = "gene_name"
    if region_query[:3] == "ENS":
        the_attrib = "gene_id"
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
    " group by fe.seqid,fe.featuretype,minmaxstart.start,minmaxstart.end,str.start,str.end,str.motif,str.strid,str.chrom;").format(region_query,the_attrib,region_query,the_attrib)

    df = ct.execute(sql2).fetchall()
    df_df = pd.DataFrame.from_records(df)
    df_df.columns = ["chrom","featuretype","gene.start","gene.end", "strid", "motif", "str.start","str.end","chrom2","period","str.length"]
    df_df["featuretype"] = "NA" # TODO set
    df_df["chrom"] = df_df["chrom"].apply(lambda x: x.replace("chr",""))
    df_df["str.length"] = df_df["str.length"].round(2)
    df_df = df_df[["chrom","str.start","str.end","motif","period","str.length","strid","featuretype"]].sort_values("str.start")
    df_df.drop_duplicates(inplace=True)
    return df_df

def GetColor(period):
    colors = ["gray","red","gold","blue","purple","green"]
    return colors[int(period)-1]

def GetGenePlotlyJSON(region_data):
    trace1 = go.Scatter(
        x = (region_data["str.start"]+region_data["str.end"])/2,
        y = [0]*region_data.shape[0],
        mode="markers",
        marker=dict(size=10, color=region_data["period"].apply(lambda x: GetColor(x)), line=dict(width=2)),
        text=region_data.apply(lambda x: x["chrom"]+":"+str(x["str.start"]), 1)
    )
    plotly_data = [trace1]
    plotly_json = json.dumps(plotly_data, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_json

