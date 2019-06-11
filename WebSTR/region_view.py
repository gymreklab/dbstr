import json
import pandas as pd
import numpy as np
import plotly
import plotly.graph_objs as go
import re
from dbutils import *

MAXREGIONSIZE = 1000000

def GetRegionData(region_query, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    colpos = region_query.find(":")
    genebuf = 0.1 # increase region width by this much
    if colpos < 0: # search is by gene
        the_attrib = "gene_name"
        if region_query[:3] == "ENS":
            the_attrib = "gene_id"
            region_query = region_query.split(".")[0]
        else: 
            the_attrib = "gene_name"
        gene_query = ("select fe.seqid,min(fe.start),max(fe.end)"
                      " from features fe, newattrib at "
                      " where at.value='{}' and at.attrib='{}' and fe.id=at.id").format(region_query, the_attrib)
        gene_df = ct.execute(gene_query).fetchall()
        if len(gene_df) == 0 or None in gene_df[0]:
            chrom = None
            start = None
            end = None
        else:
            print(gene_df)
            chrom = "chr"+gene_df[0][0].replace("chr","")
            start = int(gene_df[0][1])
            end = int(gene_df[0][2])
            buf = int((end-start)*(genebuf))
            start = start-buf
            end = end+buf
    else:
        try:
            chrom = "chr"+region_query.split(":")[0].replace("chr","")
            start = int(region_query.split(":")[1].split("-")[0])
            end = int(region_query.split(":")[1].split("-")[1])
            if (end-start)>MAXREGIONSIZE:
                chrom, start, end = None, None, None
        except:
            chrom, start, end = None,None, None
    if chrom is not None:
        region_query = ("select str.chrom,str.strid,str.motif,str.start,str.end,str.period,str.length"
                        " from"
                        " strlocmotif str"
                        " where str.chrom = '{}' and str.start >= {} and str.end <= {}").format(chrom, start, end)
        df = ct.execute(region_query).fetchall()
        if len(df) == 0: return pd.DataFrame({})
        df_df = pd.DataFrame.from_records(df)
        df_df.columns = ["chrom","strid", "motif", "str.start","str.end","period","str.length"] 
        H_query = GetHCalc(df_df.strid.unique(),DbSTRPath)
        df_df["featuretype"] = "NA"
        df_df["chrom"] = df_df["chrom"].apply(lambda x: x.replace("chr",""))
        df_df["str.length"] = df_df["str.length"].round(2)
        df_df = df_df[["chrom","str.start","str.end","motif","period","str.length","strid","featuretype"]].sort_values("str.start")
        df_df.drop_duplicates(inplace=True)
    else: df_df = pd.DataFrame({})
    return df_df

def GetHvalSeqHTML(df):
    df2 = pd.DataFrame(np.array(df).reshape(-1,2), columns= list("TR"))
    t1 = df2["R"].str.split(pat=":",n= -1, expand = True)
    df2['ret'] = None
    nrows =t1.shape[0]
    delim = [', ',';']
    for i in range(nrows):
        ret = ''
        t2ta = df2.iloc[i]
        t2 = t2ta["R"].split(":")
        for j in range(len(t2)):
            t2t = t1.iloc[i,j]
            t2=re.split(r'(?:' + '|'.join(delim) + r')', t2t )
            t2num = float(t2[1])
            if t2num == 0.0: 
                thecolor = "gray"
            elif t2num <= 0.1:
                thecolor = "lightblue"
            elif t2num <= 0.5:
                thecolor = "navy"
            else:
                thecolor = "red"
            ret += '<font size="2" color="black">' + t2[0] + '</font>'
            ret += '<font size="2" color="' + thecolor + '"> ' + t2[1] + ';</font>'
        df2.iloc[i,2] = ret
    df2.columns = ["str_id","Hvals","Hhtml"]
    return df2            



def GetHCalc(strid,DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    stridt = tuple(strid)
    gquery = (" select str_id, group_concat(name || ';' || H,':') Hvals from"
              " (select name,str_id,round(1-sum(fi*fi),1) H"
              " from"
              " (select cohort_id,str_id,copies,sum(cast(nvals as FLOAT)/cast(totvals as FLOAT)) as fi"
              " from"
              " (select af.cohort_id, af.str_id, (end-start+1+af.length)/period copies,sum(nvals) nvals, tvals.totvals from"
              " allelefreq af,"
              " strlocmotif strm,"
              " (select cohort_id,str_id,sum(nvals) totvals"
              " from allelefreq"
              " where str_id in {} "
              " group by str_id,cohort_id"
              " order by str_id,cohort_id) tvals"
              " where af.str_id = strm.strid"
              " and af.str_id in {} "
              " and af.str_id = tvals.str_id"
              " and af.cohort_id = tvals.cohort_id"
              " group by af.cohort_id, af.str_id, copies)"
              " group by cohort_id,str_id,copies) sum1,"
              " COHORTS co"
              " where sum1.cohort_id = co.cohort_id"
              " group by str_id , name ) group by str_id").format(stridt,stridt)
    df = ct.execute(gquery).fetchall()
    df2 = GetHvalSeqHTML(df)
    return df2



def GetColor(period):
    colors = ["gray","red","gold","blue","purple","green"]
    return colors[int(period)-1]

def GetGenePlotlyJSON(region_data, region_query, DbSTRPath):
    chrom = region_data["chrom"].values[0].replace("chr","")

    # Get points for each STR
    trace1 = go.Scatter(
        x = (region_data["str.start"]+region_data["str.end"])/2,
        y = [0]*region_data.shape[0],
        mode="markers",
        marker=dict(size=10, color=region_data["period"].apply(lambda x: GetColor(x)), line=dict(width=2)),
        text=region_data.apply(lambda x: x["chrom"]+":"+str(x["str.start"]) + " ("+x["motif"]+")", 1),
        hoverinfo='text'
        #name=str(region_data["period"].unique())
    )

    # Draw gene info
    gene_trace, gene_shapes, numgenes = GetGeneShapes(region_data, region_query, DbSTRPath)
    plotly_data = [trace1, gene_trace]
    plotly_layout= go.Layout(
        height=300+50*numgenes,
        hovermode= 'closest',
        showlegend= False,
        #legend=dict(orientation="h"),
        shapes=gene_shapes,
        xaxis=dict(
            title="Position (chr%s)"%chrom ,
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=True,
            showticklabels=True,
            tickformat = '.0f'
        ),
        yaxis=dict(
            fixedrange=True,
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

def make_bar_trace(X, Cohort):
    return go.Bar(
        x=X['b'],
        y=X['c'],
        name=Cohort)

def GetFreqPlotlyJSON2(freq_dist):
    data1 = pd.DataFrame(np.array(freq_dist).reshape(-1,3), columns = list("abc"))
    minx = min(data1['b'])-1
    maxx = max(data1['b'])+1
    x1=data1.loc[data1['a'] == 1]
    x2=data1.loc[data1['a'] == 2]
    x3=data1.loc[data1['a'] == 3]
    x4=data1.loc[data1['a'] == 4]
    
    #for Cohort, X in data1.groupby('a'):
    trace1 = go.Bar(
        x=x1['b'],
        y=x1['c']
    )

    trace2 = go.Bar(
        x=x2['b'],
        y=x2['c'],
        xaxis='x2',
        yaxis='y2'
    )

    trace3 = go.Bar(
        x=x3['b'],
        y=x3['c'],
        xaxis='x3',
        yaxis='y3'
    )

    trace4 = go.Bar(
        x=x4['b'],
        y=x4['c'],
        xaxis='x4',
        yaxis='y4'
    )

    data = [trace1, trace2, trace3, trace4]

    layout = go.Layout(
        showlegend=False,
        xaxis=dict(
            domain=[0, 0.24],
            title="Gtex",
            range=[minx, maxx]
        ),
        yaxis=dict(
            title="Count"
        ),
        xaxis2=dict(
            domain=[0.26, 0.5],
            anchor='y2',
            title="1000 Genomes Africa",
            range=[minx, maxx]
        ),
        yaxis2=dict(
            anchor='x2'
        ),
        xaxis3=dict(
            domain=[0.51, 0.75],
            anchor='y3',
            title="1000 Genomes East Asia",
            range=[minx, maxx]
        ),
        yaxis3=dict(
            anchor='x3'
        ),
        xaxis4=dict(
            domain=[0.76, 1],
            anchor='y4',
            title="1000 Genomes Europe",
            range=[minx, maxx]
        ),
        yaxis4=dict(
            anchor='x4'
        ))

    plotly_plot_json_datab = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)
    plotly_plot_json_layoutb = json.dumps(layout, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_plot_json_datab, plotly_plot_json_layoutb

def GetFreqPlotlyJSON(freq_dist):
    data1 = pd.DataFrame(np.array(freq_dist).reshape(-1,3), columns = list("abc"))
    print(data1)
    cohort, ppx, ppy = zip(*freq_dist)
    data=[
        go.Bar(
            x=ppx,
            y=ppy
        )
    ]


    layout = go.Layout(
        xaxis=dict(
            title="Count")
        )

    plotly_plot_json_datab = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)
    plotly_plot_json_layoutb = json.dumps(layout, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_plot_json_datab, plotly_plot_json_layoutb

exon_width = 0.3
gene_width = 0.03
gene_color = "black"
def GetGeneShapes(region_data, region_query, DbSTRPath):
    ct = connect_db(DbSTRPath).cursor()
    # First, get list of genes in this region
    genes = []
    colpos = region_query.find(":")
    if colpos < 0: # search is by gene
        if region_query[:3] == "ENS":
            gene_query = ("select at.value from newattrib at where at.id='{}' and at.attrib='gene_name'").format(region_query.split(".")[0])
            gene_df = ct.execute(gene_query).fetchall()
            genes = [item[0] for item in gene_df]
        else:
            genes.append(region_query)
    else:
        chrom = "chr"+region_query.split(":")[0].replace("chr","")
        start = int(region_query.split(":")[1].split("-")[0])
        end = int(region_query.split(":")[1].split("-")[1])
        #gene_query = ("select at.value from features fe, newattrib at where fe.seqid='{}' and fe.start>={} and fe.end<={} and fe.id=at.id and at.attrib='gene_name'").format(chrom, start, end)
        gene_query = ("select value from newattrib where attrib = 'gene_name' and id in (select id from features where  seqid='{}' and start>={} and end<={} group by id) ").format(chrom, start, end)
        gene_df = ct.execute(gene_query).fetchall()
        genes = list(set([item[0] for item in gene_df]))
    shapes = []
    # Keep track of gene info
    gene_starts = []
    gene_ends = []
    gene_strands = []
    # Then, for each gene get features 
    for i in range(len(genes)):
        gene = genes[i]
        feature_query = ("select fe.id,fe.start,fe.end,fe.strand from features fe, newattrib at where at.attrib='gene_name' and at.value='{}' and fe.id=at.id").format(gene)
        feature_df = ct.execute(feature_query).fetchall()
        if len(feature_df)==0: continue
        gene_start = min([int(item[1]) for item in feature_df])
        gene_end = max([int(item[2]) for item in feature_df])
        gene_starts.append(gene_start)
        gene_ends.append(gene_end)
        gene_strands.append(feature_df[0][3])
        # First put line for whole gene
        shape = {
            "type": "rect",
            "x0": gene_start,
            "x1": gene_end,
            "y0": (i+1)-gene_width/2,
            "y1": (i+1)+gene_width/2,
            "fillcolor": gene_color,
            "line": {"width": 0}
            }
        shapes.append(shape)
        # Then put each feature (exons only)
        for f in feature_df:
            if "CDS" in f[0]: continue
            shape = {
                "type": "rect",
                "x0": int(f[1]),
                "x1": int(f[2]),
                "y0": (i+1)-exon_width/2,
                "y1": (i+1)+exon_width/2,
                "fillcolor": gene_color,
                "line": {"width": 0}
                }
            shapes.append(shape)
    trace = go.Scatter(
        x = gene_starts,
        y = [(i+1+exon_width) for i in range(len(genes))],
        mode = "text",
        hoverinfo="none",
        textposition='middle right',
        text = [GetGeneText(genes[i], gene_strands[i]) for i in range(len(genes))],
        textfont=dict(
            family='sans serif',
            size=20,
            color='black')
    )
    return trace, shapes, len(genes)

def GetGeneText(gene, strand):
    txt = "<i>"
    if strand == "+":
        txt += gene + " " + "&#8594;"
    else:
        txt += "&#8592;" + " " + gene
    txt += " "*3
    txt += "</i>"
    return txt
