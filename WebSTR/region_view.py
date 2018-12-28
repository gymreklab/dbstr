import json
import pandas as pd
import plotly
import plotly.graph_objs as go
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
            chrom, start, end = None, None, None
    if chrom is not None:
        region_query = ("select str.chrom,str.strid,str.motif,str.start,str.end,str.period,str.length"
                        " from"
                        " strlocmotif str"
                        " where str.chrom = '{}' and str.start >= {} and str.end <= {}").format(chrom, start, end)
        df = ct.execute(region_query).fetchall()
        if len(df) == 0: return pd.DataFrame({})
        df_df = pd.DataFrame.from_records(df)
        df_df.columns = ["chrom","strid", "motif", "str.start","str.end","period","str.length"] 
        df_df["featuretype"] = "NA"
        df_df["chrom"] = df_df["chrom"].apply(lambda x: x.replace("chr",""))
        df_df["str.length"] = df_df["str.length"].round(2)
        df_df = df_df[["chrom","str.start","str.end","motif","period","str.length","strid","featuretype"]].sort_values("str.start")
        df_df.drop_duplicates(inplace=True)
    else: df_df = pd.DataFrame({})
    return df_df

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
        height=250+50*numgenes,
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
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            showticklabels=False
        )
        #annotations = [dict( text= "Setting the Custom Title Position",
        #                    showarrow = False,
        #                    x = 0,
        #                    y = -2,
        #                    font = dict(size=10)
        #)]
    )

    plotly_plot_json = json.dumps(plotly_data, cls=plotly.utils.PlotlyJSONEncoder) 
    plotly_layout_json = json.dumps(plotly_layout, cls=plotly.utils.PlotlyJSONEncoder)
    return plotly_plot_json, plotly_layout_json

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
            gene_query = ("select at.value from newattrib at where at.id='{}' and at.attrib='gene_name'").format(region_query)
            gene_df = ct.execute(gene_query).fetchall()
            genes = [item[0] for item in gene_df]
        else:
            genes.append(region_query)
    else:
        chrom = "chr"+region_query.split(":")[0].replace("chr","")
        start = int(region_query.split(":")[1].split("-")[0])
        end = int(region_query.split(":")[1].split("-")[1])
        gene_query = ("select at.value from features fe, newattrib at where fe.seqid='{}' and fe.start>={} and fe.end<={} and fe.id=at.id and at.attrib='gene_name'").format(chrom, start, end)
        gene_df = ct.execute(gene_query).fetchall()
        genes = list(set([item[0] for item in gene_df]))
    shapes = []
    # Keep track of gene info
    gene_starts = []
    gene_ends = []
    gene_strands = []
    # Then, for each gene get features 
    print(genes)
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
        y = [(i+1) for i in range(len(genes))],
        mode = "text",
        hoverinfo="none",
        textposition='middle left',
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
