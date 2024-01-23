"""
Simple Dash app for deploying on personal computer to access the PanPPI data
For more information see https://github.com/eporetsky/PanPPI
To run the app, move the app.py and functions.py files and the tmp/ folder to
    the folder containing all the analysis results (together with the input/ 
    and output/ folders). When ready, in the terminal run: "python app.py"
Developed by Elly Poretsky on 01.22.24

Helper functions are imported by the Dash app.py file when started
"""

import os
import glob
import base64
import random
import string
import pickle
from io import BytesIO
import pandas as pd
import numpy as np
import networkx as nx
import dash
from dash import Input, Output, State, dcc, html, dash_table
import dash_bootstrap_components as dbc

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')

# A function that generates tmp pickle file with the current session_id
def tmp_write(session_id, name, content):
    with open('tmp/{}_{}.pkl'.format(session_id, name), 'wb') as handle:
        pickle.dump(content, handle, protocol=pickle.HIGHEST_PROTOCOL)

# A function that reads a tmp pickle file with the current session_id
def tmp_read(session_id, name):
    with open('tmp/{}_{}.pkl'.format(session_id, name), 'rb') as handle:
        return(pickle.load(handle))

# A function that reads the phytozome gene annotation file and gets the
# arabidopsis annotation. For multiple isoforms only the first one is returned.
# Some phytozome annotations might need to be edited to get the function to work.
def get_phytozome_annotations():  
    annot = pd.DataFrame()
    for file_name in glob.glob("output/descriptions/*"):
        tmp = pd.read_csv(file_name, sep="\t")
        #tmp = tmp.drop_duplicates(subset="locusName")
        #tmp = tmp[["locusName", "Best-hit-arabi-name", "Best-hit-arabi-defline"]]
        #tmp["genotype"] = file_name.split("/")[-1].split(".")[0]
        #tmp.columns = ["geneID", "best-arabidopsis","annotation-arabidopsis"]
        annot = pd.concat([annot, tmp])
    return(annot)

# Get a list of all genotypes, including pan and core, in a given folder
def get_all_genotypes(folder):
    genotype_list = []
    enrichment_list = glob.glob("output/{}/genomes/*".format(folder))
    enrichment_list += glob.glob("output/{}/*".format(folder))
    for fl in enrichment_list:
        genotype = fl.split("/")[-1].split(".")[0]
        genotype_list.append(genotype)
    genotype_list.remove("genomes")
    return(genotype_list)

# Return a dataframe with all the clusters (p-val<0.05 and size<50) 
# containing genes of interes. The pan-gene ID dictionary is used to
# match all possible gene IDs to pan-gene IDs and return their clusters
def get_cluster_df():
    with open('input/pangenes/maizegdb.pan.pkl', 'rb') as handle:
        pan_dict = pickle.load(handle)
    pan_df = pd.DataFrame(pd.Series(pan_dict)).reset_index()
    pan_df.columns = ["geneID", "panID"]
    cluster_df = pd.DataFrame()
    genotypes = get_all_genotypes("clusters")
    for genotype in genotypes:
        if genotype in ["pan", "core"]:
            edges_df = pd.read_csv("output/networks/{}.ppi.csv".format(genotype))
            edges_df = edges_df.drop("universality", axis=1)
            cone = pd.read_csv("output/clusters/{}.cone.csv".format(genotype))
        else:
            edges_df = pd.read_csv("output/networks/genomes/{}.ppi.csv".format(genotype))
            edges_df = edges_df.drop("weight", axis=1)
            cone = pd.read_csv("output/clusters/genomes/{}.cone.csv".format(genotype))
        
        cone = cone[cone["P-value"]<0.05]
        cone = cone[cone["Size"]<50]
        tmp = {"clusterID": [], "geneID": []}
        for ix, row in cone.iterrows():
            clusterID = "{}_{}".format(genotype, row["Cluster"])
            for geneID in row["Members"].split(" "):
                tmp["clusterID"].append(clusterID)
                tmp["geneID"].append(geneID)
        cluster_df = pd.concat([cluster_df, pd.DataFrame(tmp)])
    return(cluster_df)

# Return the GO enrichment results combined dataframe of all selected clusters
def get_selected_enrichments(candidate_clusters_list, genotype_list):
    enrichment_list = glob.glob("output/enrichment/genomes/*")
    enrichment_list += glob.glob("output/enrichment/*")
    enrichment_list.remove('output/enrichment/genomes')
    df_enrichment = pd.DataFrame()
    for fl in enrichment_list:
        genotype = fl.split("/")[-1].split(".")[0]
        if genotype in genotype_list:
            tmp = pd.read_csv(fl, sep="\t")
            tmp["genotype"] = genotype
            df_enrichment = pd.concat([df_enrichment, tmp])
    df_enrichment = df_enrichment[df_enrichment["cluster"].isin(candidate_clusters_list)]
    return(df_enrichment)

# Main function for generating the content of the enrichment tab
# Returns a dropdown menu for filtering on BP, MF and CC
# Returns the dash DataTable object that will be shown on the page
def tab_content_enrichment(df, genotype_list):
    df = df[df["genotype"].isin(genotype_list)]
    df = df.drop("genotype", axis=1)
    df["p-val"] = (-np.log10(df["p-val"])).round(decimals=2)
    df["FC"] = df["FC"].round(decimals=2)
    output_table = dash_table.DataTable(
        id="state_enrichment_table",
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
        style_table={'minWidth': '100%'},
        #fixed_columns={'headers': True, 'data': 1},
        style_cell={'textAlign': 'left',
                    #'overflow': 'hidden',
                    #'textOverflow': 'ellipsis',
                    },
        editable=False,
        row_deletable=False,
        export_format="csv",)
    filter_list = df["type"].unique().tolist()
    return(dcc.Dropdown(
                id='enrichment_filter_dropdown',
                options=[{'label': genotype, 'value': genotype} for genotype in filter_list],
                value=filter_list,
                multi=True,
                searchable=False),
           output_table)

def get_candidate_pan_genes(candidate_genes):
    with open('input/pangenes/maizegdb.pan.pkl', 'rb') as handle:
        pan_dict = pickle.load(handle)
    candidate_pan_genes = []
    for gene in candidate_genes:
        try:
            candidate_pan_genes.append(pan_dict[gene])
        except:
            None
    return(candidate_pan_genes)

def fig_to_uri(in_fig, close_all=True, **save_args):
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)

def tab_content_network(session_id,
                        n_clicks_ids=0, 
                        n_clicks_coexp=0):
    candidate_clusters_df = tmp_read(session_id, "clusters")
    if len(candidate_clusters_df)==0:
        return(html.P("No clusters have been identified with the given set of protein IDs."))
    candidate_clusters_list = candidate_clusters_df["clusterID"].unique()
    
    try:
        selected_cluster = tmp_read(session_id, "selected_cluster")
    except:
        selected_cluster = candidate_clusters_list[0]
        tmp_write(session_id, "selected_cluster", selected_cluster)
    
    genotype = selected_cluster.split("_")[0]
    if genotype in ["pan", "core"]:
        edge_df = pd.read_csv("output/coexpression/{}.coexp.csv".format(genotype))
    else:
        edge_df = pd.read_csv("output/coexpression/genomes/{}.coexp.csv".format(genotype))
    edge_df = edge_df[edge_df["cluster"]==selected_cluster]
    edge_df = edge_df.drop("cluster", axis=1)
    
    show_labels = False if n_clicks_ids % 2 == 0 else True
    show_coexp = True if n_clicks_coexp % 2 == 0 else False

    return_buttons = [
            dbc.Button("Toggle protein IDs", id="btn_network_ids", color="success", n_clicks=n_clicks_ids),
            #dbc.Button("Toggle coexpression", id="btn_network_coexpression", color="warning", n_clicks=n_clicks_coexp),
        ]
    
    return_dropdown = dcc.Dropdown(
                id='dropdown_clusters',
                options=[{'label': cluster, 'value': cluster} for cluster in candidate_clusters_list],
                value=selected_cluster,
                multi=False,
                searchable=True),

    #
    G = nx.from_pandas_edgelist(edge_df, 'source', 'target', 'weight')
    # Use tmp pickle file to make the gene annotation file 
    tmp_write(session_id, "selected_ids", list(G.nodes))
    pos = nx.circular_layout(G)
    edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())
    # blue: both, red: coexpress, black: PPI
    weight_colors = {"PPI": "blue", "Coexp": "red", "Both": "black"}
    weights = [weight_colors[w] for w in weights]

    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10), constrained_layout=True)
    nx.draw(G, pos, node_color='black', edge_color=weights, ax=ax1)
    # Manually add labels with a white rectangle because otherwise
    # it is hard to read the labels infront of all edges 
    if show_labels:
        # Add a white rectangle behind the node label
        ax = plt.gca()
        for label, (x, y) in pos.items():
            text = ax.text(x, y, label, ha='center', va='center', zorder=10)
            text.set_bbox(dict(facecolor='white', edgecolor='none', pad=0.2))
    fig_path = 'tmp/{}_network.png'.format(session_id)
    fig.savefig(fig_path)
    
    return(html.Div([
                dbc.Row(return_buttons, 
                        justify="center", align="center",
                        style={"width":"60%"}),
                dbc.Row(return_dropdown, 
                        justify="center", align="center",
                        style={"width":"60%"}),
                dbc.Row([html.Img(src = fig_to_uri(fig))],
                        justify="center", align="center",
                        style={"height": "60%", "width":"60%"})
        ])
    )

def tab_content_annotation(session_id, df):
    try:
        gene_list = tmp_read(session_id, "selected_ids")
    except:
        return(html.P("No cluster has been selected"))
    df = df[df["geneID"].isin(gene_list)]
    output_table = dash_table.DataTable(
        id="state_annotation_table",
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
        style_table={'minWidth': '100%'},
        #fixed_columns={'headers': True, 'data': 1},
        style_cell={'textAlign': 'left',
                    #'overflow': 'hidden',
                    #'textOverflow': 'ellipsis',
                    },
        editable=False,
        row_deletable=False,
        export_format="csv",)
    return(output_table)
