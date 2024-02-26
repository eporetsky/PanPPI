"""
Simple Dash app for deploying on personal computer to access the PanPPI data
For more information see https://github.com/eporetsky/PanPPI
To run the app, move the app.py and functions.py files and the tmp/ folder to
    the folder containing all the analysis results (together with the input/ 
    and output/ folders). When ready, in the terminal run: "python app.py"
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

from functions import *

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP],
                suppress_callback_exceptions=True,)

# the style arguments for the sidebar. We use position:fixed and a fixed width
SIDEBAR_STYLE = {
    "position": "fixed",
    "margin-left": "20px",
    "padding": "10px 10px",
    "background-color": "#f8f9fa",
}

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "margin-left": "5px",
    "margin-right": "5px",
    "padding": "10px 10px",
}

# Using temporary random session_id Delete tmp sessions when running as localhost 
session_id = ''.join(random.SystemRandom().choice(string.ascii_letters + string.digits) for _ in range(10))
tmp_files = glob.glob('tmp/*')
for tmp in tmp_files:
    os.remove(tmp)

input_gene_list = "Zm00001eb008310\nZm00001eb014850\nZm00001eb051380\nZm00001eb059560\nZm00001eb172450\nZm00001eb103370\nZm00001eb103490\nZm00001eb109790\nZm00001eb122110\nZm00001eb143110\nZm00001eb155150\nZm00001eb258090\nZm00001eb259240\nZm00001eb319590\nZm00001eb352310\nZm00001eb265010\nZm00001eb297320\nZm00001eb414790\nZm00001eb433410"
input_go_list = "GO:0048574\nGO:2000028\nGO:0048573\nGO:0048578\nGO:0048577\nGO:0048579\nGO:0048587\nGO:0048586"

with open("pickles/sim_connected_components.pkl", 'rb') as handle:
        sim_cc_list = pickle.load(handle)
cluster_df = pd.read_pickle("pickles/cluster_df.pkl")
annots_df = pd.read_pickle("pickles/annots_df_df.pkl")

def get_sim_clust_list(selected_cluster):
    # If selected cluster is not in similarity connected component return itself
    sim_cluster_list = [selected_cluster]
    # Try to find similarity connected component for selected cluster
    for sim_cc in sim_cc_list:
        if selected_cluster in sim_cc:
            sim_cluster_list = sim_cc
    return sim_cluster_list

sidebar = html.Div(
    [
        html.H2("Menu", className="display-4"),
        html.Hr(),
        html.P("Parse the pan-interactome data using either a list of gene of interests or a selected cluster. Press the relevant Submit button to parse the data, overriding the previous results.", className="lead"),
        
        html.Hr(),
        html.P("1. Submit either genes or GO terms of interest to get relevant clusters", className="lead"),
        html.P("Submit genes of interest, one gene per row", className="lead"),
        dcc.Textarea(id="input_genes", 
                  value=input_gene_list),
        dbc.Button("Submit", id="btn_input_genes", color="primary", className="me-1"),

        html.P("Submit GO terms of interest, one term per row", className="lead"),
        dcc.Textarea(id="input_GOs", 
                  value=input_go_list),
        dbc.Button("Submit", id="btn_input_GOs", color="primary", className="me-1"),

        html.Hr(),
        html.P("2. Select which networks to include in the results", className="lead"),
        #dbc.Input(id="input_geneset", placeholder="Type something...", type="text"),
        html.Div(id="div_dropdown_genotypes"),
        html.P(""),
        dbc.Button("Select all", id="button_select_all", color="primary", className="me-1"),

    ],
    
    style={"margin-left": "10px",
           "padding": "10px 10px",
           "background-color": "#f8f9fa",},
)

tabs = html.Div(
    [
        dbc.Tabs(
            [
                dbc.Tab(label="GO Enrichment", tab_id="tab-enrichment"),
                dbc.Tab(label="Select Cluster", tab_id="tab-network"),
                dbc.Tab(label="Cluster Genes", tab_id="tab-annotation"),
                dbc.Tab(label="Similar Clusters", tab_id="tab-enrichment-sim"),
            ],
            id="tabs",
            active_tab="tab-enrichment",
        ),
        html.Div(id="content",
                 className="g-0",),
    ],
    style={"margin-left": "10px",
           "padding": "10px 10px",},
)

app.layout = html.Div([
    dbc.Row([
        dbc.Col(sidebar,
                 width=2),
        dbc.Col(tabs,
                 width={"size":10}, ),
        ]),
    ],
    className="g-0", 
)

### Tab 1: Submit candidate genes, choose genotypes and find relevant clusters
@app.callback(Output("div_dropdown_genotypes", "children"),
              Input("button_select_all", "n_clicks"),
              )
def switch_tab(button):
            genotype_list = get_all_genotypes("enrichment")
            return(dcc.Dropdown(
                id='dropdown_genotypes',
                options=[{'label': genotype, 'value': genotype} for genotype in genotype_list],
                value=genotype_list,
                multi=True,
                searchable=True),
            )


@app.callback(Output("content", "children", allow_duplicate=True),
              State("dropdown_genotypes", "value"),
              Input("dropdown_clusters", "value"),
              Input("dropdown_clusters_similairity", "value"),
              prevent_initial_call=True)
def btn_network_tab(genotype_list, selected_cluster, selected_clusters_similairity):
    context = [p["prop_id"] for p in dash.callback_context.triggered][0]
    
    # Don't allow not selecting any cluster (aka pressing the x button in the box)
    # Check if there's an option to prevent it in Dash
    if selected_cluster == None:
        selected_cluster = tmp_read(session_id, "selected_cluster")
    if selected_clusters_similairity == None:
        selected_clusters_similairity = selected_cluster
    
    # Update the selected_cluster temp variable and the similarity enrichment table
    if context == "dropdown_clusters.value":
        tmp_write(session_id, "selected_cluster", selected_cluster)
        
        sim_clust_list = get_sim_clust_list(selected_cluster)
        tmp_write(session_id, "sim_clust_list", sim_clust_list)

        all_genotype_list = get_all_genotypes("enrichment")
        selected_enrichment = get_selected_enrichments(sim_clust_list, None, all_genotype_list)
        content_enrichment = tab_content_enrichment(selected_enrichment, all_genotype_list)
        tmp_write(session_id, "enrichment_similar", content_enrichment[-1])
        #tmp_write(session_id, "sim_clust_list", content_enrichment[-1])
        
    elif context == "dropdown_clusters_similairity.value":
        tmp_write(session_id, "selected_cluster", selected_clusters_similairity)
    
    content_network = tab_content_network(session_id=session_id,
                                          n_clicks_ids=0,
                                          )
    return(content_network) 

@app.callback(Output("content", "children"),
              State("input_genes", "value"),
              State("dropdown_genotypes", "value"),
              Input("tabs", "active_tab"),
              Input("btn_input_genes", "n_clicks"),
              State("input_GOs", "value"),
              Input("btn_input_GOs", "n_clicks"),
              prevent_initial_call=True)
def switch_tab(gene_list, genotype_list, at, button, GO_list, button_GOs):
    context = [p["prop_id"] for p in dash.callback_context.triggered][0]
    
    # Check if the "Submit" button on the first tab was pressed to start process
    if context == "btn_input_genes.n_clicks":
        # Split the list of gene candidates coming from different rows
        candidate_genes = gene_list.split("\n")
        # Add to the list of gene candidates the list of associcated pan-gene IDs
        candidate_pangenes = get_candidate_pan_genes(candidate_genes, False)
        candidate_genes += candidate_pangenes.keys()
        # Get all the IDs of clusters associated with any of the candidate genes
        candidate_clusters_df = cluster_df[cluster_df["geneID"].isin(candidate_genes)]
        # Add the relevant 
        cluster_genes_dict = {}
        for unique_cluster in candidate_clusters_df["clusterID"].unique():
            tmp_cluster = candidate_clusters_df[candidate_clusters_df["clusterID"]==unique_cluster]
            tmp_genes = []
            for tmp_gene in tmp_cluster["geneID"].tolist():
                if tmp_gene in candidate_pangenes.keys():
                    tmp_genes.append(candidate_pangenes[tmp_gene])
                else:
                    tmp_genes.append(tmp_gene)
            cluster_genes_dict[unique_cluster] = ",".join(tmp_genes)
        candidate_clusters_df["candidate_genes"] = candidate_clusters_df["clusterID"]
        candidate_clusters_df = candidate_clusters_df.replace({"candidate_genes": cluster_genes_dict})
            
        # Write the cluster ID df into a temporary session file
        tmp_write(session_id, "clusters", candidate_clusters_df)

        # tmp_write(session_id, "candidate_genes", candidate_genes)
        # Try to remove/restart previous temp session files for selected clusters
        try:
            os.remove("tmp/{}_{}.pkl".format(session_id, "selected_cluster"))
        except:
            None

        # Try to write the temp df containing the enriched GO terms for relevant clusters
        try:
            # Read the tmp cluster list df used to parse clusters and get enriched GO terms
            candidate_clusters_df = tmp_read(session_id, "clusters")
            # If the tmp df is empty then return a text specidying it
            if len(candidate_clusters_df)==0:
                return(html.P("No clusters have been identified with the given set of protein IDs."))
            #candidate_clusters_list = candidate_clusters_df["clusterID"].unique()
            candidate_clusters_df = candidate_clusters_df.drop_duplicates("clusterID")
            candidate_clusters_list = candidate_clusters_df["clusterID"].tolist()
            candidate_genes_list = candidate_clusters_df["candidate_genes"].tolist()
            selected_enrichment = get_selected_enrichments(candidate_clusters_list, candidate_genes_list, genotype_list)
            content_enrichment = tab_content_enrichment(selected_enrichment, genotype_list)
            tmp_write(session_id, "enrichment", content_enrichment[-1])
            tmp_write(session_id, "selected_cluster", candidate_clusters_list[0])
            tmp_write(session_id, "sim_clust_list", get_sim_clust_list(candidate_clusters_list[0]))
        except:
            return(html.P("No protein IDs have been submitted for analysis."))


    # Check if the "Submit" button on the first tab was pressed to start process
    if context == "btn_input_GOs.n_clicks":
        # Split the list of gene candidates coming from different rows
        candidate_GOs = GO_list.split("\n")
        # Try to write the temp df containing the enriched GO terms for relevant clusters
        try:
            if len(candidate_GOs)==0:
                return(html.P("No clusters have been identified with the given set of protein IDs."))
            selected_enrichment, relevant_clusters = get_enrichments_with_GOs(candidate_GOs, genotype_list)
            if len(selected_enrichment)==0:
                return(html.P("No clusters have been identified with the given set of GO terms."))
            candidate_clusters_df = cluster_df[cluster_df["clusterID"].isin(relevant_clusters)]
            candidate_clusters_df["candidate_genes"] = [""] * len(candidate_clusters_df)
            tmp_write(session_id, "clusters", candidate_clusters_df)
            
            candidate_clusters_list = candidate_clusters_df["clusterID"].unique().tolist()
            tmp_write(session_id, "selected_cluster", candidate_clusters_list[0])
            tmp_write(session_id, "sim_clust_list", get_sim_clust_list(candidate_clusters_list[0]))
            content_enrichment = tab_content_enrichment(selected_enrichment, genotype_list)
            tmp_write(session_id, "enrichment", content_enrichment)
        except:
            return(html.P("Something went wrong with the submitted GO terms."))
        
    # Update the view based on what tab is selected
    if at == "tab-enrichment":
        content_enrichment = tmp_read(session_id, "enrichment")
        return(content_enrichment)
    elif at == "tab-network":
        content_network = tab_content_network(session_id=session_id,
                                          n_clicks_ids=0,
                                          )
    
        return content_network
    elif at == "tab-annotation":
        return(tab_content_annotation(session_id, annots_df))
    elif at == "tab-enrichment-sim":
        content_enrichment_sim = tmp_read(session_id, "enrichment_similar")
        return(content_enrichment_sim)

@app.callback(Output("state_enrichment_table", "data"),
              Input("enrichment_filter_dropdown", "value"),
              prevent_initial_call=True)
def update_enrichment_table(ls):
    content_enrichment = tmp_read(session_id, "enrichment")
    content_enrichment = tmp_read(session_id, "enrichment")
    dt = content_enrichment[-1].data
    filtered_dt = []
    if not ls:
        return(dt)
         
    for row in dt:
        if row["type"] in ls:
            filtered_dt.append(row)
    return(filtered_dt)

if __name__ == "__main__":
    app.run(host='127.0.0.1', port=8000, debug=True)