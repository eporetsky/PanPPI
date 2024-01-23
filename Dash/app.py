"""
Simple Dash app for deploying on personal computer to access the PanPPI data
For more information see https://github.com/eporetsky/PanPPI
To run the app, move the app.py and functions.py files and the tmp/ folder to
    the folder containing all the analysis results (together with the input/ 
    and output/ folders). When ready, in the terminal run: "python app.py"
Developed by Elly Poretsky on 01.22.24
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
    #"top": 0,
    #"left": 0,
    #"bottom": 0,
    "margin-left": "20px",
    #"width": "300px",
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

cluster_df = get_cluster_df()
annots_df = get_phytozome_annotations()

sidebar = html.Div(
    [
        html.H2("Menu", className="display-4"),
        html.Hr(),
        html.P("Parse the pan-interactome data using either a list of gene of interests or a selected cluster. Press the relevant Submit button to parse the data, overriding the previous results.", className="lead"),
        
        html.Hr(),
        html.P("1. Add genes of interest, one gene per row", className="lead"),
        dcc.Textarea(id="input_genes", 
                  value=input_gene_list),
        dbc.Button("Submit", id="btn_input_genes", color="primary", className="me-1"),

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
              Input("btn_network_ids", "n_clicks"),
              #Input("btn_network_coexpression", "n_clicks"),
              Input("dropdown_clusters", "value"),
              prevent_initial_call=True)
def btn_network_tab(n_clicks_ids, selected_cluster):
    context = [p["prop_id"] for p in dash.callback_context.triggered][0]
    n_clicks_ids = 0 if n_clicks_ids == None else n_clicks_ids
    #n_clicks_coexp = 0 if n_clicks_coexp == None else n_clicks_coexp
    tmp_write(session_id, "selected_cluster", selected_cluster)
    content_network = tab_content_network(session_id=session_id,
                                          n_clicks_ids=n_clicks_ids, 
                                          #n_clicks_coexp=n_clicks_coexp,  
                                          )
    return(content_network) 

@app.callback(Output("content", "children"),
              State("input_genes", "value"),
              State("dropdown_genotypes", "value"),
              Input("tabs", "active_tab"),
              Input("btn_input_genes", "n_clicks"),
              prevent_initial_call=True)
def switch_tab(gene_list, genotype_list, at, button):
    context = [p["prop_id"] for p in dash.callback_context.triggered][0]
    
    # Check if the "Submit" button on the first tab was pressed to start process
    if context == "btn_input_genes.n_clicks":
        # Split the list of gene candidates coming from different rows
        candidate_genes = gene_list.split("\n")
        # Add to the list of gene candidates the list of associcated pan-gene IDs
        candidate_genes += get_candidate_pan_genes(candidate_genes)
        # Get all the IDs of clusters associated with any of the candidate genes
        candidate_clusters_df = cluster_df[cluster_df["geneID"].isin(candidate_genes)]
        # Write the cluster ID df into a temporary session file
        tmp_write(session_id, "clusters", candidate_clusters_df)
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
        candidate_clusters_list = candidate_clusters_df["clusterID"].unique()
        selected_enrichment = get_selected_enrichments(candidate_clusters_list, genotype_list)
        content_enrichment = tab_content_enrichment(selected_enrichment, genotype_list)
        tmp_write(session_id, "enrichment", content_enrichment[-1])
    except:
        return(html.P("No protein IDs have been submitted for analysis."))

    # Update the view based on what tab is selected
    if at == "tab-enrichment":
        return(content_enrichment)
    elif at == "tab-network":
        content_network = tab_content_network(session_id)
        return content_network
    elif at == "tab-annotation":
        return(tab_content_annotation(session_id, annots_df))

@app.callback(Output("state_enrichment_table", "data"),
              Input("enrichment_filter_dropdown", "value"),
              prevent_initial_call=True)
def update_enrichment_table(ls):
    dt = tmp_read(session_id, "enrichment").data
    filtered_dt = []
    if not ls:
        return(dt)
         
    for row in dt:
        if row["type"] in ls:
            filtered_dt.append(row)
    return(filtered_dt)

if __name__ == "__main__":
    app.run_server(port=8888, debug=True)