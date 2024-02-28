import os
import glob
import pickle
import pandas as pd
import networkx as nx

os.makedirs("pickles", exist_ok=True)

# A function that generates tmp pickle file with the current session_id
def write_pickle(file_name, content):
    with open('pickles/{}.pkl'.format(file_name), 'wb') as handle:
        pickle.dump(content, handle)#, protocol=pickle.HIGHEST_PROTOCOL)


sim = pd.read_csv("output/clusters.jaccard.05.edges", sep="\t")
G = nx.from_pandas_edgelist(sim[["source", "target"]])
cc_list = [list(G.subgraph(c).copy().nodes()) for c in nx.connected_components(G)]
write_pickle("sim_connected_components", cc_list)



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
        
        cone = cone[cone["P-value"]<0.1]
        tmp = {"clusterID": [], "geneID": []}
        for ix, row in cone.iterrows():
            clusterID = "{}_{}".format(genotype, row["Cluster"])
            for geneID in row["Members"].split(" "):
                tmp["clusterID"].append(clusterID)
                tmp["geneID"].append(geneID)
        cluster_df = pd.concat([cluster_df, pd.DataFrame(tmp)])
    return(cluster_df)

cluster_df = get_cluster_df()
cluster_df.to_pickle("pickles/cluster_df.pkl")

def get_phytozome_annotations():  
    annot = pd.DataFrame()
    for file_name in glob.glob("output/descriptions/*"):
        tmp = pd.read_csv(file_name, sep="\t")
        if "pan.at11.tsv" in file_name:
            tmp = tmp.drop(["geneID", "len"], axis=1)
            tmp.columns = ["geneID", "name", "short_description", "Curator_summary", "Computational_description"]
        annot = pd.concat([annot, tmp])
    return(annot)

annots_df = get_phytozome_annotations()
annots_df.to_pickle("pickles/annots_df.pkl")