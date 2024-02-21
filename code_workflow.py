import os
import sys
import glob
import pickle
import pandas as pd
import numpy as np
from itertools import combinations
import seaborn as sns
import networkx as nx
from multiprocessing import  Pool

os.makedirs("output/", exist_ok=True)

### 1. Extract all STRING predicted physical interactions #####
def main_make_predict_interactomes():
    print("Extracting all STRING predicted physical interactions")
    os.makedirs("output/networks", exist_ok=True)
    os.makedirs("output/networks/genomes/", exist_ok=True)

    for fl in glob.glob("input/genomes/*"):
        genotype = fl.split("/")[-1].split(".")[0]
        ppi = pd.read_csv(fl, sep=" ")

        ppi = ppi[["protein1", "protein2", "experiments_transferred"]]
        ppi = ppi[ppi["experiments_transferred"]!=0]
        ppi["protein1"] = ppi["protein1"].str.split(".").str[-1]
        ppi["protein2"] = ppi["protein2"].str.split(".").str[-1]
        max_w, min_w = ppi["experiments_transferred"].max(), ppi["experiments_transferred"].min()
        ppi["experiments_transferred"] = ppi["experiments_transferred"].apply(lambda x: (x - min_w) / (max_w - min_w))
        ppi.columns = ["source", "target", "weight"]
        
        # Order all gene IDs alphabetically for easier removal of duplicates
        ppi["ordered"] = ppi.apply(lambda x: x["source"]+"-"+x["target"] if x["source"]<x["target"] else x["target"]+"-"+x["source"], axis=1)
        ppi = ppi.drop_duplicates("ordered")
        ppi["source"] = ppi["ordered"].str.split("-").str[0]
        ppi["target"] = ppi["ordered"].str.split("-").str[1]
        ppi = ppi.drop("ordered", axis=1)
        ppi = ppi[ppi["source"]!=ppi["target"]]
        ppi.to_csv("output/networks/genomes/{}.ppi.csv".format(genotype), index=None)

### 2. Convert the gene IDs to pan-gene IDs #####
def main_make_pan_gene_dict():
    # Step 1: Load all the PPI network data into a single dataframe
    pan = pd.DataFrame()
    for fl in glob.glob("output/networks/genomes/*"):
        genotype = fl.split("/")[-1].split(".")[0]
        tmp = pd.read_csv(fl, sep=",")
        tmp["genotype"] = genotype
        pan = pd.concat([pan, tmp])

    # Gene ID to pan ID conversion dictionary
    with open('input/pangenes/maizegdb.pan.pkl', 'rb') as handle:
        pan_dict = pickle.load(handle)

    # Convert to pan-gene ID if possible, otherwise NaN
    def tryconvert(x):
        global pan_dict
        try:
            return(pan_dict[x])
        except:
            return(np.NaN)

    pan["source"] = pan["source"].apply(lambda x: tryconvert(x))
    pan["target"] = pan["target"].apply(lambda x: tryconvert(x))
    pan = pan.dropna()
    pan["ordered"] = pan.apply(lambda x: x["source"]+"-"+x["target"] if x["source"]<x["target"] else x["target"]+"-"+x["source"], axis=1)

    # Step 2 - Generate the pan- and core-interactomes
    # Keep only one pair of unique PPIs per genotype
    pan = pan.drop_duplicates(["genotype", "ordered"])
    # Count number of unique edges per genotype
    pan = pan.groupby("ordered").count()
    pan = pan.reset_index()[["ordered", "source"]]
    pan = pan[pan["source"]>1]
    pan.columns = ["ordered", "universality"]
    pan["source"] = pan["ordered"].str.split("-").str[0]
    pan["target"] = pan["ordered"].str.split("-").str[1]
    pan = pan[pan["source"]!=pan["target"]]
    pan = pan.drop("ordered", axis=1)
    pan = pan[["source", "target", "universality"]]
    pan.to_csv("output/networks/pan.ppi.csv", index=None)
    core = pan[pan["universality"]==26]
    core.to_csv("output/networks/core.ppi.csv", index=None)

### 3. Geenrate all the ClusterONE networks  #####
def main_make_clusterone():
    print("3. Generate all the ClusterONE networks")
    
    os.makedirs("output/clusters", exist_ok=True)

    pan = pd.read_csv("output/networks/pan.ppi.csv")
    pan = pan[["source", "target"]]
    pan.to_csv("tmp/clusterone.tsv", index=None, header=None, sep="\t")
    os.system('java -jar cluster_one-1.0.jar tmp/clusterone.tsv -d 0.1 -f "edge_list" -F "csv" > output/clusters/pan.cone.csv')

    core = pd.read_csv("output/networks/core.ppi.csv")
    core = core[["source", "target"]]
    core.to_csv("tmp/clusterone.tsv", index=None, header=None, sep="\t")
    os.system('java -jar cluster_one-1.0.jar tmp/clusterone.tsv -d 0.2 -f "edge_list" -F "csv" > output/clusters/core.cone.csv')

    os.makedirs("output/clusters/STRING", exist_ok=True)
    for fl in glob.glob("output/networks/genomes/*"):
        
        genotype = fl.split("/")[-1].split(".")[0]
        tmp = pd.read_csv(fl, sep=",")
        tmp.to_csv("tmp/clusterone.tsv", index=None, header=None, sep="\t")
        #print(tmp.head())
        tmp = tmp[["source", "target"]]
        os.system('java -jar cluster_one-1.0.jar tmp/clusterone.tsv -d 0.2 -f "edge_list" -F "csv" > output/clusters/genomes/{}.cone.csv'.format(genotype))
    

    network_dict = {"pan": "output/clusters/pan.cone.csv",
                    "core": "output/clusters/core.cone.csv"}

    for key, val in network_dict.items():
        network = pd.read_csv("output/networks/{}.ppi.csv".format(key))
        
        cls = pd.read_csv(val)
        cls = cls[cls["P-value"]<0.1]
        cls = cls[cls["Size"]<50]
        member_list = []
        for members in cls["Members"].tolist():
            member_list += members.split(" ")
        network = network[(network["source"].isin(member_list)) & (network["target"].isin(member_list))]
        
        network.columns = ["source", "target", "universality"]
        network.to_csv("output/networks_merged/{}.merged.cone.csv".format(key), 
                                    sep=",", index=False)

    for fl in glob.glob("output/clusters/genomes/*"):
        genotype = fl.split("/")[-1].split(".")[0]
        network = pd.read_csv("output/networks/genomes/{}.ppi.csv".format(genotype))
        print(fl)
        cls = pd.read_csv(fl)
        cls = cls[cls["Size"]<50]
        cls = cls[cls["P-value"]<0.1]
        member_list = []
        for members in cls["Members"].tolist():
            member_list += members.split(" ")
        network = network[(network["source"].isin(member_list)) & (network["target"].isin(member_list))]
        network = network[["source", "target"]]
        network.to_csv("output/networks_merged/genomes/{}.merged.cone.csv".format(genotype), 
                                    sep=",", index=False)

### 4. Generate all cluster GO ernichemnt tables  #####
def main_make_cluster_enrichments():
    print("4. Generate all cluster GO ernichemnt tables")

    # First step is to convert all GO annotation files to 
    os.makedirs("output/GO", exist_ok=True)
    pan_go = {}

    # Gene ID to pan ID conversion dictionary
    with open('input/pangenes/maizegdb.pan.pkl', 'rb') as handle:
        pan_dict = pickle.load(handle)

    os.makedirs("output/GO/STRING", exist_ok=True)
    for fl in glob.glob("output/networks/genomes/*"):
        genotype = fl.split("/")[-1].split(".")[0]
        # GO terms obtained from: https://download.maizegdb.org/GeneFunction_and_Expression/Pannzer_GO_Terms/
        df = pd.read_csv("input/GO_PANZER/{}_GO.out".format(genotype), sep="\t", index_col=None)
        df["qpid"] = df["qpid"].str.split("_").str[0]
        df["goid"] = df["goid"].apply(lambda x: "GO:{}{}".format("0"*(7-len(str(x))), x))
        df = df[["qpid", "goid"]]
        df = df.drop_duplicates().dropna()
        df.columns = ["identifier", "label"]
        # All go terms are space separated in the label column
        df = df.groupby("identifier").agg({'label': lambda x: ' '.join(x)}).reset_index()
        df.to_csv("output/GO/genomes/{}.id2gos".format(genotype), sep="\t", index=None)

        for ix, row in df.iterrows():
            try:
                tmp_id = pan_dict[row['identifier']]
                tmp_go = row['label'].split(" ")
                if tmp_id in pan_go.keys():
                    pan_go[tmp_id] = list(set(pan_go[tmp_id] + tmp_go))
                else:
                    pan_go[tmp_id] = list(set(tmp_go))
            except:
                None

    # Combine GO IDs from all available genotypes into a single file
    pan_go_df = {"identifier": [], "label": []}
    for key, vals in pan_go.items():
        pan_go_df["identifier"].append(key)
        pan_go_df["label"].append(" ".join(vals))
    pd.DataFrame(pan_go_df).to_csv("output/GO/pan.id2gos", sep="\t", index=None)



    # Step 2: Conduct the GO enrichement analyses for each cluster
    os.makedirs("output/enrichment", exist_ok=True)
    os.makedirs("output/enrichment/STRING", exist_ok=True)


    from goatools.anno.idtogos_reader import IdToGosReader
    from goatools.obo_parser import GODag
    from goatools.go_enrichment import GOEnrichmentStudy

    import multiprocessing

    go = GODag("input/obo/go-basic.obo", optional_attrs={'relationship'})

    def calculate_enrichment(item):    
        global g
        
        name, genes = item
        results = g.run_study(genes)
        return_list = []
        for result in results:
            # result.enrichment is either e or p for enriched or purified
            if result.p_bonferroni < 0.05 and result.enrichment=="e":
                # https://github.com/tanghaibao/goatools/blob/main/notebooks/Enrichment_analyses_human_phenotype_ontology.ipynb
                fc = (result.ratio_in_study[0]/result.ratio_in_study[1]) / (result.ratio_in_pop[0]/result.ratio_in_pop[1])
                if fc < 2:
                    continue
                return_list.append([name, result.goterm.namespace, len(genes), result.GO, result.p_bonferroni, fc, result.goterm.name])
        return(return_list)

    gos = IdToGosReader("output/GO/pan.id2gos")
    pop = list(gos.id2gos.keys())[1:]
    assoc = {}
    for key, val in dict(gos.associations).items():
        assoc[key] = set(val.split(" "))
    g = GOEnrichmentStudy(pop, assoc, go, propagate_counts=False,
                            alpha=0.05, methods=['bonferroni'])



    for genotype in ["core", "pan"]:
        cone = pd.read_csv("output/clusters/{}.cone.csv".format(genotype))
        cone = cone[cone["P-value"]<0.1]
        cone = cone[cone["Size"]<50]
        clusters = {}
        for ix, row in cone.iterrows():
            clusters["{}_{}".format(genotype, row["Cluster"])] = row["Members"].split(" ")

        # Create a multiprocessing pool with the specified number of processes
        num_processes = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=num_processes)
        results = []
        # Parallel enrichment analysis
        for result in pool.imap(calculate_enrichment, clusters.items()):
            results += result
        # Close the multiprocessing pool
        pool.close()
        pool.join()

        results = [result for result in results if result!=None]
        results = pd.DataFrame(results)
        results.columns = ["cluster", "type", "size", "term", "p-val", "FC", "desc"]
        results.to_csv("output/enrichment/{}.enrichment.tsv".format(genotype), sep="\t", index=None)

    os.makedirs("output/GO/STRING", exist_ok=True)
    for fl in glob.glob("output/networks/genomes/*"):
        genotype = fl.split("/")[-1].split(".")[0]
        print(genotype)

        # Load the GO annotation for each genotype separately
        gos = IdToGosReader("output/GO/genomes/{}.id2gos".format(genotype))
        pop = list(gos.id2gos.keys())[1:]
        assoc = {}
        for key, val in dict(gos.associations).items():
            assoc[key] = set(val.split(" "))
        g = GOEnrichmentStudy(pop, assoc, go, propagate_counts=False,
                                alpha=0.05, methods=['bonferroni'])
        
        cone = pd.read_csv("output/clusters/genomes/{}.cone.csv".format(genotype))
        cone = cone[cone["P-value"]<0.1]
        cone = cone[cone["Size"]<50]
        clusters = {}
        for ix, row in cone.iterrows():
            clusters["{}_{}".format(genotype, row["Cluster"])] = row["Members"].split(" ")

        # Create a multiprocessing pool with the specified number of processes
        num_processes = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=num_processes)
        results = []
        # Parallel enrichment analysis
        for result in pool.imap(calculate_enrichment, clusters.items()):
            results += result
        # Close the multiprocessing pool
        pool.close()
        pool.join()

        results = [result for result in results if result!=None]
        results = pd.DataFrame(results)
        results.columns = ["cluster", "type", "size", "term", "p-val", "FC", "desc"]
        results.to_csv("output/enrichment/genomes/{}.enrichment.tsv".format(genotype), sep="\t", index=None)

### 5. Add coexpression results to clusters  #####
def main_make_coexpression_clusters():

    import numpy as np
    from scipy.special import betainc
    from scipy.stats import beta

    print("5. Add coexpression results to clusters")
    
    def corrcoef(matrix):
        r = np.corrcoef(matrix)    
        rf = r[np.triu_indices(r.shape[0], 1)]
        df = matrix.shape[1] - 2
        ts = rf * rf * (df / (1 - rf * rf))
        pf = betainc(0.5 * df, 0.5, df / (df + ts))
        p = np.zeros(shape=r.shape)
        p[np.triu_indices(p.shape[0], 1)] = pf
        p[np.tril_indices(p.shape[0], -1)] = p.T[np.tril_indices(p.shape[0], -1)]

        # Set diagonal and upper triangle to zero
        r[np.diag_indices(r.shape[0])] = np.zeros(r.shape[0])
        r *= np.tri(*r.shape)
        p[np.diag_indices(p.shape[0])] = np.zeros(p.shape[0]) #np.ones(p.shape[0])
        p *= np.tri(*p.shape)
        return r, p

    def get_cluster_with_coexpression(cluster_genes, edges_df, exp):
        exp_set = set(exp.index.tolist())
        intersecting = list(exp_set.intersection(set(cluster_genes)))
        intersecting.sort()
        pr_df = []
        df = exp.loc[intersecting]
        n=df.shape[1]
        corr=df.T.corr()
        dist = beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
        pval = 2 * dist.cdf(-abs(corr))
        pval[np.diag_indices(pval.shape[0])] = np.ones(pval.shape[0])
        corr *= np.tri(*corr.shape)

        # Reversed the orders because I removed the upper triangle, revesing the alphabetical order
        for p2, p1 in zip(*np.where((corr>0.9) & (pval < 0.01))):
            #print(corr.iloc[1,1])
            #print(pval[1][1])
            #print(corr, pval)
            pr_df.append([intersecting[p1], intersecting[p2], corr.iloc[p2,p1], 1]) # pval[p1][p2] to add pval column

        pr_df = pd.DataFrame(pr_df, columns=["source", "target", "corr", "weight"])
        # Make sure the edge pairs are also alphabetically sorted
        tmp = edges_df[(edges_df["source"].isin(cluster_genes)) & (edges_df["target"].isin(cluster_genes))]
        tmp.insert(len(tmp.columns), 'weight', -1)
        
        pr_df = pd.DataFrame(pr_df)
        pr_df = pd.concat([pr_df, tmp])

        # When the weights is -1: only PPI, 0: both, 1: only coexpressed
        pr_df = pr_df.groupby(["source", "target"]).sum()
        pr_df = pr_df.reset_index()
        pr_df["tmp"] = ""
        pr_df.loc[pr_df["weight"]==-1, 'tmp'] = 'PPI'
        pr_df.loc[pr_df["weight"]==0, 'tmp'] = 'Both'
        pr_df.loc[pr_df["weight"]==1, 'tmp'] = 'Coexp'
        pr_df["weight"] = pr_df["tmp"]
        pr_df = pr_df.drop("tmp", axis=1)
        # Change column type from object to float for rounding        
        pr_df["corr"] = pr_df["corr"].astype(float)
        pr_df["corr"] = pr_df["corr"].round(5)
        return(pr_df)
        
    os.makedirs("output/coexpression", exist_ok=True)
    os.makedirs("output/coexpression/STRING", exist_ok=True)

    for fl in glob.glob("output/networks/genomes/*"):
        genotype = fl.split("/")[-1].split(".")[0]
        #if genotype != "B73":
        #    continue
        final_df = pd.DataFrame()
        edges_df = pd.read_csv(fl)
        edges_df = edges_df.drop("weight", axis=1)
        cone = pd.read_csv("output/clusters/genomes/{}.cone.csv".format(genotype))
        cone = cone[cone["P-value"]<0.1]
        cone = cone[cone["Size"]<50]
        clusters = {}
        for ix, row in cone.iterrows():
            clusters["{}_{}".format(genotype, row["Cluster"])] = row["Members"].split(" ")
        exp = pd.read_csv("input/expression/{}.FPKM.csv.gz".format(genotype), index_col=0)
        for cluster in clusters.keys():        
            pr_df = get_cluster_with_coexpression(clusters[cluster], edges_df, exp)
            pr_df["cluster"] = cluster
            final_df = pd.concat([final_df, pr_df])
        final_df = final_df[final_df["source"]!=final_df["target"]]
        final_df.to_csv("output/coexpression/genomes/{}.coexp.csv".format(genotype), index=None)
        
        print(genotype, final_df.shape)

    with open('input/pangenes/maizegdb.pan.pkl', 'rb') as handle:
        pan_dict = pickle.load(handle)
    pan_df = pd.DataFrame(pd.Series(pan_dict)).reset_index()
    pan_df.columns = ["geneID", "panID"]

    def get_pan_cluster_with_coexpression(cluster_genes, edges_df, exp, pan_df):
        original_cluster_gene = list(cluster_genes) # copy_list
        original_cluster_gene.sort()
        
        cluster_genes = pan_df[pan_df["panID"].isin(cluster_genes)]["geneID"].unique().tolist()
        
        exp_set = set(exp.index.tolist())
        intersecting = list(exp_set.intersection(set(cluster_genes)))
        intersecting.sort()
        pr_df = []

        df = exp.loc[intersecting]
        n=df.shape[1]
        corr=df.T.corr()
        dist = beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
        pval=2*dist.cdf(-abs(corr))
        pval[np.diag_indices(pval.shape[0])] = np.ones(pval.shape[0])
        corr *= np.tri(*corr.shape)

        # Reversed the orders because I removed the upper triangle, revesing the alphabetical order
        for p2, p1 in zip(*np.where((corr>0.9) & (pval < 0.01))):
            pr_df.append([intersecting[p1], intersecting[p2], intersecting[p1], intersecting[p2], 
                        corr.iloc[p2,p1], 1])
        
        pr_df = pd.DataFrame(pr_df, columns=["source", "target", "s_ori", "t_ori", "corr", "weight"])
        #return pr_df

        source_col = []
        target_col = []
        for ix, row in pr_df.iterrows():
            tmp = [pan_dict[row["source"]], pan_dict[row["target"]]]
            tmp.sort()
            source_col.append(tmp[0])
            target_col.append(tmp[1])
        pr_df["source"] = source_col
        pr_df["target"] = target_col
        
        # In case there are multiple mapped correlated pan-genes
        pr_df = pr_df.drop_duplicates(["source", "target"])
        
        # Make sure the edge pairs are also alphabetically sorted
        tmp = edges_df[(edges_df["source"].isin(original_cluster_gene)) & (edges_df["target"].isin(original_cluster_gene))]
        tmp.insert(len(tmp.columns), 'weight', -1)
        
        pr_df = pd.DataFrame(pr_df)
        pr_df = pd.concat([pr_df, tmp])

        # When the weights is -1: only PPI, 0: both, 1: only coexpressed
        pr_df = pr_df.groupby(["source", "target"]).sum()
        pr_df = pr_df.reset_index()
        pr_df["tmp"] = ""
        pr_df.loc[pr_df["weight"]==-1, 'tmp'] = 'PPI'
        pr_df.loc[pr_df["weight"]==0, 'tmp'] = 'Both'
        pr_df.loc[pr_df["weight"]==1, 'tmp'] = 'Coexp'
        pr_df["weight"] = pr_df["tmp"]
        pr_df = pr_df.drop("tmp", axis=1)
        # Change column type from object to float for rounding        
        pr_df["corr"] = pr_df["corr"].astype(float)
        pr_df["corr"] = pr_df["corr"].round(5)
        return(pr_df)

    for network in ["core"]: #["pan"]:
        final_df = pd.DataFrame()
        for fl in glob.glob("output/networks/genomes/*"):
            genotype = fl.split("/")[-1].split(".")[0]
            edges_df = pd.read_csv(fl)
            edges_df = edges_df.drop("weight", axis=1)

            edges_df = pd.read_csv("output/networks/{}.ppi.csv".format(network))
            edges_df = edges_df.drop("universality", axis=1)
            
            cone = pd.read_csv("output/clusters/{}.cone.csv".format(network))
            cone = cone[cone["P-value"]<0.1]
            cone = cone[cone["Size"]<50]
            clusters = {}
            for ix, row in cone.iterrows():
                clusters["{}_{}".format(network, row["Cluster"])] = row["Members"].split(" ")
            
            exp = pd.read_csv("input/expression/{}.FPKM.csv.gz".format(genotype), index_col=0)
            for cluster in clusters.keys():
                pr_df = get_pan_cluster_with_coexpression(clusters[cluster], edges_df, exp, pan_df)
                pr_df["cluster"] = cluster
                final_df = pd.concat([final_df, pr_df])
        print(final_df.shape)
        final_df = final_df[final_df["source"]!=final_df["target"]]
        print(final_df.shape)
        # Sort by correlation so when dropping duplicates keep highest correlation
        final_df = final_df.sort_values("corr", ascending=False)
        final_df = final_df.drop_duplicates(["source", "target", "cluster"])
        print(final_df.shape)
        final_df.to_csv("output/coexpression/{}.coexp.csv".format(network), index=None)

### 6. Get top-DIAMON-based gene annotations
def main_make_gene_annotation():
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    os.makedirs("output", exist_ok=True)
    os.makedirs("output/descriptions", exist_ok=True)


    desc_df = pd.read_csv("input/description/araport11_descriptions.tsv", sep="\t")
    
    ref = "input/description/Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa"
    os.system("diamond makedb --in "+ref+" --db tmp/tmp.dmnd")

    for fl in glob.glob("output/coexpression/genomes/*"):

        genotype = fl.split("/")[-1].split(".")[0]
        print(genotype, fl)        
        fasta_name = "input/FASTA/{}.fa".format(genotype)
        # Run blastp while returning only the top hit and report unaligned queries
        os.system("diamond blastp -d tmp/tmp.dmnd -q "+fasta_name+" -o tmp/dmnd.out --max-target-seqs 1 --unal 1 --quiet")
        
        df = pd.read_csv("tmp/dmnd.out", sep="\t", index_col=0, header=None)
        df = pd.DataFrame(df[1])
        df.columns = [genotype]
        
        # For some unkown reason some sequences (about 100) get duplicated rows 
        # Deduplicating returns the correct number of rows though (=number of seqs in fasta)
        df = df[~df.index.duplicated(keep="first")]
        first = False
        
        df = df.reset_index()
        df.columns = ["geneID", "name"]
        df = pd.merge(df, desc_df, on="name", how="left")
        #df = df[["GeneID", "name", ]]
        df = df.drop("gene_model_type", axis=1)
        print(df.head(3))
        df.to_csv("output/descriptions/{}.at11.tsv".format(genotype), sep="\t", index=None)

    with open('input/pangenes/maizegdb.pan.pkl', 'rb') as handle:
        pan_dict = pickle.load(handle)
    pan_df = pd.DataFrame(pd.Series(pan_dict)).reset_index()
    pan_df.columns = ["geneID", "panID"]

    # Generate the pan-gene description file based on longest description across genomes
    with open('input/pangenes/maizegdb.pan.pkl', 'rb') as handle:
        pan_dict = pickle.load(handle)
    pan_df = pd.DataFrame(pd.Series(pan_dict)).reset_index()
    pan_df.columns = ["geneID", "panID"]

    desc_df = pd.DataFrame()

    for fl in glob.glob("output/descriptions/*"):
        genotype = fl.split("/")[-1].split(".")[0]
        if genotype in ["pan", "core"]:
            continue
        desc = pd.read_csv("output/descriptions/{}.at11.tsv".format(genotype), sep="\t")
        desc_df = pd.concat([desc_df, desc])

    pan_desc = pd.merge(pan_df, desc_df, on="geneID", how="left")
    pan_desc["len"] = 0
    pan_desc["len"] += pan_desc["short_description"].apply(lambda x: len(str(x)))
    pan_desc["len"] += pan_desc["Curator_summary"].apply(lambda x: len(str(x)))
    pan_desc["len"] += pan_desc["Computational_description"].apply(lambda x: len(str(x)))
    pan_desc = pan_desc.sort_values("len", ascending=False)
    pan_desc = pan_desc.drop_duplicates("panID")
    pan_desc.to_csv("output/descriptions/pan.at11.tsv", sep="\t", index=None)

#main_make_clusterone()
#main_make_gene_annotation()
main_make_coexpression_clusters()