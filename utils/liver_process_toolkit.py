import os
import anndata as ad
from matplotlib import pyplot as plt
import pandas as pd
import re
import json
import requests
from scipy.stats import zscore
import numpy as np
import scanpy as sc

from config.settings import OPENAI_API_KEY
from utils.LLM import query_llm

def get_top_differential_genes(marker_file, n_genes=10,cluster="cluster",foldchange = 'avg_log2FC',gene="gene"):
    df = pd.read_csv(marker_file)
    sorted_df = df.sort_values(by=[cluster, foldchange], ascending=[True, False])
    top_genes_dict = sorted_df.groupby(cluster).apply(lambda x: x.head(n_genes)[gene].tolist()).to_dict()
    return top_genes_dict

def truncate_eval(text):
    start = text.find('{')
    end = text.rfind('}') + 1
    json_str = text[start:end]

    # Print the extracted dictionary
    print(json_str)  

def identify_marker_genes(dotplot_data, dotplot_data_frac, exp_thresh=0.5, frac_thresh=0.5):
    """
    Identify potential marker genes for each cluster based on expression and fraction thresholds.

    Args:
    - dotplot_data: DataFrame of average expression levels.
    - dotplot_data_frac: DataFrame of fraction of cells expressing each gene.
    - exp_thresh: Minimum average expression level to consider a gene as a marker.
    - frac_thresh: Minimum fraction of cells required to express a gene for it to be a marker.

    Returns:
    - marker_genes: Dictionary where keys are clusters and values are lists of marker genes.
    """
    marker_genes = {}

    for cluster in dotplot_data.index:
        marker_genes[cluster] = []
        for gene in dotplot_data.columns:
            if dotplot_data.loc[cluster, gene] >= exp_thresh and dotplot_data_frac.loc[cluster, gene] >= frac_thresh:
                marker_genes[cluster].append(gene)

    return marker_genes

def zscore_normalize_expression(dotplot_data):
    """
    Apply Z-score normalization to gene expression levels.

    Args:
    - dotplot_data: DataFrame of average expression levels.

    Returns:
    - normalized_data: DataFrame of Z-score normalized expression levels.
    """
    normalized_data = dotplot_data.apply(zscore, axis=0)  # Normalize each gene column-wise
    return normalized_data

def min_max_scale_expression(dotplot_data):
    """
    Apply min-max scaling to gene expression levels.

    Args:
    - dotplot_data: DataFrame of average expression levels.

    Returns:
    - scaled_data: DataFrame of min-max scaled expression levels.
    """
    scaled_data = dotplot_data.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0)
    return scaled_data

def compute_logfc(dotplot_data, cluster1, cluster2):
    """
    Compute log-fold change of gene expression between two clusters.

    Args:
    - dotplot_data: DataFrame of average expression levels.
    - cluster1: First cluster for comparison.
    - cluster2: Second cluster for comparison.

    Returns:
    - logfc: Series of log-fold change values for all genes.
    """
    logfc = np.log2(dotplot_data.loc[cluster1] + 1) - np.log2(dotplot_data.loc[cluster2] + 1)
    return logfc

def identify_distinguishing_markers(dotplot_data, dotplot_data_frac, cluster1, cluster2, exp_thresh=0.5, frac_thresh=0.5, logfc_thresh=1.0):
    """
    Identify marker genes that can distinguish between two clusters.

    Args:
    - dotplot_data: DataFrame of average expression levels.
    - dotplot_data_frac: DataFrame of fraction of cells expressing each gene.
    - cluster1: First cluster for comparison.
    - cluster2: Second cluster for comparison.
    - exp_thresh: Minimum average expression level to consider a gene as a marker.
    - frac_thresh: Minimum fraction of cells required to express a gene for it to be a marker.
    - logfc_thresh: Minimum log-fold change required for a gene to be considered distinguishing.

    Returns:
    - distinguishing_markers: List of genes that distinguish between the two clusters.
    """
    logfc = compute_logfc(dotplot_data, cluster1, cluster2)
    
    distinguishing_markers = []
    
    for gene in dotplot_data.columns:
        if (dotplot_data.loc[cluster1, gene] >= exp_thresh or dotplot_data.loc[cluster2, gene] >= exp_thresh) and \
           (dotplot_data_frac.loc[cluster1, gene] >= frac_thresh or dotplot_data_frac.loc[cluster2, gene] >= frac_thresh) and \
           abs(logfc[gene]) >= logfc_thresh:
            distinguishing_markers.append(gene)
    
    return distinguishing_markers

def find_similar_cluster_pairs(dotplot_data, dotplot_data_frac, exp_diff_thresh=0.5, frac_diff_thresh=0.5, max_diff_genes=2, logfc_thresh=1.0):
    """
    Identify similar cluster pairs where most genes have similar expression profiles,
    but they differ in a small number of genes.

    Args:
    - dotplot_data: DataFrame of average expression levels.
    - dotplot_data_frac: DataFrame of fraction of cells expressing each gene.
    - exp_diff_thresh: Threshold for considering expression levels as similar.
    - frac_diff_thresh: Threshold for considering fraction of expression as similar.
    - max_diff_genes: Maximum number of genes where clusters can differ.
    - logfc_thresh: Minimum log-fold change required for a gene to be considered distinguishing.

    Returns:
    - similar_pairs: List of tuples containing similar cluster pairs and their distinguishing genes.
    """
    clusters = dotplot_data.index
    similar_pairs = []
    
    for i, cluster1 in enumerate(clusters):
        for cluster2 in clusters[i+1:]:
            diff_genes = identify_distinguishing_markers(dotplot_data, dotplot_data_frac, cluster1, cluster2,
                                                         exp_thresh=exp_diff_thresh,
                                                         frac_thresh=frac_diff_thresh,
                                                         logfc_thresh=logfc_thresh)
            # If the number of differing genes is below or equal to the threshold, consider them similar
            if len(diff_genes) <= max_diff_genes:
                similar_pairs.append((cluster1, cluster2, diff_genes))
    
    return similar_pairs

import anndata as ad
import ast

def solve_all_clusters(missing_list,marker_file,model_provider, model_name, API_key,info=None,input_dir='data/liver/input'):
    print("in auto fill in ")
    unsolved_list = missing_list
    marker_data = pd.read_csv(os.path.join(input_dir, marker_file))
    top_genes = {}
    for cluster, group in marker_data.groupby('cluster'):
        # Store the gene names in a list for each cluster
        top_genes[cluster] = list(group['gene'])
    subset_top_genes = {key: top_genes[key] for key in unsolved_list if key in top_genes}
    system_role = "You are expert in scRNA sequencing cell type annotation."
    content = f'''
            this is background information: {info}
            look at this dict: {subset_top_genes}. This is cluster number and the corresponding top differential genes of each cluster. Please provide cell type annotation for each cluster. 
            Output in text dict format just like the input dict. Keys are number of cluster, and Values are strings of cell type names. Output should be text dict, no other word should exist.
    '''
    LLM_output = query_llm(content=content,system_role=system_role,model_provider = model_provider,model_name=model_name,API_key = API_key)
    reply = LLM_output
    sanitized_str = reply.replace("```", "")
    try:
        sanitized_dict =  ast.literal_eval(sanitized_str)
    except (ValueError, SyntaxError):
        pass
    print(sanitized_dict)

    return sanitized_dict