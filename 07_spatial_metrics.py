#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings("ignore")

import os
import matplotlib.pyplot as plt
import squidpy as sq
import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
import cellcharter as cc
import scanpy as sc

# =========================
# SETTINGS
# =========================

CLUSTER_KEY = "new_cluster" ##Name of your cluster column

SAVE_DIR = "../spatial_clustering/ECM1"

ADATA_PATH = "../DF_w_the_clusters.csv"
os.makedirs(SAVE_DIR, exist_ok=True)

# =========================
# LOAD DATA
# =========================

adata = ad.read_h5ad(ADATA_PATH)
print(adata)
adata.obs['ECM_cluster'].astype('category')
adata = adata[adata.obs['ECM_cluster'].notna()]

adata_sub = adata[adata.obs["ECM_cluster"] == 1.0].copy()
print(adata_sub.obs['ECM_cluster'].unique())
# =========================
# HELPER FUNCTION
# =========================

def save_plot(filename):
    plt.savefig(
        os.path.join(SAVE_DIR, filename),
        dpi=75,
        bbox_inches="tight"
    )
    plt.close()

# =========================
# NEIGHBORHOOD ENRICHMENT
# =========================

sq.gr.nhood_enrichment(
    adata_sub,
    cluster_key=CLUSTER_KEY
)

sq.pl.nhood_enrichment(
    adata_sub,
    cluster_key=CLUSTER_KEY
)

save_plot("nhood_enrichment.svg")

# =========================
# INTERACTION MATRIX
# =========================

sq.gr.interaction_matrix(
    adata_sub,
    cluster_key=CLUSTER_KEY
)

sq.pl.interaction_matrix(
    adata_sub,
    cluster_key=CLUSTER_KEY,
    figsize=(4, 4),
    cmap="bwr"
)

save_plot("interaction_matrix.svg")

# =========================
# Neighborhood Enrichment 
# =========================

adata.obs['new_cluster'] = adata.obs['new_cluster'].astype('category')
cc.gr.nhood_enrichment(
    adata,
    cluster_key='new_cluster',
    pvalues=True,  
    n_perms=1000,
    n_jobs = 12
    
)

cc.pl.nhood_enrichment(
    adata,
    cluster_key='new_cluster',
    annotate=True,
    vmin=-1,
    vmax=1,
    fontsize=5,
)

save_plot("neighborhood_matrix.svg")

# =========================
# Diff Neighborhood MATRIX (1v3)
# =========================
adata.obs['new_cluster'] = adata.obs['new_cluster'].astype('category')

cc.gr.diff_nhood_enrichment(
    adata,
    cluster_key='new_cluster',
    condition_key='ECM_cluster',
    library_key = 'filename',
    pvalues=True,
    n_perms=1000,
    n_jobs = 12

)


cc.pl.diff_nhood_enrichment(
    adata,
    cluster_key='new_cluster',
    condition_key='ECM_cluster',
    annotate=True,
    vmin=-1,
    vmax=1,
    fontsize=5,
)
save_plot("Diff_neighborhood_matrix.svg")

# =========================
# RIPLEY L
# =========================

sq.gr.ripley(
    adata_sub,
    cluster_key=CLUSTER_KEY,
    mode="L"
)

sq.pl.ripley(
    adata_sub,
    cluster_key=CLUSTER_KEY,
    mode="L"
)

save_plot("ripley_L.svg")

# =========================
# CO-OCCURRENCE
# =========================

sq.gr.co_occurrence(
    adata_sub,
    cluster_key=CLUSTER_KEY,
    n_splits=20,
    interval=50,
    n_jobs=21
)

sq.pl.co_occurrence(
    adata_sub,
    cluster_key=CLUSTER_KEY
)

save_plot("co_occurrence.svg")

# =========================
# CENTRALITY
# =========================

sq.gr.centrality_scores(
    adata_sub,
    cluster_key=CLUSTER_KEY,
    n_jobs=21
)

sq.pl.centrality_scores(
    adata_sub,
    cluster_key=CLUSTER_KEY,
    figsize=(20, 5),
    s=500
)

save_plot("centrality.svg")

# =========================
# RIPLEY F/G
# =========================

for mode in ["F", "G"]:

    sq.gr.ripley(
        adata_sub,
        cluster_key=CLUSTER_KEY,
        mode=mode,
        n_observations=len(adata_sub)
    )

    sq.pl.ripley(
        adata_sub,
        cluster_key=CLUSTER_KEY,
        mode=mode
    )

    save_plot(f"ripley_{mode}.svg")

plt.tight_layout()
plt.show()

'''
