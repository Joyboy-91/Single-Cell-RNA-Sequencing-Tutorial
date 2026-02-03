import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import bbknn

# ============================================================
# 1. SETTINGS & DIRECTORY SETUP
# ============================================================

# Configure Scanpy aesthetics for publication-quality figures
sc.settings.verbosity = 3 # Show detailed logs during execution
sc.set_figure_params(dpi=300, fontsize=8, facecolor="white", frameon=False, vector_friendly=True)

# Define input/output paths using pathlib for cross-platform compatibility
BASE_DIR = Path.cwd()
WRITE_DIR = BASE_DIR / "write"
RESULT_DIR = Path.home() / "results"

# Automatically create the necessary folders if they don't exist
for d in [WRITE_DIR, RESULT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Create specific subdirectories for Rat and Human results
for fig in ["1_preprocessing_and_clustering_Rat", "1_preprocessing_and_clustering_Human"]:
    (RESULT_DIR / fig).mkdir(exist_ok=True)

print(f"Results will be saved to: {RESULT_DIR}")

# ============================================================
# 2. HELPER FUNCTIONS (The Analytical Core)
# ============================================================

def load_and_standardize(path, sample_id, condition, species):
    """
    Loads raw 10x matrix data and standardizes gene names (to UPPERCASE).
    Injects essential metadata (condition, sample_id, species) into the AnnData object.
    """
    try:
        adata = sc.read_10x_mtx(BASE_DIR / path, var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        adata.var_names = adata.var_names.str.upper() # Standardize for cross-species comparison
        adata.obs["sample_id"] = sample_id
        adata.obs["condition"] = condition
        adata.obs["species"] = species
        return adata
    except Exception as e:
        print(f"ERROR: Could not load {path} -> {e}")
        return None

def process_basic(adata, batch_key):
    """
    Performs standard Single-Cell RNA-seq Quality Control (QC), 
    Normalization, Batch Correction, and Clustering (t-SNE + UMAP).
    """
    # 1. QUALITY CONTROL (QC)
    # Identify mitochondrial genes (high mt-RNA indicates dying cells)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    
    print(f"Number of cells before filtering: {adata.n_obs}")
    
    # Filter cells based on count and mitochondrial thresholds.
    # CRITICAL: n_genes > 400 is deliberately chosen (instead of a stricter threshold)
    # to preserve rare cell types (e.g., specific fibroblast subtypes) that naturally have lower gene counts.
    adata = adata[
        (adata.obs.total_counts > 800) & 
        (adata.obs.total_counts < 20000) & 
        (adata.obs.pct_counts_mt < 15) & 
        (adata.obs.n_genes_by_counts > 400) 
    ].copy()
    print(f"Number of cells after filtering: {adata.n_obs}")
    
    # 2. NORMALIZATION & LOG-TRANSFORMATION
    sc.pp.normalize_total(adata, target_sum=1e4) # Normalize counts per cell to 10k
    sc.pp.log1p(adata) # Log(x+1) transformation for statistical stability
    adata.raw = adata  # Save raw data for accurate differential expression later
    
    # 3. FEATURE SELECTION (Highly Variable Genes - HVG)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
    
    # 4. REGRESS OUT TECHNICAL NOISE & SCALE
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata) # Scale to mean=0, std=1 for PCA
    
    # 5. DIMENSIONALITY REDUCTION (PCA)
    sc.pp.pca(adata, mask_var="highly_variable")
    
    # 6. BATCH EFFECT CORRECTION (BBKNN)
    # Merges control and irradiated datasets by correcting the neighborhood graph
    # directly without altering the PCA space.
    bbknn.bbknn(adata, batch_key=batch_key)
    
    # 7. EMBEDDINGS & CLUSTERING
    print("Calculating t-SNE...")
    sc.tl.tsne(adata, n_pcs=30) # t-SNE used for global atlas visualization
    
    print("Calculating UMAP...")
    sc.tl.umap(adata) # UMAP is strictly required for PAGA/Trajectory analysis in Phase 2
    
    # Cluster cells using Leiden algorithm (Resolution 2.0 to capture rare sub-clusters)
    sc.tl.leiden(adata, resolution=2) 
    
    return adata

def auto_annotate_cell_types(adata, marker_dict):
    """
    Automatically annotates Leiden clusters based on the maximum mean expression 
    of known marker gene signatures.
    """
    # Clear old scores if running multiple times
    old_scores = [c for c in adata.obs.columns if c.startswith("score_")]
    if old_scores: adata.obs.drop(columns=old_scores, inplace=True)
    
    print("\n--- Marker Gene Check ---")
    # Score each cell for each cell type's gene signature
    for cell_type, genes in marker_dict.items():
        valid_genes = [g for g in genes if g in adata.var_names]
        if len(valid_genes) > 0:
            sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f"score_{cell_type}")
        else:
            print(f"WARNING: NO marker genes found for {cell_type}!")
    
    score_cols = [col for col in adata.obs.columns if col.startswith("score_")]
    if not score_cols: return adata

    # Assign cell type to each cluster based on the highest average signature score
    cluster_scores = adata.obs.groupby("leiden")[score_cols].mean()
    cluster_mapping = {}
    for cluster in cluster_scores.index:
        best_score_col = cluster_scores.loc[cluster].idxmax()
        cell_name = best_score_col.replace("score_", "")
        cluster_mapping[cluster] = cell_name
    
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_mapping)
    return adata

def plot_split_tsne(adata, keys, key_col, color_col, save_path, figsize_per_plot=(4, 4)):
    """
    Draws condition-specific t-SNE plots side-by-side.
    Ensures all plots share the exact same coordinate limits for fair visual comparison.
    """
    n_plots = len(keys)
    fig, axes = plt.subplots(1, n_plots, figsize=(figsize_per_plot[0] * n_plots, figsize_per_plot[1]))
    
    if n_plots == 1: axes = [axes] # Convert to list if there is only one plot
    
    # Calculate global t-SNE limits with a 5% padding so points don't stick to the border
    x_min, x_max = adata.obsm['X_tsne'][:, 0].min(), adata.obsm['X_tsne'][:, 0].max()
    y_min, y_max = adata.obsm['X_tsne'][:, 1].min(), adata.obsm['X_tsne'][:, 1].max()
    pad_x = (x_max - x_min) * 0.05
    pad_y = (y_max - y_min) * 0.05

    for i, key in enumerate(keys):
        ax = axes[i]
        subset = adata[adata.obs[key_col] == key] # Extract specific condition (e.g., '7d')
        
        # Plot individual condition
        sc.pl.tsne(
            subset, 
            color=color_col, 
            ax=ax, 
            show=False, 
            title=key, 
            frameon=False,
            legend_loc="none" if i < n_plots - 1 else "right margin", 
            s=50 # Point size
        )
        
        # Lock axis limits for all plots
        ax.set_xlim(x_min - pad_x, x_max + pad_x)
        ax.set_ylim(y_min - pad_y, y_max + pad_y)
        
        # Clean axis labels 
        ax.set_xlabel("tSNE_1")
        if i == 0: ax.set_ylabel("tSNE_2")
        else: ax.set_ylabel("")
            
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved split t-SNE to: {save_path}")

# ============================================================
# MARKER LISTS
# Defined based on canonical literature markers for each species
# ============================================================
MARKERS_RAT = {
    "Keratinocytes (KC)": ["KRT1", "KRT14", "KRT5"],
    "Fibroblasts (FB)": ["DCN", "APOD", "CRABP1"],
    "Endothelial (EC)": ["TM4SF1", "CDH5", "CD93"],
    "Pericytes (PC)": ["RGS5", "DES", "PDGFRB"],
    "Schwann (SC)": ["GATM", "MPZ", "PLP1"],          
    "Smooth_Muscle (SMC)": ["MYH11", "MYL9", "TAGLN"],
    "Myoblasts (MB)": ["MYF5", "JSRP1", "CDH15"],
    "Neural (NC)": ["CMTM5", "BCHE", "AJAP1"],
    "Macrophages (Mø)": ["C1QC", "CD68", "PF4"],
    "Neutrophils (NEUT)": ["S100A8", "S100A9", "LYZ2"],
    "T_Cells (TC)": ["CD3D", "CD3E", "ICOS"],
    "B_Cells (BC)": ["TYROBP", "IFI30", "IGHM"],
    "Dendritic (DC)": ["CD207", "CD74", "RT1-DB1"]   
}

MARKERS_HUMAN = {
    "Keratinocytes (KC)": ["KRT14", "KRT5", "KRT1"],
    "Fibroblasts (FB)": ["DCN", "APOD", "COL1A1"],
    "Endothelial (EC)": ["CDH5", "CD93"],
    "Sweat_Gland (SGC)": ["DCD","AQP5"], 
    "Smooth_Muscle (SMC)": ["MYL9", "TAGLN", "ACTA2"],
    "Neural (NC)": ["CDH19", "SOX10", "PMP22"],
    "T_Cells (TC)": ["CD3D", "CD3E", "CD3G", "CCL5"],
    "NK_Cells (NK)": ["CD3D", "CD3E", "CD3G", "CCL5", "GZMB"],
    "Macrophages (Mø)": ["CD74", "AIF1", "CD68"],
    "Mast_Cells (MC)": ["TPSAB1", "TPSB2", "CPA3", "IL1RL1"],
    "Neutrophils (NEUT)": ["KRT14", "KRT5", "KRT1", "S100A8", "S100A9"]
}

# ============================================================
# 3. FIGURE 1: RAT ATLAS GENERATION (Time-series integration)
# ============================================================
print("\n--- FIGURE 1: RAT ANALYSIS ---")

# Load multiple time points to track radiation damage progression
rat_samples = [
    load_and_standardize("datas/rat/GSM5814220_con/", "Con", "Con", "rat"),
    load_and_standardize("datas/rat/GSM5814221_IR_7d/", "7d", "7d", "rat"),
    load_and_standardize("datas/rat/GSM5814222_IR_14d/", "14d", "14d", "rat"),
    load_and_standardize("datas/rat/GSM5814223_IR_28d/", "28d", "28d", "rat"),
]

# Concatenate all time-points and run the analysis pipeline
adata_rat = sc.concat([r for r in rat_samples if r is not None], label="batch")
adata_rat = process_basic(adata_rat, batch_key="batch")
adata_rat = auto_annotate_cell_types(adata_rat, MARKERS_RAT)

# --- Plot 1: Global t-SNE (Integrated) ---
sc.pl.tsne(adata_rat, color=["cell_type"], frameon=False, show=False, title="All Conditions Combined")
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Rat/tSNE_Combined.png", bbox_inches="tight")
plt.close()

# --- Plot 2: Split t-SNE (Condition specific) ---
plot_split_tsne(
    adata_rat, 
    keys=["Con", "7d", "14d", "28d"], 
    key_col="condition", 
    color_col="cell_type", 
    save_path=RESULT_DIR / "1_preprocessing_and_clustering_Rat/tSNE_Split.png"
)

# --- Plot 3: Marker Gene Dotplot ---
all_rat_markers = [g for list in MARKERS_RAT.values() for g in list]
all_rat_markers = list(dict.fromkeys(all_rat_markers)) # Remove duplicates
sc.pl.dotplot(adata_rat, all_rat_markers, groupby="cell_type", standard_scale="var", show=False)
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Rat/Dotplot.png", bbox_inches="tight")
plt.close()

# --- Plot 4: Cell Proportion Bar Chart (Shift in populations) ---
counts = adata_rat.obs.groupby(["condition", "cell_type"]).size().unstack(fill_value=0)
props = counts.div(counts.sum(axis=1), axis=0) * 100
props = props.loc[["Con", "7d", "14d", "28d"]]
props.plot(kind="bar", stacked=True, figsize=(8, 5), colormap="tab20")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("Cell Type Proportions")
plt.tight_layout()
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Rat/Proportions.png")
plt.close()

# CRITICAL STEP: Isolate Fibroblasts to pass to Phase 2 (Trajectory Analysis)
rat_fibro = adata_rat[adata_rat.obs["cell_type"].str.contains("Fibro")].copy()

# ============================================================
# 4. FIGURE 2: HUMAN ATLAS GENERATION (Control vs Treated)
# ============================================================
print("\n--- FIGURE 2: HUMAN ANALYSIS ---")

# Load Human Datasets
human_samples = [
    load_and_standardize("datas/human/GSM5821748_con/", "Con", "Con", "human"),
    load_and_standardize("datas/human/GSM5821749_IR/", "IR", "IR", "human"),
]

# Concatenate and Process
adata_human = sc.concat([h for h in human_samples if h is not None], label="batch")
adata_human = process_basic(adata_human, batch_key="batch")
adata_human = auto_annotate_cell_types(adata_human, MARKERS_HUMAN)

# --- Plot 1: Global t-SNE ---
sc.pl.tsne(adata_human, color=["cell_type"], frameon=False, show=False, title="All Conditions Combined")
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Human/tSNE_Combined.png", bbox_inches="tight")
plt.close()

# --- Plot 2: Split t-SNE ---
plot_split_tsne(
    adata_human, 
    keys=["Con", "IR"], 
    key_col="condition", 
    color_col="cell_type", 
    save_path=RESULT_DIR / "1_preprocessing_and_clustering_Human/tSNE_Split.png"
)

# --- Plot 3: Marker Gene Dotplot (Includes NUR77/NR4A1) ---
all_human_markers = [g for list in MARKERS_HUMAN.values() for g in list]
all_human_markers.append("NR4A1") 
all_human_markers = list(dict.fromkeys(all_human_markers))
sc.pl.dotplot(adata_human, all_human_markers, groupby="cell_type", standard_scale="var", show=False)
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Human/Dotplot.png", bbox_inches="tight")
plt.close()

# --- Plot 4: Cell Proportion Bar Chart ---
counts_human = adata_human.obs.groupby(["condition", "cell_type"]).size().unstack(fill_value=0)
props_human = counts_human.div(counts_human.sum(axis=1), axis=0) * 100
props_human = props_human.loc[["Con", "IR"]] 
props_human.plot(kind="bar", stacked=True, figsize=(6, 5), colormap="tab20")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("Human Cell Type Proportions")
plt.tight_layout()
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Human/Proportions.png")
plt.close()

# CRITICAL STEP: Isolate Human Fibroblasts to pass to Phase 2
human_fibro = adata_human[adata_human.obs["cell_type"].str.contains("Fibro")].copy()