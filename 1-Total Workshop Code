import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import seaborn as sns
from pathlib import Path
import bbknn

if not hasattr(scipy.sparse.csr_matrix, "ravel"):
    def ravel_fix(self, order='C'):
        return self.toarray().ravel(order)
    scipy.sparse.csr_matrix.ravel = ravel_fix

sc.settings.verbosity = 3 
sc.set_figure_params(dpi=300, fontsize=8, facecolor="white", frameon=False, vector_friendly=True)

BASE_DIR = Path.cwd()
WRITE_DIR = BASE_DIR / "write"
RESULT_DIR = Path.home() / "Adebali-Lab-Internship" / "results"

for d in [WRITE_DIR, RESULT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

for fig in ["1_preprocessing_and_clustering_Rat", "1_preprocessing_and_clustering_Human", "2_trajectory_analysis"]:
    (RESULT_DIR / fig).mkdir(exist_ok=True)

print(f"Results will be saved to: {RESULT_DIR}")

def load_and_standardize(path, sample_id, condition, species):
   try:
        adata = sc.read_10x_mtx(BASE_DIR / path, var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        adata.var_names = adata.var_names.str.upper() 
        adata.obs["sample_id"] = sample_id
        adata.obs["condition"] = condition
        adata.obs["species"] = species
        return adata
    except Exception as e:
        print(f"ERROR: Could not load {path} -> {e}")
        return None

def process_basic(adata, batch_key):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    
    print(f"Number of cells before filtering: {adata.n_obs}")

    print(f"Generating RAW QC Plots for {batch_key}...")

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False
    )
    plt.suptitle(f"Raw (Pre-Filter): {batch_key}", y=1.05) 
    plt.savefig(RESULT_DIR / f"QC_1_Raw_Violin_{batch_key}.png", bbox_inches="tight")
    plt.close()

    sc.pl.scatter(
        adata, 
        x="total_counts", 
        y="n_genes_by_counts", 
        color="pct_counts_mt", 
        show=False,
        title=f"Raw (Pre-Filter) Scatter: {batch_key}" 
    )
    plt.savefig(RESULT_DIR / f"QC_1_Raw_Scatter_{batch_key}.png", bbox_inches="tight")
    plt.close()

    adata = adata[
        (adata.obs.total_counts > 800) & 
        (adata.obs.total_counts < 20000) & 
        (adata.obs.pct_counts_mt < 15) & 
        (adata.obs.n_genes_by_counts > 400) 
    ].copy()
    print(f"Number of cells after filtering: {adata.n_obs}")

    print(f"Generating FILTERED QC Plots for {batch_key}...")

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False
    )
    plt.suptitle(f"Filtered (Post-Filter): {batch_key}", y=1.05)
    plt.savefig(RESULT_DIR / f"QC_2_Filtered_Violin_{batch_key}.png", bbox_inches="tight")
    plt.close()

    sc.pl.scatter(
        adata, 
        x="total_counts", 
        y="n_genes_by_counts", 
        color="pct_counts_mt", 
        show=False,
        title=f"Filtered (Post-Filter) Scatter: {batch_key}"
    )
    plt.savefig(RESULT_DIR / f"QC_2_Filtered_Scatter_{batch_key}.png", bbox_inches="tight")
    plt.close()

    print("Generating QC Plot: Log1P ONLY (No Normalization)...")

    adata_log_only = adata.copy()
    sc.pp.log1p(adata_log_only) 
    sc.tl.pca(adata_log_only)

    sc.pl.pca(
        adata_log_only,
        color=[batch_key, batch_key, "total_counts", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
        show=False,
        title=[
            f"LogOnly: {batch_key} (PC1-2)", 
            f"LogOnly: {batch_key} (PC3-4)", 
            "LogOnly: Counts", 
            "LogOnly: MT%"
        ]
    )
    plt.savefig(RESULT_DIR / f"QC_PCA_ScenarioA_LogOnly_{batch_key}.png", bbox_inches="tight")
    plt.close()
    del adata_log_only

    print("Generating QC Plot: Normalization ONLY (No Log1P)...")

    adata_norm_only = adata.copy()
    sc.pp.normalize_total(adata_norm_only, target_sum=1e4) 
    sc.tl.pca(adata_norm_only)

    sc.pl.pca(
        adata_norm_only,
        color=[batch_key, batch_key, "total_counts", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
        show=False,
        title=[
            f"NormOnly: {batch_key} (PC1-2)", 
            f"NormOnly: {batch_key} (PC3-4)", 
            "NormOnly: Counts", 
            "NormOnly: MT%"
        ]
    )
    plt.savefig(RESULT_DIR / f"QC_PCA_ScenarioB_NormOnly_{batch_key}.png", bbox_inches="tight")
    plt.close()
    del adata_norm_only

    sc.pp.normalize_total(adata, target_sum=1e4)
 
    sc.pp.log1p(adata)
    adata.raw = adata 
   
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata)
 
    sc.tl.pca(adata, mask_var="highly_variable")

    sc.pl.pca(
        adata,
        color=[batch_key, batch_key, "total_counts", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
        show=False,
        title=[
            f"Final: {batch_key} (PC1-2)", 
            f"Final: {batch_key} (PC3-4)", 
            "Final: Counts", 
            "Final: MT%"
        ]
    )
    plt.savefig(RESULT_DIR / f"QC_PCA_Standard_Final_{batch_key}.png", bbox_inches="tight")
    plt.close()

    bbknn.bbknn(adata, batch_key=batch_key)
    print("Calculating t-SNE...")
    sc.tl.tsne(adata, n_pcs=30)
    print("Calculating UMAP...")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=2)
    
    return adata

def auto_annotate_cell_types(adata, marker_dict):
    old_scores = [c for c in adata.obs.columns if c.startswith("score_")]
    if old_scores: adata.obs.drop(columns=old_scores, inplace=True)
    
    print("\n--- Marker Gene Check ---")

    for cell_type, genes in marker_dict.items():
        valid_genes = [g for g in genes if g in adata.var_names]
        if len(valid_genes) > 0:
            sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f"score_{cell_type}")
        else:
            print(f"WARNING: NO marker genes found for {cell_type}!")
    
    score_cols = [col for col in adata.obs.columns if col.startswith("score_")]
    if not score_cols: return adata

    cluster_scores = adata.obs.groupby("leiden")[score_cols].mean()
    cluster_mapping = {}
    for cluster in cluster_scores.index:
        best_score_col = cluster_scores.loc[cluster].idxmax()
        cell_name = best_score_col.replace("score_", "")
        cluster_mapping[cluster] = cell_name
    
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_mapping)
    return adata

def plot_split_tsne(adata, keys, key_col, color_col, save_path, custom_titles=None, figsize_per_plot=(4, 4)):
    n_plots = len(keys)
    fig, axes = plt.subplots(1, n_plots, figsize=(figsize_per_plot[0] * n_plots, figsize_per_plot[1]))
    
    if n_plots == 1: axes = [axes] 
 
    x_min, x_max = adata.obsm['X_tsne'][:, 0].min(), adata.obsm['X_tsne'][:, 0].max()
    y_min, y_max = adata.obsm['X_tsne'][:, 1].min(), adata.obsm['X_tsne'][:, 1].max()
    pad_x = (x_max - x_min) * 0.05
    pad_y = (y_max - y_min) * 0.05

    for i, key in enumerate(keys):
        ax = axes[i]
        subset = adata[adata.obs[key_col] == key] 

        if custom_titles and len(custom_titles) == len(keys):
            plot_title = custom_titles[i]
        else:
            plot_title = key
      
        sc.pl.tsne(
            subset, 
            color=color_col, 
            ax=ax, 
            show=False, 
            title=plot_title, 
            frameon=False,
            legend_loc="none" if i < n_plots - 1 else "right margin", 
            s=50 
        )
       
        ax.set_xlim(x_min - pad_x, x_max + pad_x)
        ax.set_ylim(y_min - pad_y, y_max + pad_y)

        ax.set_xlabel("tSNE_1")
        if i == 0: ax.set_ylabel("tSNE_2")
        else: ax.set_ylabel("")
            
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    print(f"Saved split t-SNE to: {save_path}")

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

print("\n--- FIGURE 1: RAT ANALYSIS ---")

rat_samples = [
    load_and_standardize("datas/rat/GSM5814220_con/", "Con", "Con", "rat"),
    load_and_standardize("datas/rat/GSM5814221_IR_7d/", "7d", "7d", "rat"),
    load_and_standardize("datas/rat/GSM5814222_IR_14d/", "14d", "14d", "rat"),
    load_and_standardize("datas/rat/GSM5814223_IR_28d/", "28d", "28d", "rat"),
]

adata_rat = sc.concat([r for r in rat_samples if r is not None], label="batch")
adata_rat = process_basic(adata_rat, batch_key="batch")
adata_rat = auto_annotate_cell_types(adata_rat, MARKERS_RAT)

sc.pl.tsne(adata_rat, color=["cell_type"], frameon=False, show=False, title="All Conditions Combined - RAT")
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Rat/tSNE_Combined_RAT.png", bbox_inches="tight")
plt.close()

plot_split_tsne(
    adata_rat, 
    keys=["Con", "7d", "14d", "28d"], 
    key_col="condition", 
    color_col="cell_type",
    custom_titles=["Control Group", "7 Days Post-IR", "14 Days Post-IR", "28 Days Post-IR"], 
    save_path=RESULT_DIR / "1_preprocessing_and_clustering_Rat/tSNE_Split_RAT.png"   
)

all_rat_markers = [g for list in MARKERS_RAT.values() for g in list]
all_rat_markers = list(dict.fromkeys(all_rat_markers)) 
sc.pl.dotplot(adata_rat, all_rat_markers, groupby="cell_type", standard_scale="var", show=False)
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Rat/Dotplot_RAT.png", bbox_inches="tight")
plt.close()

counts = adata_rat.obs.groupby(["condition", "cell_type"]).size().unstack(fill_value=0)
props = counts.div(counts.sum(axis=1), axis=0) * 100
props = props.loc[["Con", "7d", "14d", "28d"]]
props.plot(kind="bar", stacked=True, figsize=(8, 5), colormap="tab20")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("Cell Type Proportions")
plt.tight_layout()
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Rat/Proportions_RAT.png")
plt.close()

print("   > Extracting Rat Fibroblasts...")
rat_fibro = adata_rat[adata_rat.obs["cell_type"].str.contains("Fibro", na=False)].copy()
if rat_fibro.raw is not None: rat_fibro = rat_fibro.raw.to_adata()

print("\n--- FIGURE 2: HUMAN ANALYSIS ---")

human_samples = [
    load_and_standardize("datas/human/GSM5821748_con/", "Con", "Con", "human"),
    load_and_standardize("datas/human/GSM5821749_IR/", "IR", "IR", "human"),
]

adata_human = sc.concat([h for h in human_samples if h is not None], label="batch")
adata_human = process_basic(adata_human, batch_key="batch")
adata_human = auto_annotate_cell_types(adata_human, MARKERS_HUMAN)

sc.pl.tsne(adata_human, color=["cell_type"], frameon=False, show=False, title="All Conditions Combined - HUMAN")
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Human/tSNE_Combined_HUMAN.png", bbox_inches="tight")
plt.close()

plot_split_tsne(
    adata_human, 
    keys=["Con", "IR"], 
    key_col="condition", 
    color_col="cell_type", 
    custom_titles=["Control Group", "IR Treatment"], 
    save_path=RESULT_DIR / "1_preprocessing_and_clustering_Human/tSNE_Split_HUMAN.png"
)

all_human_markers = [g for list in MARKERS_HUMAN.values() for g in list] 
all_human_markers = list(dict.fromkeys(all_human_markers))
sc.pl.dotplot(adata_human, all_human_markers, groupby="cell_type", standard_scale="var", show=False)
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Human/Dotplot_HUMAN.png", bbox_inches="tight")
plt.close()

counts_human = adata_human.obs.groupby(["condition", "cell_type"]).size().unstack(fill_value=0)
props_human = counts_human.div(counts_human.sum(axis=1), axis=0) * 100
props_human = props_human.loc[["Con", "IR"]] 
props_human.plot(kind="bar", stacked=True, figsize=(6, 5), colormap="tab20")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("Human Cell Type Proportions")
plt.tight_layout()
plt.savefig(RESULT_DIR / "1_preprocessing_and_clustering_Human/Proportions_HUMAN.png")
plt.close()

print("   > Extracting Human Fibroblasts...")
human_fibro = adata_human[adata_human.obs["cell_type"].str.contains("Fibro", na=False)].copy()
if human_fibro.raw is not None: human_fibro = human_fibro.raw.to_adata()

print("\n--- FIGURE 3: FIBROBLAST DETAILED ANALYSIS ---")

def auto_tune_resolution(adata, target_clusters):
    res_low, res_high = 0.1, 2.0
    for _ in range(15): 
        current_res = (res_low + res_high) / 2
        sc.tl.leiden(adata, resolution=current_res, key_added="sub_leiden")
        n_clusters = len(adata.obs['sub_leiden'].unique())
        
        if n_clusters == target_clusters: return current_res
        if n_clusters > target_clusters: res_high = current_res 
        else: res_low = current_res 
    return current_res

def analyze_fibroblast_subtypes_branched(adata, species_name):
    
    print(f"   > Processing {species_name} Fibroblasts...")
    
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata 
    
    sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor='seurat')
    sc.pp.scale(adata)
    sc.pp.pca(adata)

    sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_pca") 
  
    target = 7 if species_name == "Rat" else 5
    res = auto_tune_resolution(adata, target)
    sc.tl.leiden(adata, resolution=res, key_added="sub_leiden")
   
    cluster_map = {str(c): f"FB{i+1}" for i, c in enumerate(sorted(adata.obs['sub_leiden'].unique().astype(int)))}
    adata.obs["sub_type"] = adata.obs["sub_leiden"].map(cluster_map)
   
    sc.tl.paga(adata, groups="sub_type")
    sc.pl.paga(adata, plot=False, threshold=0.15) 

    print(f"     > Calculating branched tree structure...")
    layout_key = 'fr' 
    if 'fa2' in locals() or 'fa2' in globals():
        try:
            sc.tl.draw_graph(adata, init_pos='paga', layout='fa', random_state=42)
            layout_key = 'fa'
        except:
            sc.tl.draw_graph(adata, init_pos='paga', layout='fr')
    else:
        sc.tl.draw_graph(adata, init_pos='paga', layout='fr')

    con_counts = adata.obs.groupby("sub_type")["condition"].apply(lambda x: (x == "Con").sum())
    root_clust = con_counts.idxmax()
    root_cells = np.flatnonzero(adata.obs["sub_type"] == root_clust)
    adata.uns["iroot"] = root_cells[0] if len(root_cells) > 0 else 0 

    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
  
    sc.tl.tsne(adata, n_pcs=20)
    
    return adata

rat_fibro = analyze_fibroblast_subtypes_branched(rat_fibro, "Rat")
human_fibro = analyze_fibroblast_subtypes_branched(human_fibro, "Human")

print("\n--- Generating Figure 3 Plots ---")

RAT_TFS = ["BACH1", "ETS1", "NFKB1", "RFX5", "SMAD1", "STAT2", "ATF3", "FOS", "FOSB", 
           "JUN", "JUND", "IRF7", "IRF1", "CEBPE", "HOXD4", "DLX5", "GRHL2", "NFE2L1", "POU3F1", "SREBF2"]

HUMAN_TFS = ["EHF", "EMX2", "FOXA1", "HOXB3", "HOXB4", "HOXC6", "HOXC8", "HOXC9", "IRF6", 
             "MESP1", "FOSL1", "NFE2L2", "NR2F1", "RARB", "SREBF2", "TEAD4", "MAFA", "NFATC1", "LHX9", "ETV3"]

def plot_fig3a(adata, genes, filename, order, plot_title):
    valid = [g for g in genes if g in adata.raw.var_names]
    present = [c for c in order if c in adata.obs["condition"].unique()]
    adata.obs["condition"] = adata.obs["condition"].astype("category").cat.reorder_categories(present)
    
    if valid:
        sc.pl.dotplot(adata, valid, groupby="condition", standard_scale="var", show=False, title=plot_title)
        plt.savefig(RESULT_DIR / f"2_trajectory_analysis/{filename}", bbox_inches="tight"); plt.close()

plot_fig3a(rat_fibro, RAT_TFS, "Rat_TF.png", ["Con", "7d", "14d", "28d"], "Rat Fibroblast TF Enrichment")
plot_fig3a(human_fibro, HUMAN_TFS, "Human_TF.png", ["Con", "IR"], "Human Fibroblast TF Enrichment")

tsne_rat = sc.pl.tsne(rat_fibro, color="sub_type", legend_loc="on data", palette="tab10", 
                      title="Rat Fibroblasts", return_fig=True, show=False)
tsne_rat.savefig(RESULT_DIR / "2_trajectory_analysis/Rat_tSNE.png", bbox_inches="tight")
plt.close(tsne_rat)

plot_split_tsne(rat_fibro, keys=["Con", "7d", "14d", "28d"], key_col="condition", color_col="sub_type",
                custom_titles=["Control Group", "7 Days Post-IR", "14 Days Post-IR", "28 Days Post-IR"], 
                save_path=RESULT_DIR / "2_trajectory_analysis/Rat_tSNE_Split.png")

tsne_human = sc.pl.tsne(human_fibro, color="sub_type", legend_loc="on data", palette="tab10", 
                        title="Human Fibroblasts", return_fig=True, show=False)
tsne_human.savefig(RESULT_DIR / "2_trajectory_analysis/Human_tSNE.png", bbox_inches="tight")
plt.close(tsne_human)

plot_split_tsne(human_fibro, keys=["Con", "IR"], key_col="condition", color_col="sub_type", 
                custom_titles=["Control Group", "IR Treatment"], 
                save_path=RESULT_DIR / "2_trajectory_analysis/Human_tSNE_Split.png")

def plot_trajectory_branched(adata, name, label):
    layout_key = 'fa' if 'X_draw_graph_fa' in adata.obsm else 'fr'
    print(f"   > Plotting {name} branched trajectory using '{layout_key}'...")
    
    gene = "NR4A1" if "NR4A1" in adata.raw.var_names else "NUR77"

    sc.pl.draw_graph(adata, color="dpt_pseudotime", layout=layout_key, 
                     frameon=False, show=False, cmap="viridis", 
                     edges=True, edges_width=0.2, 
                     title=f"{name} Pseudotime (Trajectory)")
    plt.savefig(RESULT_DIR / f"2_trajectory_analysis/{label}_{name}_Pseudotime_Trajectory.png", bbox_inches="tight")
    plt.close()
    
    if gene in adata.raw.var_names:
        raw_data = adata.raw[:, gene].X
        if hasattr(raw_data, "toarray"): raw_data = raw_data.toarray()
        raw_data = raw_data.flatten()
  
        counts = np.expm1(raw_data)
        custom_val = np.log10(counts + 0.1)
    
        plot_col_name = f"{gene}_log10"
        adata.obs[plot_col_name] = custom_val

        sc.pl.draw_graph(adata, color=plot_col_name, layout=layout_key, 
                         frameon=False, show=False, cmap="viridis", 
                         edges=True, edges_width=0.2, 
                         title=f"{name} {gene} (log10(val+0.1))")
        plt.savefig(RESULT_DIR / f"2_trajectory_analysis/{label}_{name}_{gene}_Trajectory.png", bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(6, 4))
        sc.pl.scatter(adata, x='dpt_pseudotime', y=plot_col_name, color='sub_type', 
                      show=False, title=f"{name} {gene} across Pseudotime")
        plt.xlabel("Pseudotime (Start -> End)")
        plt.ylabel(f"{gene} Expression (log10(val+0.1))")
        plt.savefig(RESULT_DIR / f"2_trajectory_analysis/{label}_{name}_{gene}_vs_Pseudotime_Scatter.png", bbox_inches="tight")
        plt.close()

    sc.pl.draw_graph(adata, color="sub_type", layout=layout_key, 
                     frameon=False, show=False, palette="tab10", 
                     edges=True, edges_width=0.3, 
                     title=f"{name} Subtype Branches")
    plt.savefig(RESULT_DIR / f"2_trajectory_analysis/{label}_{name}_Subtypes_Trajectory.png", bbox_inches="tight"); plt.close()

plot_trajectory_branched(rat_fibro, "Rat", "Fig3E")
plot_trajectory_branched(human_fibro, "Human", "Fig3F")
    
print("\n--- ALL ANALYSES COMPLETE ---")
