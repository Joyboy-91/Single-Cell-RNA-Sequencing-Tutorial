# 🧬 Comprehensive Single-Cell RNA-Seq Analysis Pipeline

![Python](https://img.shields.io/badge/python-3.10.19-blue.svg)
![Scanpy](https://img.shields.io/badge/Scanpy-1.11+-green.svg)

This repository provides an end-to-end analytical pipeline for **Single-Cell RNA Sequencing (scRNA-seq)** data. It includes scripts for quality control, batch correction, dimensionality reduction, cell type annotation, and advanced trajectory inference (pseudotime) to study cellular differentiation and dynamics.

---

## 📖 Table of Contents
1. [Introduction: What is scRNA-seq?](#1-introduction)
2. [Dataset & Inputs](#2-dataset-inputs)
3. [Requirements & Setup (Environment)](#3-environment)
4. [Libraries Used](#4-libraries)
5. [Code Workflow and Usage](#5-workflow)
6. [Results and Outputs](#6-outputs)
7. [Notes & Acknowledgements](#7-notes)

---

## 🔬 1. Introduction: What is Single-Cell RNA Sequencing? <a name="1-introduction"></a>

Traditional RNA sequencing (Bulk RNA-seq) gives the average gene expression of millions of cells in a tissue (like a "fruit smoothie"). **Single-Cell RNA Sequencing (scRNA-seq)**, on the other hand, tags and analyzes the genetic material of each cell individually (like identifying each "fruit" separately).



**Why use this pipeline?**
* To discover **rare cell types** and hidden sub-populations in complex tissues.
* To perform robust **batch correction** across multiple samples or experimental conditions (e.g., Control vs. Treated).
* To track the developmental or disease-driven changes of cells over time (**Pseudotime/Trajectory Analysis**).

---

## 📂 2. Dataset & Inputs <a name="2-dataset-inputs"></a>

This pipeline is compatible with standard 10x Genomics output formats. To run the code, your raw data files must be placed inside the `data/` directory.

**Supported Input Formats:**
* **10x CellRanger Outputs:** `matrix.mtx`, `barcodes.tsv`, `features.tsv` (or `genes.tsv`).
* **AnnData:** Pre-processed `.h5ad` files.

You can customize sample names and condition labels directly within the configuration section of the scripts.

---

## 💻 3. Requirements & Setup (Environment) <a name="3-environment"></a>

Using a **Conda** environment is highly recommended for the reproducibility of the analyses. Python 3.10 or higher is required.

```bash
# 1. Create a new conda environment
conda create -n sc_pipeline_env python=3.10
conda activate sc_pipeline_env

# 2. Install necessary libraries
pip install scanpy pandas numpy matplotlib seaborn bbknn
pip install fa2-modified  # Critical for robust Trajectory tree plots (ForceAtlas2)
```

---

## 🛠️ 4. Libraries Used <a name="4-libraries"></a>

**Scanpy (v1.11.5)** The core analytical toolkit used for the entire single-cell pipeline. Handled quality control, normalizations, Leiden clustering, and advanced trajectory inference tools including **PAGA** (topological skeleton) and **Diffusion Maps** (pseudotime calculation).

**BBKNN (v1.6.0)** A fast and lightweight graph-based batch effect correction algorithm. It integrates control and irradiated datasets by altering the neighborhood graph without altering the PCA space.

**leidenalg (v0.11.0)** The backend algorithm used by Scanpy to perform graph-based community detection, allowing for high-resolution sub-clustering of specific cell types (e.g., Fibroblasts).

**igraph (v1.0.0)** A high-performance network analysis library acting as the computational backend for both **PAGA** and **ForceAtlas2** to calculate complex topological structures and cellular connectivity.

**fa2-modified (v0.4)** A highly optimized graph layout engine used to visualize the PAGA-derived branching trajectories. *Note: The pipeline uses a custom `sys.modules` monkey patch to dynamically load the more stable `fa2_modified` branch, ensuring compatibility with the latest Scanpy versions, while keeping standard `fa2` as a fallback.*

**Pandas (v2.2.3) & NumPy (v2.2+)** Core Python libraries used for metadata manipulation, sparse matrix operations, and numerical computations (e.g., root cell identification for pseudotime).

**Matplotlib (v3.10.0) & Seaborn (v0.13.2)** Used to generate all publication-ready visualizations, including t-SNE split plots, gene expression dot plots, and pseudotime distribution scatter plots.

**pathlib (Built-in, Python 3.10)** A standard Python library used for object-oriented filesystem path manipulation. It ensures cross-platform compatibility (Windows, macOS, Linux) when loading raw data and automatically creates the necessary `results/` output directories.

**sys (Built-in, Python 3.10)** Used for system-specific parameters and functions. In this pipeline, it is specifically utilized for "monkey patching" `sys.modules` to dynamically redirect standard `fa2` calls to the `fa2_modified` library.

---

## 🚀 5. Code Workflow and Usage <a name="5-workflow"></a>

The workflow is divided into **two sequential Python scripts**.

### 🔹 Phase 1: Preprocessing & Global Cell Clustering

**File:** [`1_preprocessing_and_clustering.py`](1_preprocessing_and_clustering.py) [`2_trajectory_analysis.py`](2_trajectory_analysis.py)

**Operation:** The pipeline loads raw single-cell data, performs rigorous quality control, normalizes expression, and corrects batch effects using **BBKNN**. It maps global cell populations via high-resolution **Leiden** clustering and **t-SNE/UMAP** embeddings, auto-annotating cell types with marker genes. After isolating a specific target subset, the workflow executes deep sub-clustering and infers developmental/disease trajectories from a root state using **PAGA** and **Diffusion Maps (Pseudotime)**. Finally, it visualizes cellular transition dynamics and gene expression changes through **ForceAtlas2** branched tree structures and condition-split scatter plots.

**Run:**
```bash
python 1_preprocessing_and_clustering.py
```

```bash
python 2_trajectory_analysis.py
```

--- 
## 📊 6. Results and Outputs <a name="6-outputs"></a>
Successful execution will automatically generate a results/ folder containing the following:

**Global Cell Atlas**
**t-SNE / UMAP Embeddings:** Visualizations showing the overall cellular diversity.

**Condition Split Plots:** Side-by-side comparisons of different experimental groups (e.g., Control vs. Treated).

**Marker Gene Dotplots:** A summary of top differentially expressed genes across all identified clusters.

**Developmental Trajectory**
**Transcription Factor (TF) Enrichment:** Plots showing active regulatory mechanisms.

**Expression Dynamics:** Violin plots showing the expression of your Genes of Interest across different conditions.

**Pseudotime Tree:** A continuous map (ForceAtlas2) showing how cells transition from State A to State B, colored by pseudotime gradients and gene expression.

--- 

## 📝 7. Notes & Acknowledgements <a name="7-notes"></a>

**Auto-Tune Resolution:** The sub-clustering algorithm features a dynamic "binary search" logic to automatically find the optimal Leiden resolution for a target number of clusters.

**Rare Cell Preservation:** Filtering parameters are optimized to preserve rare cell types that might otherwise be lost during strict quality control.

**Contact / Author:** Pipeline developed and maintained by [Onur Yıldırım](https://github.com/Joyboy-91) during an internship at the Adebali Lab, Sabancı University.
