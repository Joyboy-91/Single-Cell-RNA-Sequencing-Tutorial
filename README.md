Metninizdeki bozulan başlıkları, madde işaretlerini ve özellikle tüm kod bloklarını GitHub Markdown formatına uygun şekilde (siyah kod kutusu içinde görünecek şekilde) düzenledim.

Aşağıdaki metni olduğu gibi kopyalayıp README.md dosyanıza yapıştırabilirsiniz.

Markdown

# 🧬 Comprehensive Single-Cell RNA-Seq Analysis Pipeline

![Python](https://img.shields.io/badge/python-3.9-blue.svg)
![Scanpy](https://img.shields.io/badge/Scanpy-1.9+-green.svg)

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

Using a **Conda** environment is highly recommended for the reproducibility of the analyses. Python 3.8 or higher is required.

```bash
# 1. Create a new conda environment
conda create -n sc_pipeline_env python=3.9
conda activate sc_pipeline_env

# 2. Install necessary libraries
pip install scanpy pandas numpy matplotlib seaborn bbknn
pip install fa2  # Critical for robust Trajectory tree plots (ForceAtlas2)
```

---

## 🛠️ 4. Libraries Used <a name="4-libraries"></a>

**Scanpy**  
The core library for scRNA-seq data analysis (quality control, normalization, clustering, visualization).

**BBKNN**  
Used for batch effect correction to integrate multiple datasets while preserving biological variation.

**fa2 (ForceAtlas2)**  
Graph-based layout engine for visualizing developmental trajectories and pseudotime trees.

**Pandas & NumPy**  
Libraries for data manipulation, matrix operations, and numerical computation.

**Matplotlib & Seaborn**  
Used to generate high-quality, publication-ready visualizations.

---

## 🚀 5. Code Workflow and Usage <a name="5-workflow"></a>

The workflow is divided into **two sequential Python scripts**.

---

### 🔹 Phase 1: Preprocessing & Global Cell Clustering

**File:** `1_preprocessing_and_clustering.py`

**Operation:**  
Loads raw data, performs quality control (mitochondrial ratio, gene counts), normalizes expression values, applies **BBKNN** batch correction, performs clustering (t-SNE / UMAP), and generates marker gene dot plots for automated cell type annotation.

**Run:**
```bash
python 1_preprocessing_and_clustering.py
```
--- 
📊 6. Results and Outputs <a name="6-outputs"></a>
Successful execution will automatically generate a results/ folder containing the following:

Global Cell Atlas
t-SNE / UMAP Embeddings: Visualizations showing the overall cellular diversity.

Condition Split Plots: Side-by-side comparisons of different experimental groups (e.g., Control vs. Treated).

Marker Gene Dotplots: A summary of top differentially expressed genes across all identified clusters.

Developmental Trajectory
Transcription Factor (TF) Enrichment: Plots showing active regulatory mechanisms.

Expression Dynamics: Violin plots showing the expression of your Genes of Interest across different conditions.

Pseudotime Tree: A continuous map (ForceAtlas2) showing how cells transition from State A to State B, colored by pseudotime gradients and gene expression.

--- 

📝 7. Notes & Acknowledgements <a name="7-notes"></a>
Auto-Tune Resolution: The sub-clustering algorithm features a dynamic "binary search" logic to automatically find the optimal Leiden resolution for a target number of clusters.

Rare Cell Preservation: Filtering parameters are optimized to preserve rare cell types that might otherwise be lost during strict quality control.
