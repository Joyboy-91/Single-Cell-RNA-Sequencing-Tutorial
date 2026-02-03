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
* **To identify rare cell types:** Uses high-resolution clustering and tailored QC thresholds to discover hidden sub-populations (e.g., specific Fibroblast subtypes) in complex tissues.
* **To integrate complex experimental designs:** Performs robust **batch effect correction** across multiple species (Human/Rat) and diverse conditions (e.g., Control vs. Irradiated / Time-series data).
* **To track cellular fate decisions:** Uses advanced **Pseudotime and Trajectory Analysis** to map how cells differentiate or respond to disease states over time.

##How is the Raw Data Generated? (The 10x Genomics Workflow)##
Before running this code, the biological samples go through a fascinating micro-engineering process:
1. **Dissociation:** Complex tissues are enzymatically digested into a single-cell suspension.
2. **Droplet Encapsulation:** Each cell is trapped inside a microscopic oil droplet alongside a specialized Gel Bead.
3. **Barcoding:** RNA molecules are tagged with a unique **Cell Barcode** (to identify which cell they belong to) and a **UMI** (to count individual transcripts).
4. **Sequencing:** The tagged RNAs are amplified and sequenced via Illumina platforms.
5. **Matrix Generation:** Finally, the *CellRanger* pipeline maps these reads to a reference genome, producing the final `matrix.mtx`, `barcodes.tsv`, and `features.tsv` files that this Python pipeline uses.

---

## 📂 2. Dataset & Inputs <a name="2-dataset-inputs"></a>

This pipeline is fully compatible with standard 10x Genomics output formats. To run the code successfully, your raw data folders must be placed inside the `datas/` directory.

### 📥 Download the Dataset
The raw dataset used in this tutorial (515 MB uncompressed) is hosted in the GitHub Releases section to keep the repository lightweight. 

**Option A: One-Liner for Terminal (Linux / Mac / WSL) - Recommended**
Run the following commands in your terminal to download and extract the dataset automatically:
```bash
wget [https://github.com/Joyboy-91/Single-Cell-RNA-Sequencing-Tutorial/releases/download/v1.0.0/datas.tar.gz](https://github.com/Joyboy-91/Single-Cell-RNA-Sequencing-Tutorial/releases/download/v1.0.0/datas.tar.gz)
tar -xzvf datas.tar.gz
```

**Option B: Manual Download (Windows) `datas.tar.gz`]
(https://github.com/Joyboy-91/Single-Cell-RNA-Sequencing-Tutorial/releases/download/v1.0.0/datas.tar.gz)**

Place the downloaded file in your project folder and extract it so that the `datas/` folder is visible.

**To extract via Terminal (Linux/Mac/WSL/Git Bash for Windows):**
```bash
tar -xzvf datas.tar.
```

**Supported Input Formats:**
* **10x CellRanger Outputs:** The standard trio of `matrix.mtx`, `barcodes.tsv`, and `features.tsv` (or `genes.tsv`).
* **AnnData:** Pre-processed `.h5ad` files can also be integrated directly.

**Configuring Your Data:**
The pipeline uses a custom helper function (`load_and_standardize`) that automatically tags your datasets with specific metadata. You can easily define your own sample names, experimental conditions (e.g., "Control", "7d", "IR"), and species directly within the script's data loading section.

---

## 💻 3. Requirements & Setup (Environment) <a name="3-environment"></a>

Using a **Conda** environment is highly recommended for the reproducibility of the analyses. Python 3.10 or higher is required.

```bash
# 1. Create a new conda environment (Python 3.10 is recommended for package compatibility)
conda create -n sc_pipeline_env python=3.10
conda activate sc_pipeline_env

# 2. Install necessary libraries with specific versions for perfect reproducibility

# Option A: Pure PIP Installation
pip install scanpy==1.11.5 pandas==2.2.3 numpy==2.2.6 matplotlib==3.10.0 seaborn==0.13.2 bbknn==1.6.0 leidenalg==0.11.0 igraph==1.0.0
pip install fa2-modified==0.4 

# Option B: Conda + PIP Installation (Recommended for better stability with C++ backend)
conda install -c conda-forge scanpy==1.11.5 pandas==2.2.3 numpy==2.2.6 matplotlib==3.10.0 seaborn==0.13.2 leidenalg==0.11.0 python-igraph==1.0.0 -y
pip install bbknn==1.6.0 fa2-modified==0.4
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

**Files:** [`1_preprocessing_and_clustering.py`](1_preprocessing_and_clustering.py) [`2_trajectory_analysis.py`](2_trajectory_analysis.py)

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

Successful execution will automatically generate a `results/` folder containing the following high-quality, publication-ready visualizations:

### Global Cell Atlas
* **t-SNE / UMAP Embeddings:** Visualizations showing overall cellular diversity across all integrated datasets.
  
* **Condition-Split Plots:** Side-by-side t-SNE comparisons of different experimental groups (e.g., Control vs. Treated) on the same coordinate scale.
  
* **Cellular Proportions:** Stacked bar charts quantifying the percentage shifts in cell populations between conditions.
  
* **Marker Gene Dotplots:** A comprehensive summary of key differentially expressed genes used for automated cluster annotation.

### Developmental Trajectory
* **Transcription Factor (TF) Enrichment:** Dotplots showing active regulatory mechanisms and transcription factors across different time points/conditions.
  
* **Expression Dynamics:** Violin plots showing the expression levels of your *Genes of Interest* across specific cell subtypes.
  
* **Pseudotime Tree:** A continuous, branched map (ForceAtlas2) showing how cells transition from a Root state to terminal states, colored by pseudotime gradients and gene expression.
  
* **Expression vs. Pseudotime Scatter:** Dynamic scatter plots tracking how the expression of a specific gene changes as cells progress through the developmental timeline.

--- 

## 📝 7. Notes & Acknowledgements <a name="7-notes"></a>

* **Auto-Tune Resolution:** The sub-clustering algorithm (Phase 2) features a dynamic "binary search" logic to automatically find the optimal Leiden resolution required to achieve a user-defined target number of clusters.
  
* **Rare Cell Preservation:** During Phase 1, the Leiden clustering resolution is intentionally kept high (res=2.0), and specific QC thresholds are optimized to preserve rare cell types that might otherwise be filtered out.
  
* **Trajectory Branching Optimization:** The neighborhood graph parameters (e.g., `n_neighbors=10`) are fine-tuned to ensure the PAGA/ForceAtlas2 algorithms can accurately calculate structural branches rather than forcing cells into a single linear trajectory.

**Contact / Author:** Pipeline developed and maintained by [Onur Yıldırım](https://github.com/Joyboy-91) during an internship at the Adebali Lab, Sabancı University.
