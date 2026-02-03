# Single-Cell-RNA-Sequencing-Tutorial
A comprehensive single-cell RNA-seq (scRNA-seq) analysis pipeline in Python. Includes quality control, batch correction (BBKNN), subclustering, and trajectory inference (PAGA/FA2) using Scanpy.

🧬 Single-Cell RNA-Seq Analysis: Radiation-Induced Skin Fibrosis & Nur77/NR4A1 Dynamics
This repository contains the analysis of Single-Cell RNA Sequencing (scRNA-seq) data obtained from radiation-exposed rat and human skin samples. The study investigates the trajectory changes of fibroblast subtypes during the development of radiation-induced fibrosis and the role of the Nur77 (NR4A1) gene in this process.

📖 Table of Contents
Introduction: What is scRNA-seq?

Dataset (Data Availability)

Requirements & Setup (Environment)

Libraries Used

Code Workflow and Usage

Results and Outputs

🔬 1. Introduction: What is Single-Cell RNA Sequencing (scRNA-seq)? <a name="what-is-scrna-seq"></a>
Traditional RNA sequencing (Bulk RNA-seq) gives the average gene expression of millions of cells in a tissue (like a "fruit smoothie"). Single-Cell RNA Sequencing (scRNA-seq), on the other hand, tags and analyzes the genetic material of each cell individually (like identifying each "fruit" separately).

Why use it?

To discover rare cell types in the tissue.

To track the developmental changes of cells over time (Pseudotime/Trajectory).

As in this project, to see at cellular resolution how radiation causes fibroblasts in the skin tissue to differentiate.

📂 2. Dataset (Data Availability) <a name="dataset"></a>
The analyses were performed on publicly available GEO (Gene Expression Omnibus) datasets generated using 10x Genomics technology. For the code to run, the data must be downloaded into the datahr/ folder.

Rat Data: Control, 7-Day, 14-Day, and 28-Day radiation (IR) samples.

GSE193836 (GSM5814220-GSM5814223)

Human Data: Control and Radiation (IR) samples.

GSE194121 (GSM5821748-GSM5821749)

💻 3. Requirements & Setup (Environment) <a name="environment"></a>
Using a Conda environment is highly recommended for the reproducibility of the analyses. Python 3.8 or higher is required.

Bash

# 1. Create a new conda environment
conda create -n sc_nur77_env python=3.9
conda activate sc_nur77_env

# 2. Install necessary libraries
pip install scanpy pandas numpy matplotlib seaborn bbknn
pip install fa2  # Critical for Trajectory tree plots (ForceAtlas2)!
🛠️ 4. Libraries Used <a name="libraries-used"></a>
Scanpy: The main library for scRNA-seq data analysis (Quality control, normalization, dimensionality reduction).

BBKNN: For "Batch Effect" correction between different biological samples (e.g., prevents control and irradiated samples from artificially separating on the graph).

fa2 (ForceAtlas2): To visualize the biological trajectories of cells as a branched tree in fibroblast branching analyses.

Pandas & Numpy: For data manipulation and statistical operations.

Matplotlib & Seaborn: For data visualization (t-SNE, UMAP, DotPlot, Violin plots).

🚀 5. Code Workflow and Usage <a name="code-workflow"></a>
The project consists of two main Python files. They should be run sequentially from the command line.

Phase 1: Identification of General Cell Types
File: Nur77 Figures 1&2.py

Operation: Loads Rat and Human data, performs quality control (QC) (<15% MT gene ratio, 800-20000 total counts). Applies batch correction with BBKNN and clusters cells using t-SNE algorithms. Performs auto-annotation using gene marker lists (Keratinocytes, Fibroblasts, T Cells, etc.).

Run: python "Nur77 Figures 1&2.py"

Phase 2: Fibroblast Subtypes and Trajectory (Developmental Process)
File: Nur77 Figure 3.py

Operation: Extracts the Fibroblast cells isolated in the first phase. Identifies their subtypes. Calculates how cells differentiate from the control group towards the radiation group (Pseudotime) using PAGA and Diffusion Map. Examines Transcription Factors (TF).

Run: python "Nur77 Figure 3.py"

📊 6. Results and Outputs <a name="results"></a>
When the codes run successfully, the following figures will be generated in your output folder (default: ~/Adebali-Lab-Internship/):

Figure 1 & 2: Rat and Human Skin Atlas
General maps (t-SNE) where cells are grouped by their biological identities. The change rates in cell populations post-radiation are graphed.

Figure 1D/2_Split: Side-by-side comparison of control and radiation groups.

Dotplot Graphs: Expression levels of genes specific to each cell type (e.g., KRT14 for Keratinocytes, DCN for Fibroblasts).

Figure 3: Branching Analysis of Fibroblasts
This is the core part of the study. It examines how fibroblasts transform into fibrotic cells with radiation.

3A (TF Enrichment): Representation of Transcription Factors (BACH1, ETS1, etc.) activated at different time points.

3D (Nur77/NR4A1 Expression): Violin plots showing a significant increase of the Nur77 gene in irradiated fibroblasts.

3E & 3F (Pseudotime Trajectory): "Cell Developmental Tree" drawn with the ForceAtlas2 algorithm. This graph shows how cells start from the Control state (root) and differentiate into new subtypes over time (Day 7 -> Day 14 -> Day 28) under the effect of radiation. The expression of the Nur77 gene in this process is reflected as a gradient.

📝 Notes
Auto-Tune Resolution: The fibroblast subtyping in Figure 3 works with a dynamic "binary search" logic to allow the algorithm to automatically find the desired number of clusters (Rat: 7, Human: 5).

Preservation of Rare Cells: The Leiden clustering resolution is kept high (res=2.0) so that rare cell types (Mast cells, neutrophils, etc.) are not lost during filtering.

Contact / Author: These codes were optimized by [Your Name/Project Team] as part of the [Lab Name] internship.

💡 (Gemini Tip for GitHub): Embedding Images (Optional)
Visuals make GitHub READMEs much easier to understand. Once the analysis is complete, I recommend adding a few of the generated PNG files directly to your README.

Example Markdown code for adding an image:

Markdown

### Example Output: Rat Fibroblast Trajectory Analysis
![Rat Fibroblast Trajectory](Fig3_Fibroblast/Fig3E_Rat_Pseudotime_Trajectory.png)
