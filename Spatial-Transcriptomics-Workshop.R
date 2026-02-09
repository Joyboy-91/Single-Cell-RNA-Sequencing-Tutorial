library(NanoStringNCTools)  # Core package for NanoString objects
library(GeomxTools)         # Specialized tools for GeoMx DSP data processing
library(dplyr)              # Data manipulation (filter, select, mutate)
library(ggplot2)            # Visualization
library(umap)               # Dimensionality reduction for visualization
library(lmerTest)           # Linear Mixed Models (critical for spatial data statistics)
library(ggrepel)            # Non-overlapping text labels in plots
library(reshape2)           # Reshaping data (wide to long format)

# -------------------------------------------------------------------------
# 1. FILE PATHS AND METADATA LOADING
# -------------------------------------------------------------------------
# PURPOSE: Define where raw data (DCC) and probe info (PKC) are located.
base_dir <- "C:/Users/ben10/Desktop"

dcc_dir <- file.path(base_dir, "GSE298108.dcc")
pkc_file <- file.path(base_dir, "GSE298108.pkc", "GSM9008473_Mm_R_NGS_WTA_v1.0.pkc")
meta_file <- file.path(base_dir, "GSE298108.metadata", "GSE298108.metadata.tsv")

# Load file lists and metadata table
DCCFiles <- dir(dcc_dir, pattern = ".dcc$", full.names = TRUE)
phenoData <- read.table(meta_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# --- NTC Correction ---
# METHOD: Negative Template Control (NTC) Handling.
# PURPOSE: Identify "background" noise. NTC samples measure non-specific binding.
# We explicitly label them to ensure they aren't mistaken for biological tissue.
phenoData[["slide name"]] <- "Slide1"
ntc_index <- which(phenoData$tissue == "No Template Control" | 
                     phenoData$genotype == "No Template Control")

if(length(ntc_index) > 0) {
  phenoData[ntc_index, "slide name"] <- "No Template Control"
  print(paste("NTC Sample Defined:", phenoData[ntc_index, "gsm"]))
}

# Create a unique identifier matching the raw file names
phenoData$dcc_filename <- paste(phenoData$gsm, phenoData$title, sep = "_")

# Check/Control
print(head(phenoData[, c("gsm", "title", "slide name", "dcc_filename")]))

# -------------------------------------------------------------------------
# 2. CREATING GEOMX OBJECT
# -------------------------------------------------------------------------

# METHOD: Object Construction
# PURPOSE: Aggregates Raw Counts (DCC), Gene Info (PKC), and Metadata (phenoData)
# into a single S4 object (GeoMxSet) for downstream analysis.
demoData <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = pkc_file,
                         phenoData = phenoData,
                         phenoDataDccColName = "dcc_filename",
                         protocolDataColNames = NULL,
                         experimentDataColNames = NULL
  )

# Check if successful
print(paste("Raw Data Sample Count:", dim(demoData)[2]))

# METHOD: Filtering Irrelevant Tissue
# PURPOSE: We remove 'adipose' tissue because comparing fat tissue to skin (epidermis)
# creates too much biological noise. We want to focus on the skin/immune response.
non_adipose <- !grepl("adipose", pData(demoData)$genotype, ignore.case = TRUE)
demoData <- demoData[, non_adipose]

print(paste("Remaining Sample Count After Removing Adipose:", dim(demoData)[2]))

# -------------------------------------------------------------------------
# 3. QUALITY CONTROL (QC)
# -------------------------------------------------------------------------
# METHOD: Technical QC Thresholding
# PURPOSE: Ensure data reliability. We filter out Low-Quality Regions of Interest (ROIs).
QC_params <-
  list(minSegmentReads = 1000,    # Minimum sequencing depth per ROI
       percentTrimmed = 80,       # Sequencing quality (adapter trimming)
       percentStitched = 80,      # Paired-end read merging success
       percentAligned = 80,       # Alignment to genome reference
       percentSaturation = 75,    # Library complexity (did we sequence enough?)
       maxNTCCount = 9000,        # Max background noise allowed
       minNegativeCount = 0,      # Min background allowed
       minNuclei = 20,            # Minimum cells per ROI (biological sufficiency)
       minArea = 0                # Minimum surface area
  )

# Apply QC flags
demoData <- setSegmentQCFlags(demoData, qcCutoffs = QC_params)
QCResults <- protocolData(demoData)[["QCFlags"]]
QCResults$QCStatus <- apply(QCResults, 1L, function(x) ifelse(sum(x) == 0L, "PASS", "WARNING"))
protocolData(demoData)[["QCFlags"]] <- QCResults

# Summarizing QC results
QC_Summary <- data.frame(Pass = colSums(!QCResults[, -ncol(QCResults)]), 
                         Warning = colSums(QCResults[, -ncol(QCResults)]))

print("QC SUMMARY:")
print(QC_Summary)

# METHOD: Exclusion
# PURPOSE: Remove any sample that failed QC to prevent skewing statistical results.
demoData <- demoData[, QCResults$QCStatus == "PASS"]
print(paste("Remaining Sample Count After QC:", dim(demoData)[2]))

# -------------------------------------------------------------------------
# 4. GENE FILTERING & NORMALIZATION (Q3)
# -------------------------------------------------------------------------
# Filter 1: Keep only Endogenous genes (ignore control probes)
demoData <- demoData[featureData(demoData)$CodeClass == "Endogenous", ]
# Filter 2: Drop "Ghost" genes (genes with 0 counts in all samples)
demoData <- demoData[rowSums(exprs(demoData)) > 0, ]

print(paste("Gene Count Entering Analysis:", dim(demoData)[1]))

# METHOD: Q3 (Upper Quartile) Normalization
# PURPOSE: GeoMx ROIs vary in size and cell count. Standard TPM/RPM is not ideal.
# Q3 normalization scales samples based on the 75th percentile of expression,
# which is robust against outliers and differences in ROI area.
demoData <- normalize(demoData, norm_method = "quant", desiredQuantile = .75, toElt = "q_norm")

# METHOD: Log2 Transformation
# PURPOSE: Stabilize variance. Gene expression follows a power law; log-transform makes
# the distribution more normal (Gaussian), which is required for Linear Models (LMM).
assayDataElement(demoData, "log_q") <- log2(assayDataElement(demoData, "q_norm") + 1)

# Check
print("Normalization Done. First 5 rows:")
print(head(assayDataElement(demoData, "log_q")[,1:5]))

# -------------------------------------------------------------------------
# 5. VISUALIZATION: UMAP (Figure 1B)
# -------------------------------------------------------------------------
# METHOD: UMAP (Uniform Manifold Approximation and Projection)
# PURPOSE: Dimensionality Reduction. Collapses 18,000 genes into 2D coordinates.
# It helps visualize how samples cluster. We expect samples to cluster by Treatment or Genotype.

# STEP 1: Double Filtering
# We focus ONLY on "Epidermis" and "PanCK" (Epithelial) cells to reduce heterogeneity.
panck_epidermis_index <- which(
  grepl("epidermis", pData(demoData)$genotype, ignore.case = TRUE) & 
    pData(demoData)$cell_type == "PanCK"
)

# Control: How many found?
print(paste("Found PanCK+Epidermis Samples:", length(panck_epidermis_index)))

if(length(panck_epidermis_index) > 0) {
  
  # Prepare Data and Matrix
  all_matris <- assayDataElement(demoData, "log_q")
  # Get only selected indices
  target_matris <- all_matris[, panck_epidermis_index]
  
  # Calculate UMAP
  umap_out <- umap(t(target_matris))
  
  # Prepare Plot Data
  plot_data <- pData(demoData)[panck_epidermis_index, ]
  plot_data$UMAP1 <- umap_out$layout[, 1]
  plot_data$UMAP2 <- umap_out$layout[, 2]
  
  # PLOTTING
  # Color = Donor (Batch), Shape = Class (Treatment)
  # Interpretation: If points group by Color, we have a Batch Effect.
  # If they group by Shape, the Treatment is working.
  print(
    ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = batch, shape = treatment)) +
      geom_point(size = 5, alpha = 0.9) +
      theme_bw() +
      labs(title = "Figure 1B: UMAP of PanCK Epidermis Samples",
           x = "UMAP1", y = "UMAP2",
           color = "Donor", shape = "Class") +
      theme(panel.grid = element_line(color = "grey92"), 
            plot.title = element_text(hjust = 0.5, face = "bold"))
  )
} else {
  print("WARNING: No samples found matching PanCK and Epidermis criteria. Check 'cell_type' column.")
}

# -------------------------------------------------------------------------
# 6. DIFFERENTIAL EXPRESSION ANALYSIS - LMM (Linear Mixed Model)
# -------------------------------------------------------------------------
# METHOD: Linear Mixed Model (LMM) via 'lmerTest'
# PURPOSE: Unlike a standard T-test, GeoMx data has nested structures. 
# We often have multiple ROIs from the SAME mouse.
# Fixed Effect: Treatment (What we want to test: 0h vs 24h)
# Random Effect: (1 | batch) -> Accounts for correlation within the same mouse/donor.
run_DE_analysis_LMM <- function(object, cell_type_subset, p_thresh = 0.1, fc_thresh = 0.7) {
  
  message(paste("LMM Analysis Starting (This process may take 1-3 mins depending on gene count):", cell_type_subset))
  
  # 1. DATA FILTERING
  # Epidermis only, relevant cell type and 0h/24h timepoints
  keep_indices <- which(pData(object)$cell_type == cell_type_subset & 
                          pData(object)$treatment %in% c("0 h", "24 h") &
                          grepl("epidermis", pData(object)$genotype, ignore.case = TRUE))
  
  sub_obj <- object[, keep_indices]
  
  if(ncol(sub_obj) < 2) {
    warning("Insufficient sample size.")
    return(NULL)
  }
  
  # Set factors (Reference: 0 h) - Essential for calculating Fold Change direction
  pData(sub_obj)$treatment <- droplevels(factor(pData(sub_obj)$treatment, levels = c("0 h", "24 h")))
  
  # Get Data Matrix and Metadata
  expr_mat <- assayDataElement(sub_obj, "log_q")
  metadata <- pData(sub_obj)
  
  # 2. LMM LOOP (Gene by gene iteration)
  # Model Formula: Expression ~ Treatment + (1 | MouseID)
  
  results_list <- lapply(rownames(expr_mat), function(gene) {
    tryCatch({
      # Create temporary data table
      tmp_data <- metadata
      tmp_data$expression <- expr_mat[gene, ]
      
      # FIT THE MODEL
      # 'suppressMessages' keeps the console clean from convergence warnings
      fit <- suppressMessages(suppressWarnings(lmer(expression ~ treatment + (1 | batch), data = tmp_data)))
      
      # Extract coefficients
      coefs <- summary(fit)$coefficients
      
      # Find "treatment24 h" row (The comparison vs 0h)
      target_row <- grep("24 h", rownames(coefs))
      
      if(length(target_row) > 0) {
        logFC <- coefs[target_row, "Estimate"]    # The magnitude of change
        P.Value <- coefs[target_row, "Pr(>|t|)"]  # Statistical significance
        return(c(logFC = logFC, P.Value = P.Value))
      } else {
        return(c(logFC = NA, P.Value = NA))
      }
    }, error = function(e) {
      
      return(c(logFC = NA, P.Value = NA))
    })
  })
  
  # 3. MERGING RESULTS
  results <- do.call(rbind, results_list)
  results <- as.data.frame(results)
  rownames(results) <- rownames(expr_mat)
  
  # Clean NA rows (failed genes)
  results <- na.omit(results)
  
  # P-value Correction (FDR - Benjamini-Hochberg)
  # PURPOSE: Correct for multiple testing errors (testing 18,000 genes increases false positives).
  results$adj.P.Val <- p.adjust(results$P.Value, method = "BH")
  
  # Matching Gene Names (RTS IDs -> Readable Gene Name)
  rts_ids <- rownames(results)
  gene_names <- featureData(object)$TargetName[match(rts_ids, rownames(featureData(object)))]
  gene_names[is.na(gene_names)] <- rts_ids[is.na(gene_names)]
  results$Gene <- gene_names
  
  # Color Status (Paper Thresholds)
  # Categorizing genes into Up-regulated, Down-regulated, or Non-significant
  results$color_status <- "no"
  results$color_status[results$adj.P.Val < p_thresh & results$logFC > fc_thresh] <- "up"
  results$color_status[results$adj.P.Val < p_thresh & results$logFC < -fc_thresh] <- "dn"
  results$color_status <- factor(results$color_status, levels = c("dn", "no", "up"))
  
  return(results)
}

# --- Volcano Plot Function ---
# PURPOSE: Visualizing DEGs. 
# X-axis: Biological Impact (Fold Change). Y-axis: Statistical Confidence (-log10 P-value).
plot_volcano <- function(results, title_text, highlight_genes = NULL) {
  paper_colors <- c("dn" = "#2b59a8", "no" = "grey70", "up" = "#d61c1c")
  
  if(!is.null(highlight_genes)) {
    label_data <- subset(results, Gene %in% highlight_genes)
  } else {
    # If no list, label top 20 significant genes with lowest P-value
    label_data <- head(subset(results, color_status != "no"), 20)
  }
  
  ggplot(results, aes(x = logFC, y = -log10(P.Value), color = color_status)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = paper_colors) +
    geom_vline(xintercept = c(-0.7, 0.7), lty = 2, col = "black", alpha = 0.5) + 
    geom_hline(yintercept = -log10(0.1), lty = 2, col = "black", alpha = 0.5) + 
    
    # ggrepel ensures text labels do not overlap with each other
    geom_label_repel(data = label_data, 
                     aes(label = Gene), 
                     size = 3.5, 
                     color = "black",        
                     fill = "white",     
                     box.padding = 0.5,    
                     point.padding = 0.5,    
                     min.segment.length = 0,
                     max.overlaps = Inf,     
                     force = 7,              
                     show.legend = FALSE) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
    labs(title = title_text, 
         x = "Log2 Fold Change", 
         y = "-Log10 Adjusted P-value")
}

# --- Run Analysis ---
print("--- LMM Differential Expression Analysis Starting ---")

# Calling LMM function for PanCK (Epithelial) and CD45 (Immune) cells
res_panck <- run_DE_analysis_LMM(demoData, "PanCK", p_thresh = 0.1, fc_thresh = 0.7)
res_cd45 <- run_DE_analysis_LMM(demoData, "CD45", p_thresh = 0.1, fc_thresh = 0.7)

# Highlighted Genes in Paper (Validating known markers)
panck_genes <- c("Usp19", "Adam25", "Krt18", "Phlda1", "Adam8", "Ifitm1", 
                 "Krt17", "Cdh1", "Lor", "Krt77") 

cd45_genes <- c("Cxcl2", "Cxcl3", "Ccl4", "Il1b", "S100a9", "Ifitm1", 
                "Apoe", "Ptma", "Cd14", "Cst6")

# Plotting Graphs
print("--- Plotting Figure 1C & 1D ---")
if(!is.null(res_panck)) print(plot_volcano(res_panck, "Figure 1C: PanCK (LMM: 24h vs 0h)", panck_genes))
if(!is.null(res_cd45)) print(plot_volcano(res_cd45, "Figure 1D: CD45 (LMM: 24h vs 0h)", cd45_genes))

# -------------------------------------------------------------------------
# 7. VIOLIN PLOTS (Figure 1E & 1F)
# -------------------------------------------------------------------------
# PURPOSE: Visualize distribution of expression for specific candidate genes.
target_genes <- c("Krt18", "Phlda1")
target_rts <- rownames(demoData)[match(target_genes, featureData(demoData)$TargetName)]

if(length(target_rts) > 0) {
  
  # Fetch Data and reshape to "Long Format" for ggplot
  expr_data <- assayDataElement(demoData, "q_norm")[target_rts, , drop=FALSE]
  df <- as.data.frame(t(expr_data))
  colnames(df) <- target_genes
  df$sample_id <- rownames(df)
  
  df_meta <- pData(demoData)
  df_meta$sample_id <- rownames(df_meta)
  
  # Merge Expression + Metadata
  plot_df <- merge(df, df_meta, by="sample_id")
  # Melt: Converts wide table (genes as columns) to long table (1 row per gene per sample)
  plot_long <- melt(plot_df, id.vars = colnames(df_meta), measure.vars = target_genes, 
                    variable.name="Gene", value.name="Expression")
  
  # === CRITICAL FILTERS ===
  # To replicate paper figures, we must isolate the exact tissue context.
  plot_long <- plot_long %>% 
    # 1. Remove Adipose
    filter(!grepl("adipose", genotype, ignore.case = TRUE)) %>%
    # 2. SELECT PanCK CELLS ONLY (To match paper results)
    filter(cell_type == "PanCK") 
  
  # Genotype Name Cleaning for cleaner plot labels
  # "PanCK epidermis" -> "Epidermis"
  plot_long$genotype <- gsub("PanCK |CD45 ", "", plot_long$genotype)
  
  # Factor Ordering to control the X-axis order
  plot_long$genotype <- factor(plot_long$genotype, levels = c("deep dermis", "epidermis"))
  
  # NA Cleaning
  plot_long <- plot_long[!is.na(plot_long$genotype), ]
  
  # Log Transformation for visualization
  plot_long$LogExpression <- log2(plot_long$Expression + 1)
  
  # Time Ordering (0h -> 6h -> 24h)
  avail_times <- intersect(c("0 h", "6 h", "24 h"), unique(plot_long$treatment))
  plot_long$treatment <- factor(plot_long$treatment, levels = avail_times)
  
  # Plotting Function
  plot_violin <- function(data, gene, title) {
    # Colors in image: 0h(Blue), 6h(Yellow), 24h(Grey)
    fills <- c("0 h" = "#1f3b75", "6 h" = "#f7e64f", "24 h" = "#808080")
    
    ggplot(data, aes(x=treatment, y=LogExpression, fill=treatment)) +
      # Violin: Shows probability density of the data
      geom_violin(trim=FALSE, scale="width", color="black", alpha=0.9) + 
      # Jitter: Shows individual data points to reveal sample size
      geom_jitter(width=0.15, size=1.5, color="black", alpha=0.7) +
      
      scale_fill_manual(values=fills) + 
      facet_wrap(~genotype, scales="fixed") +
      
      theme_bw() + 
      theme(strip.background = element_rect(fill = "grey90"),
            strip.text = element_text(face="bold", size=12),
            plot.title = element_text(hjust = 0.5, face="bold", size=14),
            legend.position = "right") + # Legend on right
      labs(title=title, y="Expression", x="Class")
  }
  
  print("--- Plotting Figure 1E & 1F ---")
  print(plot_violin(subset(plot_long, Gene=="Krt18"), "Figure 1E: Krt18 (PanCK)", "Figure 1E"))
  print(plot_violin(subset(plot_long, Gene=="Phlda1"), "Figure 1F: Phlda1 (PanCK)", "Figure 1F"))
}