# Dynamic Execution of Level 1 Tasks

## Overview

This document demonstrates the dynamic execution framework for Level 1 (low complexity) bioinformatics tasks using the OmixBenchR package. Level 1 tasks represent the most straightforward category in our three-tier complexity classification system, characterized by:

- Simple analytical workflows (typically 1-3 computational steps)
- Basic statistical methods (descriptive statistics, simple hypothesis tests, or standard single-model analyses)
- Minimal data integration (usually single data source or simple two-file merging)
- Low computational intensity (no extensive iteration, cross-validation, or nested loops)

For detailed classification criteria, see the [Task Complexity Classification System](https://github.com/SolvingLab/OmixBench/blob/main/Task%20Complexity%20Classification.md).

These tasks serve as foundational benchmarks to assess whether large language models can handle routine bioinformatics operations that require domain knowledge but limited procedural complexity. While conceptually straightforward, successful execution still demands correct package usage, appropriate statistical methodology, and proper data handling—skills essential for any AI-assisted bioinformatics workflow.

## Framework Components

The execution framework utilizes two core R packages:

- **[OmixBenchR](https://github.com/SolvingLab/OmixBench/tree/main/OmixBenchR)**: Core benchmarking framework providing task execution, automated error correction, and result validation
- **[llmhelper](https://github.com/SolvingLab/OmixBench/tree/main/llmhelper)**: R interface for seamless large language model integration, enabling programmatic LLM-assisted code generation

The framework implements automatic error detection and iterative correction (up to 10 attempts per task), mirroring real-world interactive debugging workflows while maintaining reproducibility through containerized execution environments.

## Representative Level 1 Tasks

The examples below span multiple omics domains (genomics, transcriptomics, epigenomics, metabolomics, microbiome) to demonstrate the breadth of Level 1 task coverage.

### Example 1: Transcriptomics - Clustering Association

**Task Description:** Perform hierarchical clustering on top 500 most variable miRNAs and test association between resulting clusters and patient gender using Chi-squared test.

**Complexity Profile:**
- Number of data files: 2 (expression matrix and clinical data)
- Preprocessing complexity: Moderate (variance-based filtering and data merging)
- Analysis level: Standard (hierarchical clustering)
- Computational steps: 3 (filter features, perform clustering, statistical testing)
- Output complexity: Simple (single p-value)

### Example 2: Epigenomics - Methylation Profiling

**Task Description:** Calculate the highest chromosome 17 methylation level in GSE149282 dataset using ChAMP analysis pipeline and EPIC array annotation.

**Complexity Profile:**
- Number of data files: 1 (GEO series matrix)
- Preprocessing complexity: Moderate (ChAMP filtering and normalization)
- Analysis level: Standard (methylation profiling)
- Computational steps: 3 (data loading, normalization, calculation)
- Output complexity: Simple (single numeric value)

### Example 3: Microbiome - Community Dissimilarity

**Task Description:** Perform community dissimilarity analysis between habitat types using ANOSIM statistical test on GlobalPatterns microbial dataset.

**Complexity Profile:**
- Number of data files: 1 (built-in phyloseq dataset)
- Preprocessing complexity: Simple (direct distance calculation)
- Analysis level: Standard (ANOSIM)
- Computational steps: 2 (distance computation, statistical testing)
- Output complexity: Simple (single R statistic)

### Example 4: Metabolomics - Feature Importance

**Task Description:** Identify top discriminatory metabolite using PLS-DA analysis with Variable Importance in Projection (VIP) scoring on nutrimouse metabolomics dataset.

**Complexity Profile:**
- Number of data files: 1 (built-in dataset)
- Preprocessing complexity: Simple (direct extraction)
- Analysis level: Standard (PLS-DA)
- Computational steps: 2 (model fitting, VIP score extraction)
- Output complexity: Simple (single metabolite identifier)


## Installation

```r
# Install required packages
# devtools::install_github("SolvingLab/OmixBench/OmixBenchR")
# devtools::install_github("SolvingLab/OmixBench/llmhelper")
```

## Execute Level 1 Tasks

```r
rm(list = ls())

# Load libraries
library(jsonlite)
library(OmixBenchR)
library(llmhelper)

# Load example tasks
load('task_for_examples.rda')
current_dir <- getwd()

# Configure LLM client
llm_client <- llm_provider(
  base_url = Sys.getenv('DS_BASE_URL'),
  api_key = Sys.getenv('DS_API_KEY'),
  model = 'deepseek-chat',
  stream = TRUE
)

# Task display function
print_task_banner <- function(omics, task_prompt) {
  prompt_preview <- substr(task_prompt, 1, 100)
  if (nchar(task_prompt) > 100) {
    prompt_preview <- paste0(prompt_preview, "...")
  }
  
  separator <- paste(rep("=", 80), collapse = "")
  cat("\n", separator, "\n", sep = "")
  cat(sprintf("  Omics: %s\n", omics))
  cat(sprintf("  Task: %s\n", prompt_preview))
  cat(separator, "\n")
}

# Execute Level 1 tasks
resultL <- list()
for (idx in seq_along(level1_tasks)) {
  task_dir <- file.path(current_dir, 'OmixQA_Tasks', level1_tasks[idx])
  task_name <- level1_tasks[idx]
  
  tryCatch({
    # Change to task directory
    old_wd <- current_dir
    setwd(task_dir)
    
    # Load task metadata
    meta <- jsonlite::fromJSON('task_prompt.json')
    
    # Display task info
    print_task_banner(
      omics = meta$Omics,
      task_prompt = meta$Task_prompt
    )
    
    # Execute task with automatic error correction
    result <- Execute_Task(
      task_prompt = meta$Task_prompt,
      llm_client = llm_client,
      timeout_sec = 1200,
      max_interactions = 10
    )
    
    # Restore directory
    setwd(old_wd)
    
    # Store result
    resultL[[task_name]] <- result
    cat("✓ Task completed successfully\n")
    
  }, error = function(e) {
    cat("✗ Task failed:", task_name, "\n")
    cat("Error:", e$message, "\n")
    
    resultL[[task_name]] <- list(
      status = "failed",
      error = e$message,
      timestamp = Sys.time()
    )
    
    if (exists("old_wd")) setwd(old_wd)
  })
  
  cat("\n")
}
```

```
# ================================================================================
#   Omics: Transcriptomics
# Task: Task: Perform hierarchical clustering on top 500 most variable miRNAs and test association between r...
# ================================================================================ 
#   ⚡ Executing task: task_20251027_144233 
# 🔧 Execution engine: tidyprompt 
# --- Sending request to LLM provider (deepseek-chat): ---
#   Task: Perform hierarchical clustering on top 500 most variable miRNAs and test association between resulting clusters and patient gender using Chi-squared test.
# 
# Data: File tcga_luad_mirna_rpm_matrix.csv contains miRNA RPM expression matrix (first column: sample IDs, subsequent columns: miRNAs). File tcga_luad_clin.csv contains clinical data with 'sampleid' column and gender.
# 
# Package: readr, dplyr, stats
# 
# Requirements:
#   1. Load datasets, convert miRNA expression to numeric, remove NAs, identify top 500 most variable miRNAs by standard deviation
# 2. Perform hierarchical clustering using Euclidean distance and complete linkage, cut tree into 2 clusters, merge cluster assignments with gender data
# 3. Perform Chi-squared test between cluster assignments and gender using chisq.test(), extract p-value, round to 4 decimal places
# 
# **IMPORTANT**: Return result as final expression (don't use cat() or print())
# Output: Single numeric p-value from Chi-squared test for cluster-gender association (e.g., 0.0000)
# 
# IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).
# 
# You must code in the programming language 'R' to answer this prompt.
# You may not install or load any additional packages.
# --- Receiving response from LLM provider: ---
# ```r
# # Load required packages
# library(readr)
# library(dplyr)
# library(stats)
# 
# # Load miRNA expression data
# mirna_data <- readr::read_csv("tcga_luad_mirna_rpm_matrix.csv")
# 
# # Load clinical data
# clin_data <- readr::read_csv("tcga_luad_clin.csv")
# 
# # Convert miRNA expression to numeric and remove NAs
# mirna_numeric <- mirna_data %>%
#   dplyr::select(-1) %>%  # Remove sample ID column
#   dplyr::mutate_all(as.numeric) %>%
#   na.omit()
# 
# # Calculate standard deviation for each miRNA and select top 500 most variable
# variances <- apply(mirna_numeric, 2, stats::sd, na.rm = TRUE)
# top_500_indices <- order(variances, decreasing = TRUE)[1:500]
# top_500_mirna <- mirna_numeric[, top_500_indices]
# 
# # Perform hierarchical clustering
# distance_matrix <- stats::dist(t(top_500_mirna), method = "euclidean")
# hc_result <- stats::hclust(distance_matrix, method = "complete")
# cluster_assignments <- stats::cutree(hc_result, k = 2)
# 
# # Prepare cluster data for merging
# cluster_df <- data.frame(
#   miRNA = names(cluster_assignments),
#   cluster = as.factor(cluster_assignments)
# )
# 
# # Since we clustered miRNAs (columns), not samples, we need to transpose approach
# # Actually, let me reconsider - we need to cluster SAMPLES based on top 500 miRNAs
# 
# # Transpose the data to have samples as rows and miRNAs as columns
# sample_data <- t(mirna_numeric[, top_500_indices])
# colnames(sample_data) <- mirna_data[[1]]  # Set sample IDs as column names
# 
# # Perform hierarchical clustering on samples
# sample_distance <- stats::dist(t(sample_data), method = "euclidean")
# sample_hc <- stats::hclust(sample_distance, method = "complete")
# sample_clusters <- stats::cutree(sample_hc, k = 2)
# 
# # Create data frame with sample clusters
# sample_cluster_df <- data.frame(
#   sampleid = names(sample_clusters),
#   cluster = as.factor(sample_clusters)
# )
# 
# # Merge with clinical data
# merged_data <- dplyr::inner_join(sample_cluster_df, clin_data, by = "sampleid")
# 
# # Perform Chi-squared test
# chi_test <- stats::chisq.test(merged_data$cluster, merged_data$gender)
# 
# # Extract and round p-value
# round(chi_test$p.value, 4)
# ```
# [1] "New names:\n• `` -> `...1`\nRows: 538 Columns: 1001\n── Column specification ────────────────────────────────────────────────────────\nDelimiter: \",\"\nchr    (1): ...1\ndbl (1000): hsa-mir-21, hsa-mir-143, hsa-mir-375, hsa-mir-148a, hsa-mir-99b,...\n\nℹ Use `spec()` to retrieve the full column specification for this data.\nℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.\nRows: 538 Columns: 24\n── Column specification ────────────────────────────────────────────────────────\nDelimiter: \",\"\nchr (14): sampleid, barcode, sample_type, tissue_type, ajcc_pathologic_stage...\ndbl (10): age, cigarettes_per_day, os_event, dss_event, dfs_event, pfs_event...\n\nℹ Use `spec()` to retrieve the full column specification for this data.\nℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.\n"
# [1] "0.779"
# ✅ Task completed successfully!
# ⏱️  Duration: 26.38 seconds
# 🔄 Interactions: 2 
# Task completed successfully.
# 
# 
# ================================================================================
#   Omics: Epigenomics
#   Task: Task: Calculate the highest chromosome 17 mean methylation level per sample in GSE149282 EPIC array ...
# ================================================================================ 
# ⚡ Executing task: task_20251027_144300 
# 🔧 Execution engine: tidyprompt 
# --- Sending request to LLM provider (deepseek-chat): ---
# Task: Calculate the highest chromosome 17 mean methylation level per sample in GSE149282 EPIC array dataset
# 
# Data: Load GEO methylation dataset GSE149282_series_matrix.txt.gz from local file using GEOquery package (EPIC methylation array beta values with genomic annotations)
# 
# Processing Steps:
# 1. Load GEO dataset with getGPL=FALSE and extract beta values matrix
# 2. Load EPIC array annotation and filter for chromosome 17 probes
# 3. Subset beta values for chromosome 17 probes
# 4. Calculate mean methylation level for each sample across all chromosome 17 probes
# 5. Identify the maximum mean methylation value across all samples
# 6. Return value formatted to 5 decimal places
# 
# Output: Single numeric value with 5 decimal places (e.g., '0.85432')
# 
# Packages: GEOquery, Biobase, IlluminaHumanMethylationEPICanno.ilm10b4.hg19
# 
# Critical: Use getGPL=FALSE in getGEO(); use Biobase::exprs() to extract beta values; do NOT use ChAMP functions (champ.filter, champ.norm) as they require unavailable manifest data; load annotation using data(Locations, package='IlluminaHumanMethylationEPICanno.ilm10b4.hg19') or library then getAnnotation(); filter for chr=='chr17'; subset beta matrix using probe names; calculate colMeans with na.rm=TRUE; use sprintf('%.5f', max(...)) for formatting; return as final expression without cat() or print()
# 
# IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).
# 
# You must code in the programming language 'R' to answer this prompt.
# You may not install or load any additional packages.
# --- Receiving response from LLM provider: ---
# ```r
# # Load required packages
# library(GEOquery)
# library(Biobase)
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# 
# # Load GEO dataset
# gse <- GEOquery::getGEO(filename = "GSE149282_series_matrix.txt.gz", getGPL = FALSE)
# 
# # Extract beta values matrix
# beta_matrix <- Biobase::exprs(gse)
# 
# # Load EPIC array annotation and filter for chromosome 17 probes
# data(Locations, package = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
# chr17_probes <- rownames(Locations)[Locations$chr == "chr17"]
# 
# # Subset beta values for chromosome 17 probes
# chr17_beta <- beta_matrix[rownames(beta_matrix) %in% chr17_probes, ]
# 
# # Calculate mean methylation level for each sample across all chromosome 17 probes
# sample_means <- colMeans(chr17_beta, na.rm = TRUE)
# 
# # Identify the maximum mean methylation value across all samples and format
# sprintf('%.5f', max(sample_means, na.rm = TRUE))
# ```
# [1] "0.56605"
# ✅ Task completed successfully!
# ⏱️  Duration: 28.32 seconds
# 🔄 Interactions: 2 
# Task completed successfully.
# 
# 
# ================================================================================
#   Omics: Microbiome
#   Task: Task: Perform community dissimilarity analysis between habitat types using ANOSIM statistical test o...
# ================================================================================ 
# ⚡ Executing task: task_20251027_144328 
# 🔧 Execution engine: tidyprompt 
# --- Sending request to LLM provider (deepseek-chat): ---
# Task: Perform community dissimilarity analysis between habitat types using ANOSIM statistical test on GlobalPatterns microbial dataset.
# 
# Data: Load built-in GlobalPatterns dataset from phyloseq package containing microbiome data with OTU abundance table, sample metadata, and taxonomic information. Dataset includes sample information accessible via sample_data function with SampleType column indicating habitat types for group comparison, OTU abundance counts for dissimilarity calculations, and taxonomic classifications. Data structure represents microbial community composition across multiple environmental sample types with count data for analysis of similarities using Bray-Curtis dissimilarity between habitat groups.
# 
# Package: phyloseq, vegan
# 
# Requirements:
# 1. Load GlobalPatterns dataset using data function, calculate Bray-Curtis distance matrix using distance function with method='bray' parameter on phyloseq object
# 2. Extract SampleType grouping variable from sample metadata using sample_data function and convert to factor using as.factor function
# 3. Apply ANOSIM test using anosim function with Bray-Curtis distance matrix and grouping factor to test differences between habitat types
# 4. Extract R statistic from ANOSIM results using $statistic component, round to 3 decimal places using round function
# 5. Return R value representing degree of separation between groups where 0 indicates no separation and 1 indicates complete separation
# 
# **IMPORTANT**: Return result as final expression (don't use cat() or print())
# Output: Single numeric R statistic value rounded to 3 decimal places (e.g., '0.765' or '0.234')
# 
# IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).
# 
# You must code in the programming language 'R' to answer this prompt.
# You may not install or load any additional packages.
# --- Receiving response from LLM provider: ---
#   ```r
# library(phyloseq)
# library(vegan)
# 
# data(GlobalPatterns)
# 
# bray_dist <- phyloseq::distance(GlobalPatterns, method = "bray")
# sample_types <- as.factor(phyloseq::sample_data(GlobalPatterns)$SampleType)
# anosim_result <- vegan::anosim(bray_dist, sample_types)
# round(anosim_result$statistic, 3)
# ```
# [1] "0.99"
# ✅ Task completed successfully!
#   ⏱️  Duration: 6.98 seconds
# 🔄 Interactions: 2 
# Task completed successfully.
# 
# 
# ================================================================================
#   Omics: Microbiome
# Task: Task: Calculate phylogenetic distance between specific samples using weighted UniFrac metric on esop...
# ================================================================================ 
#   ⚡ Executing task: task_20251027_144335 
# 🔧 Execution engine: tidyprompt 
# --- Sending request to LLM provider (deepseek-chat): ---
#   Task: Calculate phylogenetic distance between specific samples using weighted UniFrac metric on esophagus microbial dataset.
# 
# Data: Load built-in esophagus dataset from phyloseq package containing microbial community data with OTU abundance table, sample metadata, taxonomic information, and phylogenetic tree. Dataset includes sample identifiers accessible for pairwise distance calculations, OTU abundance counts across samples, taxonomic classifications, and phylogenetic relationships required for UniFrac calculations. Data structure represents microbial community composition with phylogenetic information enabling evolutionary distance-based community comparisons between samples.
# 
# Package: phyloseq
# 
# Requirements:
#   1. Load esophagus dataset using data function, compute weighted UniFrac distance matrix using phyloseq::distance function with method='unifrac' and weighted=TRUE parameters
# 2. Convert distance object to matrix format using as.matrix function to enable indexing for specific sample pairs
# 3. Extract pairwise distance between samples 'B' and 'C' using matrix indexing with sample identifiers
# 4. Return distance value representing phylogenetic dissimilarity incorporating both abundance and evolutionary divergence
# 
# **IMPORTANT**: Return result as final expression (don't use cat() or print())
# Output: Single numeric distance value representing weighted UniFrac distance with 6 decimal places (e.g., '0.123456' or '0.876543')
# 
# IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).
# 
# You must code in the programming language 'R' to answer this prompt.
# You may not install or load any additional packages.
# --- Receiving response from LLM provider: ---
# ```r
# data(esophagus, package = "phyloseq")
# dist_matrix <- phyloseq::distance(esophagus, method = "unifrac", weighted = TRUE)
# dist_matrix <- as.matrix(dist_matrix)
# dist_matrix["B", "C"]
# ```
# [1] "0.203542384088409"
# ✅ Task completed successfully!
# ⏱️  Duration: 5.55 seconds
# 🔄 Interactions: 2 
# Task completed successfully.
# 
# 
# ================================================================================
#   Omics: Metabolomics
#   Task: Task: Identify top discriminatory metabolite using PLS-DA analysis with Variable Importance in Proje...
# ================================================================================ 
# ⚡ Executing task: task_20251027_144341 
# 🔧 Execution engine: tidyprompt 
# --- Sending request to LLM provider (deepseek-chat): ---
# Task: Identify top discriminatory metabolite using PLS-DA analysis with Variable Importance in Projection scoring on nutrimouse metabolomics dataset.
# 
# Data: Load built-in nutrimouse dataset from mixOmics package containing nutrigenomic study data with lipid metabolite abundance matrix and genotype classification labels. Dataset includes nutrimouse$lipid component with metabolite abundance measurements across samples and nutrimouse$genotype component with group labels for supervised classification analysis. Data structure represents targeted lipidomics profiles requiring Partial Least Squares Discriminant Analysis for feature importance identification and metabolite ranking.
# 
# Package: mixOmics
# 
# Requirements:
# 1. Load nutrimouse dataset using data function from mixOmics package, extract metabolite abundance matrix using nutrimouse$lipid and group labels using nutrimouse$genotype for supervised analysis
# 2. Perform Partial Least Squares Discriminant Analysis using plsda function with parameters: X=metabolite matrix, Y=group labels, ncomp=2 for two-component model
# 3. Calculate Variable Importance in Projection scores using vip function on PLS-DA result object to quantify metabolite contributions to group discrimination
# 4. Extract VIP scores for first component using column indexing [,1] on VIP results matrix, identify metabolite with highest VIP score using which.max function
# 5. Extract metabolite name using names function on VIP scores vector with maximum VIP index to obtain top discriminatory metabolite identifier
# 
# **IMPORTANT**: Return result as final expression (don't use cat() or print())
# Output: Single metabolite name string representing highest discriminatory metabolite (e.g., 'MetaboliteName' or 'CompoundID')
# 
# IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).
# 
# You must code in the programming language 'R' to answer this prompt.
# You may not install or load any additional packages.
# --- Receiving response from LLM provider: ---
#   ```r
# data("nutrimouse", package = "mixOmics")
# metabolite_matrix <- nutrimouse$lipid
# group_labels <- nutrimouse$genotype
# plsda_result <- mixOmics::plsda(X = metabolite_matrix, Y = group_labels, ncomp = 2)
# vip_scores <- mixOmics::vip(plsda_result)
# top_metabolite <- names(which.max(vip_scores[,1]))
# top_metabolite
# ```
# [1] "C16.0"
# ✅ Task completed successfully!
#   ⏱️  Duration: 6.66 seconds
# 🔄 Interactions: 2 
# Task completed successfully.
# 
# 
# ================================================================================
#   Omics: Genomics
# Task: Task: Analyze individual genetic diversity patterns and identify the most heterozygous individual in...
# ================================================================================ 
#   ⚡ Executing task: task_20251027_144348 
# 🔧 Execution engine: tidyprompt 
# --- Sending request to LLM provider (deepseek-chat): ---
#   Task: Analyze individual genetic diversity patterns and identify the most heterozygous individual in DROPS dataset.
# 
# Data: Built-in DROPS marker dataset from statgenGWAS package with 246 individuals and 41,722 SNPs. Genotype calls encoded as numeric values (0=homozygous reference, 1=heterozygous, 2=homozygous alternative). Individual genetic diversity measured by overall heterozygosity rate across genome.
# 
# Package: statgenGWAS
# 
# Requirements:
#   1. Load DROPS marker data using data function with 'dropsMarkers' from statgenGWAS package
# 2. Remove 'Ind' column and calculate individual heterozygosity rates using rowSums function: count SNPs with value=1 divided by total non-missing SNPs per individual
# 3. Create histogram of heterozygosity distribution across individuals using hist function
# 4. Calculate additional diversity metrics: proportion of rare genotypes (MAF < 0.1) per individual using apply function
# 5. Identify individual with highest overall heterozygosity rate using which.max function and return individual ID
# 
# **IMPORTANT**: Return result as final expression (don't use cat() or print())
# Output: Single individual identifier representing most genetically diverse individual (e.g., 'A310')
# 
# IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).
# 
# You must code in the programming language 'R' to answer this prompt.
# You may not install or load any additional packages.
# --- Receiving response from LLM provider: ---
# ```r
# # Load the DROPS marker data from statgenGWAS package
# data("dropsMarkers", package = "statgenGWAS")
# 
# # Remove the 'Ind' column and work with the genotype matrix
# geno_data <- dropsMarkers[, -1]
# 
# # Calculate heterozygosity rate for each individual
# # Count heterozygous calls (value = 1) divided by total non-missing SNPs per row
# het_rates <- rowSums(geno_data == 1, na.rm = TRUE) / rowSums(!is.na(geno_data))
# 
# # Create histogram of heterozygosity distribution
# graphics::hist(het_rates, main = "Distribution of Individual Heterozygosity Rates", 
#                xlab = "Heterozygosity Rate", col = "lightblue")
# 
# # Calculate proportion of rare genotypes (MAF < 0.1) per individual
# # First calculate MAF for each SNP
# maf <- colSums(geno_data, na.rm = TRUE) / (2 * colSums(!is.na(geno_data)))
# maf <- pmin(maf, 1 - maf)  # Take the minor allele frequency
# 
# # Identify SNPs with MAF < 0.1
# rare_snps <- maf < 0.1
# 
# # Calculate proportion of rare genotypes per individual
# rare_geno_prop <- apply(geno_data[, rare_snps], 1, function(x) {
#   sum(!is.na(x)) / length(x)
# })
# 
# # Identify individual with highest heterozygosity rate
# max_het_ind <- dropsMarkers$Ind[which.max(het_rates)]
# 
# # Return the individual with highest heterozygosity as final expression
# max_het_ind
# ```
# [1] "B113"
# ✅ Task completed successfully!
# ⏱️  Duration: 16.99 seconds
# 🔄 Interactions: 2 
# Task completed successfully.
```

## Citation

```bibtex
@article{liu2025systematic,
  title={Systematic Evaluation and Strategic Optimization of Large Language Models for Multi-omics Analysis},
  author={Liu, Zaoqu and Wu, Yushuai and Yang, Jingkuan and others},
  year={2025}
}
```
