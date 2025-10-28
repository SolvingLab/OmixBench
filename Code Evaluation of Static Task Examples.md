# Code Evaluation of Static Task Examples

## Overview

This document demonstrates the complete workflow for generating and evaluating LLM responses for bioinformatics tasks. We showcase the static evaluation framework using three state-of-the-art reasoning models: **DeepSeek-R1**, **OpenAI O3**, and **Claude Opus 4**.

## Evaluation Framework

Our evaluation employs a **7-dimension scoring system** (100 points total):

| Dimension | Points | Description |
|-----------|--------|-------------|
| Problem Solving | 20 | Biological appropriateness and analytical strategy |
| Technical Implementation | 25 | Code quality and computational rigor |
| Data Handling | 15 | Input/output management and preprocessing |
| Statistical Rigor | 15 | Statistical methods and validation |
| Domain Knowledge | 10 | Bioinformatics expertise and interpretation |
| Robustness | 10 | Error handling and edge cases |
| Documentation | 5 | Code clarity and usability |

### Evaluation Prompt

The complete evaluation prompt that defines the scoring criteria and assessment methodology is available at:

**[Bioinformatics Code Static Evaluation Prompt](https://github.com/SolvingLab/OmixBench/blob/main/OmixBenchR/inst/Bioinformatics%20Code%20Static%20Evaluation%20Prompt.md)**

---

## Installation

```r
# Install required packages
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install llmhelper for LLM integration
devtools::install_github('SolvingLab/OmixBench', subdir = 'llmhelper')

# Install OmixBench for evaluation
devtools::install_github('SolvingLab/OmixBench', subdir = 'OmixBenchR')
```

---

## Evaluation Example: TMB Analysis

### Task Definition

**Domain:** Genomics

**Task:** Identify tumor mutational burden (TMB) in a melanoma patient cohort by analyzing VCF files from whole exome sequencing, and generate a summary table of TMB values per patient.

### Step 1: Generate LLM Response

```r
library(llmhelper)
library(glue)

# Configure LLM client (example with o4-mini)
response_client <- llm_provider(
  model = 'o4-mini',
  temperature = 0.2,
  max_tokens = 20000,
  timeout = 500,
  stream = FALSE,
  verbose = FALSE
)

# Define task
task <- 'Identify tumor mutational burden (TMB) in a melanoma patient cohort by analyzing VCF files from whole exome sequencing, and generate a summary table of TMB values per patient.'

# Generate response
response <- get_llm_response(
  set_prompt(
    system = "You are an AI assistant specialized in bioinformatics. Provide complete R code solutions.",
    user = glue::glue("Please provide R code to solve the following bioinformatics task:\n\n{task}")
  ),
  response_client,
  max_retries = 1,
  verbose = TRUE
)
```

### Step 2: Configure Evaluator Models

```r
library(OmixBenchR)

# DeepSeek-R1 evaluator
llm_ds <- llm_provider(
  api_key = Sys.getenv('DeepSeek_API_KEY'),
  model = 'deepseek-reasoner',
  temperature = 0.2,
  timeout = 500,
  stream = TRUE,
  verbose = FALSE
)

# OpenAI O3 evaluator
llm_o3 <- llm_provider(
  api_key = Sys.getenv('O3_API_KEY'),
  model = 'o3',
  temperature = 0.2,
  timeout = 500,
  stream = TRUE,
  verbose = FALSE
)

# Claude Opus 4 evaluator
llm_opus4t <- llm_provider(
  api_key = Sys.getenv('OPUS4T_API_KEY'),
  model = 'claude-opus-4-20250514-thinking',
  temperature = 0.2,
  timeout = 500,
  stream = TRUE,
  verbose = FALSE
)
```

### Step 3: Evaluate Response

```r
# Evaluate with DeepSeek-R1
ds_res <- OmixBenchR::CodeEval(
  task, response, 
  llm_client = llm_ds,
  schema_type = 'text-based',
  verbose = FALSE
)

# Evaluate with OpenAI O3
o3_res <- OmixBenchR::CodeEval(
  task, response, 
  llm_client = llm_o3,
  schema_type = 'text-based',
  verbose = FALSE
)

# Evaluate with Claude Opus 4
opus4t_res <- OmixBenchR::CodeEval(
  task, response, 
  llm_client = llm_opus4t,
  schema_type = 'text-based',
  verbose = FALSE
)
```

### Step 4: Compare Results

```r
# Calculate final averaged score
final_score <- mean(c(
  ds_res$total_score, 
  o3_res$total_score, 
  opus4t_res$total_score
))

cat(sprintf("Final Score: %.2f/100\n", final_score))
```

---

## Evaluation Results

**Average Score: 85.33/100**

### Detailed Breakdown

#### DeepSeek-R1 Evaluation
```r
$total_score
[1] 92

$breakdown
  $problem_solving: 20/20
  $technical_implementation: 25/25
  $data_handling: 15/15
  $statistical_rigor: 15/15
  $domain_knowledge: 10/10
  $robustness: 5/10
  $documentation: 2/5
```

#### Claude Opus 4 Evaluation
```r
$total_score
[1] 85

$breakdown
  $problem_solving: 20/20
  $technical_implementation: 25/25
  $data_handling: 13/15
  $statistical_rigor: 15/15
  $domain_knowledge: 10/10
  $robustness: 2/10
  $documentation: 0/5
```

#### OpenAI O3 Evaluation
```r
$total_score
[1] 79

$breakdown
  $problem_solving: 18/20
  $technical_implementation: 21/25
  $data_handling: 11/15
  $statistical_rigor: 15/15
  $domain_knowledge: 8/10
  $robustness: 3/10
  $documentation: 3/5
```

---

## Generated Code Sample

<details>
<summary>View Complete TMB Analysis Code</summary>

```r
# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

pkgList <- c("VariantAnnotation", "GenomicFeatures", "GenomicRanges", "dplyr", "data.table")
for (pkg in pkgList) {
    if (!requireNamespace(pkg, quietly=TRUE))
        BiocManager::install(pkg, ask=FALSE)
}

library(VariantAnnotation)
library(GenomicFeatures)
library(GenomicRanges)
library(dplyr)
library(data.table)

# Load a TxDb for human coding regions
txdb <- makeTxDbFromUCSC(genome="hg19", tablename="knownGene")
exons_gr <- exons(txdb)

# Calculate total exome size in megabases
total_exome_bp <- sum(width(reduce(exons_gr)))
total_exome_mb <- total_exome_bp / 1e6
message("Total exome size (MB): ", round(total_exome_mb, 2))

# Point to your directory of VCFs
vcf_dir <- "path/to/vcf_directory"
vcf_files <- list.files(vcf_dir, pattern="\\.vcf(\\.gz)?$", full.names=TRUE)

results <- data.table(sample=character(), n_mutations=integer(), TMB=numeric())

# Loop through VCFs
for (vcf_file in vcf_files) {
    sample_name <- tools::file_path_sans_ext(basename(vcf_file))
    vcf <- readVcf(vcf_file, genome="hg19")
    vcf <- vcf[fixed(vcf)$FILTER == "PASS", ]
    is_snv <- width(ref(vcf)) == 1 & width(unlist(alt(vcf))) == 1
    vcf <- vcf[is_snv, ]
    vr <- rowRanges(vcf)
    hits <- findOverlaps(vr, exons_gr)
    vr_exonic <- vr[unique(queryHits(hits))]
    n_mut <- length(vr_exonic)
    tmb_val <- n_mut / total_exome_mb
    
    results <- rbind(results, data.table(
        sample = sample_name,
        n_mutations = n_mut,
        TMB = round(tmb_val, 2)
    ))
    
    message(sprintf("Processed %s: %d mutations, TMB=%.2f", 
                    sample_name, n_mut, tmb_val))
}

# Save summary table
out_file <- "TMB_summary_per_sample.csv"
fwrite(results, file=out_file)
message("Summary written to: ", out_file)
```

</details>

---


## Multi-Evaluator Workflow

### Complete Evaluation Pipeline

```r
# Define task list (10 omics domains)
task_examples <- c(
  "Identify tumor mutational burden (TMB) in a melanoma patient cohort...",
  "Perform WGCNA on RNA-seq data from a breast cancer patient cohort...",
  "Identify enriched transcription factor motifs in ATAC-seq peak regions..."
)

# Initialize results storage
evaluation_results <- list()

# Evaluate all tasks
for (i in seq_along(task_examples)) {
  task <- task_examples[i]
  
  # Generate response
  response <- get_llm_response(
    set_prompt(
      system = "You are an AI assistant specialized in bioinformatics. Provide complete R code solutions.",
      user = glue::glue("Please provide R code to solve the following bioinformatics task:\n\n{task}")
    ),
    response_client,
    max_retries = 1,
    verbose = FALSE
  )
  
  # Evaluate with three models
  ds_score <- CodeEval(task, response, llm_ds, schema_type = 'text-based', verbose = FALSE)
  o3_score <- CodeEval(task, response, llm_o3, schema_type = 'text-based', verbose = FALSE)
  opus4t_score <- CodeEval(task, response, llm_opus4t, schema_type = 'text-based', verbose = FALSE)
  
  # Store results
  evaluation_results[[i]] <- list(
    task = task,
    response = response,
    scores = list(
      deepseek = ds_score,
      o3 = o3_score,
      claude = opus4t_score,
      average = mean(c(ds_score$total_score, o3_score$total_score, opus4t_score$total_score))
    )
  )
  
  cat(sprintf("Task %d/%d completed. Average score: %.2f\n", i, length(task_examples), 
              evaluation_results[[i]]$scores$average))
}

# Summary statistics
avg_scores <- sapply(evaluation_results, function(x) x$scores$average)
cat(sprintf("\nOverall Performance:\n"))
cat(sprintf("  Mean Score: %.2f ± %.2f\n", mean(avg_scores), sd(avg_scores)))
cat(sprintf("  Min Score:  %.2f\n", min(avg_scores)))
cat(sprintf("  Max Score:  %.2f\n", max(avg_scores)))
```

---

## Visualization

### Score Distribution

```r
library(ggplot2)
library(tidyr)

# Prepare data
scores_df <- data.frame(
  Task = paste0("Task ", 1:10),
  DeepSeek = sapply(evaluation_results, function(x) x$scores$deepseek$total_score),
  O3 = sapply(evaluation_results, function(x) x$scores$o3$total_score),
  Claude = sapply(evaluation_results, function(x) x$scores$claude$total_score)
) %>%
  pivot_longer(cols = c(DeepSeek, O3, Claude), names_to = "Evaluator", values_to = "Score")

# Create plot
ggplot(scores_df, aes(x = Task, y = Score, fill = Evaluator)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 85, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Static Evaluation Scores Across 10 Omics Tasks",
       subtitle = "Comparison of three reasoning-enhanced evaluator models",
       x = "Task", y = "Score (0-100)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

---

## Resources

### Package Documentation

- **llmhelper**: https://github.com/SolvingLab/OmixBench/tree/main/llmhelper
- **OmixBench**: https://github.com/SolvingLab/OmixBench/tree/main/OmixBenchR
- **Evaluation Prompt**: https://github.com/SolvingLab/OmixBench/blob/main/OmixBenchR/inst/Bioinformatics%20Code%20Static%20Evaluation%20Prompt.md


### Citation

```bibtex
@article{liu2025systematic,
  title={Systematic Evaluation and Strategic Optimization of Large Language Models for Multi-omics Analysis},
  author={Liu, Zaoqu and Wu, Yushuai and Yang, Jingkuan and others},
  journal={In preparation},
  year={2025}
}
```

---

## Contact

For questions or technical support:
- Email: liuzaoqu@163.com
- GitHub Issues: https://github.com/SolvingLab/OmixBench/issues
