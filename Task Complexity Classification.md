# Task Complexity Classification System

## Overview

This document describes the multi-dimensional framework used to classify 405 bioinformatics tasks into three complexity levels (Level 1-3). The classification was developed by an expert panel of three doctoral-level researchers with complementary bioinformatics expertise, with all assignments validated through independent task-solving attempts.

---

## Six Evaluation Dimensions

Each task is systematically evaluated across six dimensions:

### 1. Number of Input Data Files (`N_files`)
**Range**: 1-5

Quantifies data integration requirements by counting distinct input files.

| Score | Description | Examples |
|-------|-------------|----------|
| 1 | Single data source | One expression matrix |
| 2 | Two files requiring merging | Expression + clinical data |
| 3 | Three files | Multi-omics pair + metadata |
| 4-5 | Complex multi-omics integration | 4+ data types requiring harmonization |

---

### 2. Data Preprocessing Complexity (`C_preprocess`)
**Range**: 0-2

Assesses sophistication of data manipulation before primary analysis.

| Score | Complexity | Operations |
|-------|------------|------------|
| 0 | **Simple** | Basic file loading, single column extraction, straightforward filtering |
| 1 | **Moderate** | Dataset merging by IDs, multi-step filtering, feature selection, binary encoding |
| 2 | **Complex** | Multi-omics integration (3+ types), imputation strategies, batch effect correction |

---

### 3. Advanced Analysis Level (`L_analysis`)
**Range**: 0-2

Evaluates statistical and computational sophistication of primary methods.

| Score | Level | Methods |
|-------|-------|---------|
| 0 | **Basic** | Descriptive statistics, simple hypothesis tests (t-test, correlation, ANOVA) |
| 1 | **Standard** | Regression models, dimensionality reduction, basic clustering (Cox, limma, LASSO, PCA, k-means) |
| 2 | **Advanced** | Network analysis, multi-step integrated pipelines, ensemble methods |

---

### 4. Number of Computational Steps (`N_steps`)
**Range**: 1-10 → **Transformed to `C_workflow` (0-2)**

Critical dimension reflecting cognitive load and potential for cascading errors.

**Step Counting Criteria**:
- Data loading/preprocessing with complex merging: 1 step
- Feature selection/filtering: 1 step
- Each distinct statistical model: 1 step each
- Multiple testing correction: 1 step
- Pathway/enrichment analysis: 1-2 steps
- Network construction/analysis: 1-2 steps
- Survival analysis: 1 step
- Complex result extraction: 1 step

**Transformation to Workflow Complexity**:

| `N_steps` | `C_workflow` | Category |
|-----------|--------------|----------|
| 1-3 | 0 | Simple workflows |
| 4-6 | 1 | Moderate workflows |
| 7+ | 2 | Complex workflows |

---

### 5. Output Complexity (`C_output`)
**Range**: 0-2

Characterizes structural complexity of expected outputs.

| Score | Type | Examples |
|-------|------|----------|
| 0 | **Simple** | Single numeric value, single string, simple count |
| 1 | **Moderate** | Two components with delimiter (e.g., `GENE:0.0234`) |
| 2 | **Complex** | Three+ components requiring interpretation (e.g., `GENE:HR:CI:PATHWAY`) |

---

### 6. Computational Intensity (`I_compute`)
**Range**: 0-2

Assesses overall computational demand based on algorithmic complexity.

| Score | Intensity | Characteristics |
|-------|-----------|-----------------|
| 0 | **Low** | Simple single-pass algorithms, no iteration (e.g., summary statistics) |
| 1 | **Medium** | Standard iteration patterns, per-gene testing, basic cross-validation |
| 2 | **High** | Extensive iteration, nested loops, combinatorial searches, hyperparameter tuning, iterative optimization |

---

## Composite Scoring Formula

The final complexity score combines all dimensions with differential weighting:

```
S_complexity = (N_files/5 × 10) + (C_preprocess/2 × 10) + (L_analysis/2 × 30) + 
               (C_workflow/2 × 30) + (C_output/2 × 10) + (I_compute/2 × 10)
```

**Weight Distribution**:
- 🔷 **30 points**: Advanced analysis level (`L_analysis`)
- 🔷 **30 points**: Workflow complexity (`C_workflow`)
- 🔹 **10 points**: Number of data files (`N_files`)
- 🔹 **10 points**: Preprocessing complexity (`C_preprocess`)
- 🔹 **10 points**: Output complexity (`C_output`)
- 🔹 **10 points**: Computational intensity (`I_compute`)

**Total Range**: 0-100 points

> **Rationale**: Analytical sophistication and workflow complexity receive higher weights as they most directly capture cognitive and computational demands.

---

## Complexity Level Assignment

Tasks are stratified into three levels using **tertile-based cutoffs** from the empirical distribution:

| Level | Percentile Range | `S_complexity` | Characteristics |
|-------|------------------|----------------|-----------------|
| **Level 1** (Low) | 0-33.33% | Lower tertile | Straightforward procedures, limited integration, basic statistics |
| **Level 2** (Moderate) | 33.33-66.67% | Middle tertile | Multi-step workflows, standard regression/clustering |
| **Level 3** (High) | 66.67-100% | Upper tertile | Advanced methods, multi-omics integration, iterative optimization |

This ensures **balanced representation** (~135 tasks per level) while maintaining discriminatory power.


## Methodological Considerations

### Limitations
- Complexity assessment involves expert judgment; alternative weighting schemes could yield different stratifications
- Dimension weights determined through iterative refinement and expert consensus (not formal optimization)
- Future work could explore ML-based complexity prediction or psychometric validation

### Reproducibility
All dimension ratings, composite scores, and level assignments are documented in our repository alongside complete task descriptions, enabling:
- Transparent evaluation of classification decisions
- Future refinement of the framework
- Replication and extension by other researchers

---

## Citation

```bibtex
@article{liu2025systematic,
  title={Systematic Evaluation and Strategic Optimization of Large Language Models for Multi-omics Analysis},
  author={Liu, Zaoqu and Wu, Yushuai and Yang, Jingkuan and others},
  year={2025}
}
```

---

For detailed examples, see task documentation in `OmixQA_Tasks/`.
