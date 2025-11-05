# AutoFlow - Step by Step Example

This document demonstrates how AutoFlow works internally by running each step separately.

---

## Setup

Load required libraries and prepare the LLM clients:

```r
# Run AutoFlow step by step
library(llmflow)
library(tidyverse)
library(jsonlite)

# Load task data
load("task.rda")

# Configure LLMs: separate models for RAG and ReAct
rag_llm <- chat_anthropic(model = 'claude-opus-4-20250514-thinking')  # For documentation retrieval
react_llm <- chat_deepseek(model = 'deepseek-chat')                   # For task execution
```

---

## Step 1: RAG - Retrieve Function Documentation

Use the RAG LLM to identify relevant functions and retrieve their documentation:

```r
# Retrieve documentation for functions needed by the task
ret_docs_prompt <- tryCatch(
  retrieve_docs(
    rag_llm,  # Use RAG model to identify relevant functions
    prompt = task_prompt,
    skip_undocumented = FALSE,
    use_llm_fallback = TRUE
  ),
  error = function(e) {
    if (verbose) warning("⚠️  ", e$message, call. = FALSE)
    ""
  }
)

#The LLM only outputs the JSON content as shown below:
#{
#  CellChat::createCellChat
#  CellChat::computeCommunProb
#  CellChat::subsetCommunication
#}
#Then, the retrieve_docs function will antomatically extract the docs of these functions.
```

### RAG Output - Function Identification

The RAG LLM analyzes the task and identifies critical functions:

```
# Looking at this cell-cell communication analysis task, I need to identify only the CellChat-specific functions that have complex 
# parameters or non-obvious usage patterns that would be difficult for an LLM to use correctly without documentation.
# 
# After careful analysis:
#   
#   {"functions": ["CellChat::createCellChat", "CellChat::computeCommunProb", "CellChat::subsetCommunication"]}
# 
# These three functions are truly critical because:
#   - `createCellChat`: Requires specific input structure and grouping variable setup
#   - `computeCommunProb`: Has the critical `population.size=FALSE` parameter that affects calculations
#   - `subsetCommunication`: Requires specific syntax for extracting directional cell-cell interactions (DC->NK)
# 
# All other operations (data loading, filtering NAs, basic subsetting) use standard R/Seurat functions that any LLM would handle correctly 
# without special documentation.

# Generating examples for CellChat::computeCommunProb ...
# Example 1: Basic computation without considering population size
# Typical for analyzing communication strength independent of cell abundance
cellchat <- computeCommunProb(cellchat, 
                              population.size = FALSE)  # Don't weight by cell numbers

# Example 2: Advanced computation with population weighting and custom parameters
# Useful when cell type abundance matters for biological interpretation
cellchat <- computeCommunProb(cellchat, 
                              type = "truncatedMean",   # Use truncated mean method
                              trim = 0.1,               # Trim 10% of extreme values
                              population.size = TRUE,   # Weight by cell population size
                              raw.use = FALSE)          # Use normalized data (default)
```

### View Retrieved Documentation

```r
# Print the complete documentation prompt
cat(ret_docs_prompt)
```

**Output:**

```
# Functions for reference
You may use the following functions for this task (**If the functions are not listed quite correctly, please ignore the functional reference below.**):
- CellChat::createCellChat
- CellChat::computeCommunProb
- CellChat::subsetCommunication

# Function Documentation & Examples

Below are detailed examples and usage patterns for each function. Use these to understand how to properly implement the analysis task:

## CellChat::createCellChat

## Not run: 
##D Create a CellChat object from single-cell transcriptomics data
##D # Input is a data matrix
##D ## create a dataframe consisting of the cell labels
##D meta = data.frame(labels = cell.labels, row.names = names(cell.labels))
##D cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
##D 
##D # input is a Seurat object
##D ## use the default cell identities of Seurat object
##D cellChat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
##D ## use other meta information as cell groups
##D cellChat <- createCellChat(object = seurat.obj, group.by = "seurat.clusters")
##D 
##D # input is a SingleCellExperiment object
##D cellChat <- createCellChat(object = sce.obj, group.by = "sce.clusters")
##D 
##D Create a CellChat object from spatial imaging data
##D # Input is a data matrix
##D cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
##D                            datatype = "spatial", coordinates = coordinates, scale.factors = scale.factors)
##D 
##D # input is a Seurat object
##D cellChat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "SCT",
##D                            datatype = "spatial", scale.factors = scale.factors)
##D 
## End(Not run)

## CellChat::computeCommunProb

# LLM-Generated Examples (official documentation not available)

# Example 1: Basic computation without considering population size
# Typical for analyzing communication strength independent of cell abundance
cellchat <- computeCommunProb(cellchat, 
                              population.size = FALSE)  # Don't weight by cell numbers

# Example 2: Advanced computation with population weighting and custom parameters
# Useful when cell type abundance matters for biological interpretation
cellchat <- computeCommunProb(cellchat, 
                              type = "truncatedMean",   # Use truncated mean method
                              trim = 0.1,               # Trim 10% of extreme values
                              population.size = TRUE,   # Weight by cell population size
                              raw.use = FALSE)          # Use normalized data (default)

## CellChat::subsetCommunication

## Not run: 
##D # access all the inferred cell-cell communications
##D df.net <- subsetCommunication(cellchat)
##D 
##D # access all the inferred cell-cell communications at the level of signaling pathways
##D df.net <- subsetCommunication(cellchat, slot.name = "netP")
##D 
##D # Subset to certain cells with sources.use and targets.use
##D df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
##D 
##D # Subset to certain signaling, e.g., WNT and TGFb
##D df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
## End(Not run)

# The examples above are for reference only. If they are not helpful, please disregard them entirely!
```

---

## Step 2: ReAct - Execute Task

Use the ReAct LLM to execute the task with the retrieved documentation:

```r
# Run ReAct workflow with documentation-augmented prompt
react_r(
  chat_obj = react_llm,  # Use ReAct model for task execution
  task = if (nchar(ret_docs_prompt) > 0) {
    paste0(task_prompt, "\n\n\n", ret_docs_prompt)  # Combine task with documentation
  } else {
    task_prompt
  }
)
```

---

## All-in-One: AutoFlow

The complete workflow above can be simplified into a single AutoFlow call:

```r
# AutoFlow automatically handles RAG + ReAct workflow
AutoFlow(
  react_llm = react_llm,      # Model for task execution
  task_prompt = task_prompt,  # Your analysis task
  rag_llm = rag_llm          # Model for documentation retrieval
)
```

---

## Summary

**AutoFlow Architecture:**

1. **RAG Phase** (`rag_llm`): 
   - Analyzes task requirements
   - Identifies relevant functions
   - Retrieves documentation and examples

2. **ReAct Phase** (`react_llm`):
   - Receives task + documentation
   - Executes reasoning and action loops
   - Returns final result
