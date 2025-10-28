#' Retrieve and Format R Function Documentation for LLM Consumption
#'
#' @param chat_obj LLM chat client object
#' @param prompt Task description
#' @param schema JSON schema (default: package_function_schema())
#' @param schema_strict Strict schema validation
#' @param skip_undocumented Skip functions without docs
#' @param use_llm_fallback Use LLM to generate examples
#' @param example_count Number of examples per function
#' @param warn_skipped Show warning for skipped functions
#' @return Formatted documentation string
#' @export
retrieve_docs <- function(chat_obj,
                          prompt,
                          schema = package_function_schema(), # 现在默认 min=0
                          schema_strict = TRUE,
                          skip_undocumented = TRUE,
                          use_llm_fallback = FALSE,
                          example_count = 2,
                          warn_skipped = TRUE) {
  # Extract functions using LLM
  function_response <- response_as_json(
    chat_obj = chat_obj,
    prompt = package_extraction_prompt(prompt),
    schema = schema,
    schema_strict = schema_strict
  )

  functions_vec <- function_response$functions

  # Handle empty result - return empty string
  if (length(functions_vec) == 0) {
    if (warn_skipped) {
      message("ℹ️  No specialized functions identified for this task. Using general R knowledge.")
    }
    return("")
  }

  # Parse package::function format
  functions_df <- data.frame(
    package = sapply(strsplit(functions_vec, "::"), `[`, 1),
    function_name = sapply(strsplit(functions_vec, "::"), `[`, 2),
    stringsAsFactors = FALSE
  )

  # Extract examples for all functions
  examples_list <- vector("list", nrow(functions_df))
  for (i in seq_len(nrow(functions_df))) {
    examples_list[[i]] <- extract_function_examples(
      package_name = functions_df$package[i],
      function_name = functions_df$function_name[i]
    )
  }

  # Filter out undocumented functions if requested
  if (skip_undocumented) {
    has_docs <- sapply(examples_list, function(x) {
      !is.na(x) && nchar(trimws(x)) > 0
    })

    skipped_functions <- functions_vec[!has_docs]

    if (length(skipped_functions) > 0 && warn_skipped) {
      message(sprintf(
        "ℹ️  Skipped %d function(s) without documentation: %s",
        length(skipped_functions),
        paste(skipped_functions, collapse = ", ")
      ))
    }

    functions_vec <- functions_vec[has_docs]
    functions_df <- functions_df[has_docs, ]
    examples_list <- examples_list[has_docs]

    # Return empty if no documented functions remain
    if (length(functions_vec) == 0) {
      if (warn_skipped) {
        message("ℹ️  No functions with documentation found. Using general R knowledge.")
      }
      return("")
    }
  }

  # Use LLM fallback for missing examples if enabled
  if (use_llm_fallback) {
    for (i in seq_along(examples_list)) {
      needs_generation <- is.na(examples_list[[i]]) ||
        nchar(trimws(examples_list[[i]])) == 0

      if (needs_generation) {
        func_spec <- paste0(functions_df$package[i], "::", functions_df$function_name[i])

        tryCatch(
          {
            cat("Generating examples for", func_spec, "...\n")

            llm_prompt <- paste0(
              "Generate ", example_count, " clear R code examples for: ",
              func_spec, "\n\n",
              "Requirements:\n",
              "- Provide realistic, working code examples\n",
              "- Include brief inline comments\n",
              "- Show different use cases if applicable\n",
              "- Use proper R syntax\n",
              "- Keep examples concise but complete\n\n",
              "Format as clean R code with comments."
            )

            llm_response <- as.character(chat_obj$chat(llm_prompt))

            examples_list[[i]] <- paste0(
              "# LLM-Generated Examples (official documentation not available)\n\n",
              llm_response
            )
          },
          error = function(e) {
            warning("LLM generation failed for ", func_spec, ": ", e$message)
            examples_list[[i]] <<- paste0(
              "No official usage examples available for ", func_spec
            )
          }
        )
      }
    }
  } else {
    # Replace NA or empty with standard message
    for (i in seq_along(examples_list)) {
      if (is.na(examples_list[[i]]) || nchar(trimws(examples_list[[i]])) == 0) {
        func_spec <- paste0(functions_df$package[i], "::", functions_df$function_name[i])
        examples_list[[i]] <- paste0("No official usage examples available for ", func_spec)
      }
    }
  }

  # Build final prompt text
  prompt_text <- paste0(
    "# Functions for reference\n",
    "You may use the following functions for this task (**If the functions are not listed quite correctly, please ignore the functional reference below.**):\n",
    paste0("- ", functions_vec, collapse = "\n"), "\n\n",
    "# Function Documentation & Examples\n\n",
    "Below are detailed examples and usage patterns for each function. ",
    "Use these to understand how to properly implement the analysis task:\n\n",
    paste0(
      sapply(seq_along(functions_vec), function(i) {
        paste0(
          "## ", functions_vec[i], "\n\n",
          examples_list[[i]], "\n"
        )
      }),
      collapse = "\n"
    ),
    "The examples above are for reference only. If they are not helpful, please disregard them entirely!"
  )

  return(prompt_text)
}


#' Create JSON Schema for Package Function Validation
#'
#' @param min_functions Integer. Minimum number of functions (default: 0 to allow empty)
#' @param max_functions Integer. Maximum number of functions (default: 10)
#' @param description Character string. Custom description
#' @return List containing JSON schema
#' @export
package_function_schema <- function(min_functions = 0, # 改为 0
                                    max_functions = 10,
                                    description = "Array of ONLY the most critical, domain-specific functions that truly require documentation (exclude common functions like read.csv, order, mean). Return empty array if no specialized functions are needed.") {
  # Validate inputs
  if (!is.numeric(min_functions) || min_functions < 0) {
    stop("min_functions must be a non-negative integer (>= 0)")
  }
  if (!is.numeric(max_functions) || max_functions < min_functions) {
    stop("max_functions must be >= min_functions")
  }
  if (max_functions > 50) {
    warning("max_functions > 50 may lead to poor documentation quality and context overflow")
  }

  schema <- list(
    type = "object",
    properties = list(
      functions = list(
        type = "array",
        items = list(
          type = "string",
          pattern = "^[a-zA-Z][a-zA-Z0-9._]*::[a-zA-Z][a-zA-Z0-9._]*$"
        ),
        description = description,
        maxItems = as.integer(max_functions)
      )
    ),
    required = c("functions"),
    additionalProperties = FALSE
  )

  # Only add minItems if > 0
  if (min_functions > 0) {
    schema$properties$functions$minItems <- as.integer(min_functions)
  }

  return(schema)
}



#' Generate Function Extraction Prompt for LLM Analysis
#'
#' Creates a highly refined prompt that guides LLMs to identify ONLY the most
#' documentation-critical, domain-specific R functions from a task description.
#' The prompt uses sophisticated filtering criteria to exclude common, well-known
#' functions (like read.csv, mean, order) that any LLM can use correctly without
#' explicit documentation, focusing instead on specialized functions where examples
#' truly add value.
#'
#' This function applies a "documentation necessity test": only include functions
#' where a proficient LLM would struggle without explicit documentation and examples.
#' This dramatically improves output quality and reduces token waste.
#'
#' @param task_description Character string. Detailed description of the R task
#'   or analysis workflow that needs to be performed. Should include:
#'   - Data types and sources involved
#'   - Analytical objectives and methods
#'   - Expected outputs or deliverables
#'   - Domain-specific context (e.g., bioinformatics, spatial analysis)
#'   The more domain-specific the description, the better the function selection.
#' @param include_criteria Character vector. Additional inclusion criteria
#'   beyond the defaults. Specify domain-specific requirements or function
#'   characteristics that should be documented. Default is NULL (use standard criteria).
#' @param exclude_criteria Character vector. Additional exclusion criteria
#'   beyond the defaults. Specify function types or patterns that should be
#'   skipped (e.g., "Basic ggplot2 themes", "Standard dplyr verbs"). Default is NULL.
#' @param prioritization_factors Character vector. Additional factors for
#'   prioritizing functions beyond the defaults. Specify what makes certain
#'   functions more important to document. Default is NULL (use standard priorities).
#' @param emphasis Character string. Additional emphasis or context to guide
#'   the extraction process. Use this to highlight specific aspects of the task
#'   or to emphasize certain types of functions. Default is NULL.
#'
#' @return Character string containing the complete extraction prompt with:
#'   \itemize{
#'     \item Clear documentation necessity principle
#'     \item Strict inclusion criteria for domain-specific functions
#'     \item Comprehensive exclusion rules with concrete examples
#'     \item Four-question decision heuristic for each function
#'     \item Concrete good/bad examples from multiple domains
#'     \item Prioritization by domain specialization and complexity
#'     \item Quality-over-quantity guidance
#'   }
#'
#' @details The enhanced prompt applies a rigorous "documentation necessity test"
#'   with four key questions:
#'
#'   1. Would a proficient LLM struggle without documentation?
#'   2. Is this function domain-specific or universally known?
#'   3. Does it use specialized terminology or workflows?
#'   4. Would examples significantly improve usage accuracy?
#'
#'   **Automatic exclusions** (common functions that waste tokens):
#'   - Data I/O: read.csv, write.csv, readLines
#'   - Basic operations: order, sort, subset, head, tail
#'   - Simple statistics: mean, median, sd, sum
#'   - Core structures: c, list, data.frame
#'   - Well-known tidyverse: simple dplyr::filter, dplyr::mutate
#'   - Basic control flow: if, for, while
#'   - Common utilities: paste, grep, unique
#'
#'   **What gets included** (documentation-critical functions):
#'   - Domain-specific methods (clusterProfiler::enrichGO for GO analysis)
#'   - Complex statistical procedures (DESeq2::DESeq)
#'   - Specialized transformations (sf::st_transform for spatial data)
#'   - Functions with many non-obvious parameters
#'   - Methods where wrong usage produces plausible but incorrect results
#'
#'   This approach ensures that "GO enrichment analysis" returns
#'   clusterProfiler functions, NOT read.csv or order.
#'
#' @examples
#' # Basic usage with improved prompt engineering
#' # Now automatically excludes common functions like read.csv, order, mean
#' prompt <- package_extraction_prompt(
#'   "Perform GO enrichment analysis on differentially expressed genes"
#' )
#' # Will focus on: clusterProfiler::enrichGO, org.Hs.eg.db::org.Hs.eg.db
#' # Will exclude: read.csv, dplyr::filter, order, mean
#'
#' # Compare: general data analysis (less domain-specific)
#' prompt1 <- package_extraction_prompt(
#'   "Clean and analyze time series sales data using tidyverse"
#' )
#' # Will likely return fewer functions since tidyverse basics are well-known
#'
#' # Compare: highly specialized domain
#' prompt2 <- package_extraction_prompt(
#'   "Spatial analysis of seismic events using terra and sf packages"
#' )
#' # Will focus on: terra::extract, sf::st_transform, sf::st_join
#' # Will exclude: basic ggplot2, dplyr operations
#'
#' # Add domain-specific guidance
#' prompt <- package_extraction_prompt(
#'   task_description = "Single-cell RNA-seq analysis with Seurat",
#'   include_criteria = c(
#'     "Seurat-specific normalization and scaling methods",
#'     "Dimensionality reduction functions unique to scRNA-seq"
#'   ),
#'   exclude_criteria = c(
#'     "Standard dplyr data manipulation",
#'     "Basic ggplot2 visualization"
#'   ),
#'   emphasis = "Focus on Seurat functions that require understanding of scRNA-seq
#'               methodology, not general R programming constructs."
#' )
#'
#' # Use with retrieve_docs - complete workflow
#' docs <- retrieve_docs(
#'   chat_obj = llm,
#'   prompt = package_extraction_prompt(
#'     task_description = "Perform differential expression and pathway enrichment",
#'     emphasis = "Exclude basic data manipulation, focus on statistical methods"
#'   )
#' )
#' # Result: Only domain-critical functions with examples
#'
#' @seealso \code{\link{retrieve_docs}} for using this prompt in documentation extraction
#'
#' @export
package_extraction_prompt <- function(task_description,
                                      include_criteria = NULL,
                                      exclude_criteria = NULL,
                                      prioritization_factors = NULL,
                                      emphasis = NULL) {
  # Build base prompt
  prompt <- paste0(
    "Extract ONLY the task-specific, domain-specialized R functions that truly require documentation.\n\n",
    "CRITICAL PRINCIPLE: Focus on functions where a proficient LLM would struggle WITHOUT explicit documentation.\n",
    "Exclude common, general-purpose functions that any LLM can use correctly without context.\n\n",
    "# TASK DESCRIPTION\n",
    task_description,
    "\n\n"
  )

  # Add inclusion criteria
  prompt <- paste0(
    prompt,
    "# INCLUSION CRITERIA\n",
    "ONLY include functions that meet ALL of these:\n",
    "1. Domain-specific functions unique to specialized packages (bioinformatics, spatial analysis, etc.)\n",
    "2. Functions with complex parameter configurations that require explicit examples\n",
    "3. Functions with non-obvious defaults or side effects specific to the domain\n",
    "4. Advanced statistical/analytical methods not commonly known\n",
    "5. Functions where incorrect usage would silently produce wrong results\n"
  )

  # Add custom inclusion criteria if provided
  if (!is.null(include_criteria) && length(include_criteria) > 0) {
    custom_include <- paste0(
      seq_along(include_criteria) + 5, ". ",
      include_criteria,
      collapse = "\n"
    )
    prompt <- paste0(prompt, custom_include, "\n")
  }

  prompt <- paste0(prompt, "\n")

  # Add exclusion criteria with MUCH more specificity
  prompt <- paste0(
    prompt,
    "# EXCLUSION CRITERIA (VERY IMPORTANT)\n",
    "EXCLUDE these common function categories:\n\n",
    "## Data I/O & Basic Operations\n",
    "- File reading/writing: read.csv, write.csv, readLines, read.table, etc.\n",
    "- Basic data manipulation: subset, transform, merge, head, tail, str\n",
    "- Basic sorting/ordering: sort, order, rank, rev\n",
    "- Basic subsetting: [, [[, $, which, subset\n\n",
    "## Basic Statistics & Math\n",
    "- Descriptive statistics: mean, median, sd, var, sum, min, max, range, quantile\n",
    "- Basic arithmetic: +, -, *, /, ^, %%, %/%\n",
    "- Common distributions: rnorm, dnorm, pnorm, runif, etc. (unless parameters are domain-critical)\n\n",
    "## Core Data Structures & Types\n",
    "- Constructors: c, list, data.frame, matrix, array, factor\n",
    "- Type checking: is.numeric, is.character, as.numeric, as.character, class, typeof\n",
    "- Length/dimension: length, nrow, ncol, dim, names, colnames, rownames\n\n",
    "## Common tidyverse Basics (well-known patterns)\n",
    "- Simple dplyr verbs: select, rename, arrange (unless used with complex selectors)\n",
    "- Basic %>% or |> piping\n",
    "- Simple ggplot2: ggplot, aes, geom_point, geom_line (unless very domain-specific styling)\n\n",
    "## Control Flow & Functions\n",
    "- if, else, for, while, repeat, break, next\n",
    "- function, return, stop, warning, message\n",
    "- apply, lapply, sapply, mapply (unless domain-specific usage pattern)\n\n",
    "## Other Common Utilities\n",
    "- paste, paste0, sprintf, cat, print\n",
    "- grep, grepl, sub, gsub, strsplit (basic string operations)\n",
    "- unique, duplicated, table, with, within\n"
  )

  # Add custom exclusion criteria if provided
  if (!is.null(exclude_criteria) && length(exclude_criteria) > 0) {
    custom_exclude <- paste0("- ", exclude_criteria, collapse = "\n")
    prompt <- paste0(prompt, "\n## Custom Exclusions\n", custom_exclude, "\n")
  }

  prompt <- paste0(prompt, "\n")

  # Add concrete examples
  prompt <- paste0(
    prompt,
    "# CONCRETE EXAMPLES\n\n",
    "## GOOD - These SHOULD be included:\n",
    "For GO enrichment analysis task:\n",
    "[OK] clusterProfiler::enrichGO (domain-specific, complex parameters)\n",
    "[OK] clusterProfiler::setReadable (domain-specific transformation)\n",
    "[OK] org.Hs.eg.db::org.Hs.eg.db (specialized annotation database)\n",
    "[OK] DESeq2::DESeq (complex statistical method)\n\n",
    "For spatial analysis task:\n",
    "[OK] sf::st_join (spatial-specific join with complex options)\n",
    "[OK] sf::st_transform (coordinate system transformation)\n",
    "[OK] terra::extract (spatial extraction with many parameters)\n\n",
    "## BAD - These should NOT be included:\n",
    "[X] read.csv (basic I/O, universally known)\n",
    "[X] dplyr::filter (too common, LLM knows well)\n",
    "[X] dplyr::mutate (too common, LLM knows well)\n",
    "[X] order (basic R, universally known)\n",
    "[X] mean (basic statistics)\n",
    "[X] ggplot2::ggplot (too common, unless very specific theming needed)\n\n"
  )

  # Add decision heuristic
  prompt <- paste0(
    prompt,
    "# DECISION HEURISTIC\n",
    "Ask yourself these questions:\n",
    "1. Would GPT-4 or Claude struggle to use this function correctly WITHOUT seeing its documentation?\n",
    "   - If NO -> EXCLUDE it\n",
    "   - If YES -> Consider including\n\n",
    "2. Is this function used across many domains with consistent behavior?\n",
    "   - If YES -> EXCLUDE it (it's too general)\n",
    "   - If NO -> Consider including\n\n",
    "3. Does this function have domain-specific terminology or workflows?\n",
    "   - If YES -> INCLUDE it\n",
    "   - If NO -> Likely exclude\n\n",
    "4. Would seeing code examples significantly improve usage accuracy for THIS specific task?\n",
    "   - If NO -> EXCLUDE it\n",
    "   - If YES -> Include it\n\n"
  )

  # Add prioritization
  prompt <- paste0(
    prompt,
    "# PRIORITIZATION\n",
    "When selecting functions, prioritize by:\n",
    "1. Domain specialization (how unique is this to the specific field?)\n",
    "2. Parameter complexity (how many non-obvious parameters?)\n",
    "3. Potential for silent errors (will wrong usage produce plausible but incorrect results?)\n",
    "4. Documentation necessity (will examples dramatically improve usage?)\n",
    "5. Task criticality (is this essential for the core analysis?)\n",
    "6. Uniqueness (are there no common equivalents across packages?)\n"
  )

  # Add custom prioritization factors if provided
  if (!is.null(prioritization_factors) && length(prioritization_factors) > 0) {
    custom_priority <- paste0(
      seq_along(prioritization_factors) + 6, ". ",
      prioritization_factors,
      collapse = "\n"
    )
    prompt <- paste0(prompt, custom_priority, "\n")
  }

  prompt <- paste0(prompt, "\n")

  # Add emphasis if provided
  if (!is.null(emphasis) && nchar(trimws(emphasis)) > 0) {
    prompt <- paste0(
      prompt,
      "# SPECIAL EMPHASIS\n",
      emphasis,
      "\n\n"
    )
  }

  # Add output format
  prompt <- paste0(
    prompt,
    "# OUTPUT FORMAT\n",
    "Return ONLY domain-specific, documentation-critical functions as a JSON object.\n",
    "Format: {\"functions\": [\"package::function\", ...]}\n\n",
    "Remember: Quality over quantity. 3-5 truly critical functions > 15 mixed functions.\n",
    "Each function you include should pass ALL the decision heuristic questions above.\n\n",
    "Example outputs:\n",
    "For \"GO enrichment analysis\":\n",
    '{\"functions\": [\"clusterProfiler::enrichGO\", \"clusterProfiler::enrichKEGG\", \"org.Hs.eg.db::org.Hs.eg.db\"]}\n\n',
    "For \"spatial point pattern analysis\":\n",
    '{\"functions\": [\"spatstat::ppp\", \"spatstat::Kest\", \"spatstat::envelope\"]}'
  )

  return(prompt)
}
