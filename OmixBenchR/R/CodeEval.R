#' Evaluate Bioinformatics R Code Using LLM
#'
#' This function evaluates bioinformatics R code solutions using a Language Learning Model (LLM)
#' based on scientific correctness, technical quality, and professional standards. It returns
#' a structured JSON assessment with detailed scoring across multiple dimensions.
#'
#' @param task A character string containing the bioinformatics task or question that the code
#'   is supposed to solve. This provides context for the evaluation and helps the LLM assess
#'   whether the code addresses the biological question appropriately.
#'
#' @param response A character string containing the R code to be evaluated. This should be
#'   the complete code solution that attempts to solve the given bioinformatics task.
#'
#' @param llm_client An LLM provider object created by functions like \code{llm_openai()} or
#'   \code{llm_ollama()}. This object contains the configuration for connecting to and
#'   communicating with the specific LLM service that will perform the evaluation.
#'
#' @param max_retries Integer. Maximum number of retry attempts if the LLM fails to provide
#'   a valid JSON response (default: 3). The function will retry if JSON parsing fails or
#'   if the response doesn't match the expected schema structure.
#'
#' @param schema_strict Logical. Whether to enforce strict JSON schema validation (default: TRUE).
#'   When TRUE, the LLM response must exactly match the predefined evaluation schema with no
#'   additional properties allowed. When FALSE, allows more flexible JSON structure.
#'
#' @param schema_type Character. Method for enforcing JSON response format (default: "auto").
#'   Options include:
#'   \itemize{
#'     \item "auto": Automatically detect best method based on LLM provider
#'     \item "text-based": Add JSON instructions to prompt (works with any provider)
#'     \item "openai": Use OpenAI's native JSON mode (requires compatible OpenAI API)
#'     \item "ollama": Use Ollama's native JSON mode (requires compatible Ollama model)
#'     \item "openai_oo": OpenAI mode without schema enforcement in API
#'     \item "ollama_oo": Ollama mode without schema enforcement in API
#'   }
#'
#' @param verbose Logical. Whether to print detailed interaction logs to console (default: FALSE).
#'   When TRUE, shows the prompt being sent, the LLM's response, and any retry attempts.
#'   Useful for debugging and monitoring the evaluation process.
#'
#' @return A named list containing the evaluation results, or a list with an error message if
#'   evaluation fails. The successful return structure includes:
#'   \itemize{
#'     \item \code{total_score}: Integer (0-100) representing the overall evaluation score
#'     \item \code{breakdown}: Named list with individual dimension scores:
#'       \itemize{
#'         \item \code{problem_solving}: Score 0-20 for addressing the biological question
#'         \item \code{technical_implementation}: Score 0-25 for appropriate methods and workflow
#'         \item \code{data_handling}: Score 0-15 for preprocessing and quality control
#'         \item \code{statistical_rigor}: Score 0-15 for statistical methods and corrections
#'         \item \code{domain_knowledge}: Score 0-10 for biological context understanding
#'         \item \code{robustness}: Score 0-10 for error handling and validation
#'         \item \code{documentation}: Score 0-5 for code clarity and documentation
#'       }
#'   }
#'   If evaluation fails, returns \code{list(error = "error message")}.
#'
#' @details
#' This function implements a comprehensive evaluation framework specifically designed for
#' bioinformatics code assessment. The evaluation covers three main areas:
#'
#' \strong{Core Scientific Validity (45 points):}
#' \itemize{
#'   \item Problem solving approach and biological question addressing
#'   \item Technical implementation with appropriate bioinformatics methods
#'   \item Complete analysis workflows with logical flow
#' }
#'
#' \strong{Technical Quality (30 points):}
#' \itemize{
#'   \item Data handling including preprocessing and quality control
#'   \item Statistical rigor with proper corrections and methods
#' }
#'
#' \strong{Professional Excellence (25 points):}
#' \itemize{
#'   \item Domain knowledge and biological interpretation
#'   \item Code robustness and error handling
#'   \item Documentation and usability
#' }
#'
#' The function uses a positive scoring system where points are awarded for good practices
#' rather than deducted for shortcomings. This encourages comprehensive evaluation and
#' recognizes partial solutions that demonstrate scientific understanding.
#'
#' The evaluation leverages the \code{llmhelper::get_llm_response()} function with JSON
#' schema enforcement to ensure consistent, structured output that can be easily processed
#' and analyzed programmatically.
#'
#' @examples
#' \dontrun{
#' # Set up LLM client
#' library(llmhelper)
#' client <- llm_openai(
#'   base_url = "https://api.openai.com/v1/chat/completions",
#'   api_key = Sys.getenv("OPENAI_API_KEY"),
#'   model = "gpt-4",
#'   temperature = 0.2
#' )
#'
#' # Define a bioinformatics task
#' task <- "Analyze RNA-seq data to identify differentially expressed genes"
#'
#' # Example R code to evaluate
#' code <- "
#' library(DESeq2)
#' library(ggplot2)
#'
#' # Load count data
#' counts <- read.csv('counts.csv', row.names=1)
#' coldata <- read.csv('coldata.csv', row.names=1)
#'
#' # Create DESeq2 object
#' dds <- DESeqDataSetFromMatrix(countData = counts,
#'                               colData = coldata,
#'                               design = ~ condition)
#'
#' # Run DESeq2 analysis
#' dds <- DESeq(dds)
#' res <- results(dds)
#'
#' # Apply multiple testing correction
#' res$padj <- p.adjust(res$pvalue, method='BH')
#'
#' # Filter significant genes
#' sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
#' print(summary(res))
#' "
#'
#' # Evaluate the code
#' evaluation <- CodeEval(task, code, client, verbose = TRUE)
#'
#' # Print results
#' cat("Total Score:", evaluation$total_score, "/100\n")
#' cat("Problem Solving:", evaluation$breakdown$problem_solving, "/20\n")
#' cat("Technical Implementation:", evaluation$breakdown$technical_implementation, "/25\n")
#'
#' # Batch evaluation example
#' tasks <- c("Perform GO enrichment analysis", "Create a volcano plot")
#' codes <- c("library(clusterProfiler)...", "library(ggplot2)...")
#'
#' results <- mapply(CodeEval, tasks, codes,
#'   MoreArgs = list(llm_client = client),
#'   SIMPLIFY = FALSE
#' )
#' }
#'
#' @seealso
#' \code{\link[llmhelper]{get_llm_response}} for the underlying LLM communication function,
#' \code{\link[llmhelper]{llm_openai}} for creating OpenAI-compatible LLM providers,
#' \code{\link[llmhelper]{llm_ollama}} for creating Ollama LLM providers
#'
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
CodeEval <- function(task,
                     response,
                     llm_client,
                     max_retries = 3,
                     schema_strict = TRUE,
                     schema_type = "auto",
                     verbose = FALSE) {
  # Evaluation prompt template
  evaluation_prompt_template <- "## Evaluation Instructions

You will be provided with:
1. A bioinformatics task/question
2. R Code that attempts to solve this task

Evaluate the R code using the framework below. This is a positive-scoring system - award points for good practices, do not subtract.

## Scoring Framework (100 points total)

### Part 1: Core Scientific Validity (45 points)

#### 1.1 Problem Solving (20 points)
- [ ] Code addresses the biological question (+10)
- [ ] Solution approach is scientifically sound (+10)

#### 1.2 Technical Implementation (25 points)
- [ ] Appropriate methods for data type (+15)
  - RNA-seq: DESeq2/edgeR/limma or proper transformation
  - Genomics: GenomicRanges or appropriate coordinate handling
  - Clustering: suitable algorithm for data characteristics
  - General: method matches the biological question
- [ ] Complete analysis workflow (+10)
  - All necessary steps present
  - Logical flow from input to conclusion

### Part 2: Technical Quality (30 points)

#### 2.1 Data Handling (15 points)
- [ ] Appropriate preprocessing/normalization (+8)
  - Log transformation for expression data
  - Scaling for multi-omics integration
  - Batch effect consideration
- [ ] Data quality control (+7)
  - Checks for NA/missing values
  - Dimension validation
  - Sample/feature filtering

#### 2.2 Statistical Rigor (15 points)
- [ ] Multiple testing correction when needed (+8)
  - FDR/BH/Bonferroni for multiple comparisons
  - **If not applicable (e.g., clustering only): award 8 points**
- [ ] Appropriate statistical methods (+7)
  - Correct test for data distribution
  - Proper handling of paired/grouped data
  - **If only visualization/processing: award 7 points**

### Part 3: Professional Excellence (25 points)

#### 3.1 Domain Knowledge & Interpretation (10 points)
- [ ] Shows understanding of biological context (+5)
  - Explains why methods suit the biological data
  - Mentions relevant biological considerations
- [ ] Results are biologically interpretable (+5)
  - Output can be understood by biologists
  - Includes relevant biological annotation

#### 3.2 Robustness & Completeness (10 points)
- [ ] Error handling or input validation (+5)
  - Try-catch blocks or if-statements for edge cases
  - Informative messages or warnings
- [ ] Analysis validation or quality checks (+5)
  - Parameter optimization (e.g., choosing k for clustering)
  - Stability/reproducibility considerations (e.g., set.seed)

#### 3.3 Documentation & Usability (5 points)
- [ ] Clear workflow and outputs (+3)
  - Can follow the analysis logic
  - Results are saved or displayed
- [ ] Adequate documentation (+2)
  - Key steps explained
  - Parameters justified

## Evaluation Guidelines

1. **Focus on scientific merit** - Correct biology > elegant code
2. **Recognize good practices** - Award points for any professional elements
3. **Consider context** - Simple working solutions can score well if scientifically sound

## Task to evaluate:
{task}

## R Code to evaluate:
{response}

Please evaluate the code and return your assessment as a JSON object with total_score and breakdown of individual dimension scores."

  # Build the complete prompt
  full_prompt <- glue::glue(evaluation_prompt_template,
    task = task,
    response = response
  )

  # Get evaluation schema
  eval_schema <- get_evaluation_schema()

  # Call LLM with JSON schema enforcement
  tryCatch(
    {
      result <- llmhelper::get_llm_response(
        prompt = full_prompt,
        llm_client = llm_client,
        max_retries = max_retries,
        json_schema = eval_schema,
        schema_strict = schema_strict,
        schema_type = schema_type,
        verbose = verbose,
        return_mode = "only_response"
      )

      # Validate and adjust total score if needed
      if (!is.null(result) && !is.null(result$breakdown)) {
        breakdown_sum <- sum(unlist(result$breakdown))
        if (is.null(result$total_score) || result$total_score != breakdown_sum) {
          result$total_score <- breakdown_sum
        }
      }

      return(result)
    },
    error = function(e) {
      message("Error during evaluation: ", e$message)
      return(list(error = paste("Error during evaluation:", e$message)))
    }
  )
}

#' @noRd
#' @title Define the evaluation schema for JSON output
get_evaluation_schema <- function() {
  list(
    name = "bioinformatics_code_evaluation",
    description = "Schema for bioinformatics code evaluation results",
    schema = list(
      type = "object",
      properties = list(
        total_score = list(
          type = "integer",
          minimum = 0,
          maximum = 100,
          description = "Total evaluation score"
        ),
        breakdown = list(
          type = "object",
          properties = list(
            problem_solving = list(type = "integer", minimum = 0, maximum = 20),
            technical_implementation = list(type = "integer", minimum = 0, maximum = 25),
            data_handling = list(type = "integer", minimum = 0, maximum = 15),
            statistical_rigor = list(type = "integer", minimum = 0, maximum = 15),
            domain_knowledge = list(type = "integer", minimum = 0, maximum = 10),
            robustness = list(type = "integer", minimum = 0, maximum = 10),
            documentation = list(type = "integer", minimum = 0, maximum = 5)
          ),
          required = c(
            "problem_solving", "technical_implementation", "data_handling",
            "statistical_rigor", "domain_knowledge", "robustness", "documentation"
          ),
          additionalProperties = FALSE
        )
      ),
      required = c("total_score", "breakdown"),
      additionalProperties = FALSE
    )
  )
}
