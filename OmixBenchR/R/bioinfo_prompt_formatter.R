#' Generate standardized bioinformatics analysis prompt from existing materials
#'
#' This function takes an existing prompt, expert code, and expected answer to generate
#' a standardized prompt following the bioinformatics analysis task format.
#'
#' @param current_prompt Character string. The original/current prompt text
#' @param gold_code Character string. The expert-written code that correctly solves the task
#' @param answer Character string. The expected answer (used for reference but not revealed in output)
#' @param llm_client LLM provider object created by llm_openai() or llm_ollama().
#' @return Character string containing the standardized prompt
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
#' @examples
#' \dontrun{
#' # Example usage
#' standardized_prompt <- bioinfo_prompt_formatter(
#'   current_prompt = "Find the most upregulated gene...",
#'   gold_code = "library(limma)\n...",
#'   answer = "GENE1"
#' )
#' }
bioinfo_prompt_formatter <- function(current_prompt,
                                     gold_code,
                                     answer,
                                     llm_client) {
  # Validate inputs
  if (!is.character(current_prompt) || nchar(trimws(current_prompt)) == 0) {
    stop("current_prompt must be a non-empty character string")
  }

  if (!is.character(gold_code) || nchar(trimws(gold_code)) == 0) {
    stop("gold_code must be a non-empty character string")
  }

  if (!is.character(answer) || nchar(trimws(answer)) == 0) {
    stop("answer must be a non-empty character string")
  }

  # Clean up line endings in gold_code
  gold_code <- gsub("\r\n", "\n", gold_code)

  # Construct the meta-prompt
  meta_prompt_template <- "You are a bioinformatics prompt engineer. Your task is to reformat existing prompts into a standardized format that will guide an LLM to:
1. Generate executable R code for the analysis
2. Return ONLY the final answer value (no code, no explanation)

## Input Information:
1. Original Prompt: {current_prompt}
2. Expert Code (Gold Code):
<EXPERT_CODE>
{gold_code}
</EXPERT_CODE>
3. Expected Answer Format: {answer} [EXAMPLE ONLY - DO NOT INCLUDE THIS VALUE]

## Output Format Requirements:

Your reformatted prompt must follow this EXACT structure:

### TASK: [One clear sentence describing what value to find/calculate]

### DATA SPECIFICATION:
- Dataset: [Dataset name and source]
- File path: '{filepath}filename.ext' (for local files)
- Data type: [e.g., gene expression matrix, methylation data]
- Platform: [e.g., Illumina 450K, RNA-seq] (if applicable)

### ANALYSIS PROCEDURE:
1. DATA LOADING:
   - [How to load the data]
   - [Required format/structure]

2. DATA PREPROCESSING:
   - [Key preprocessing steps]
   - [Essential quality control measures]

3. CORE ANALYSIS:
   - [Main analytical method]
   - [Statistical approach to use]

4. RESULT EXTRACTION:
   - [How to obtain the final answer]
   - [Selection/ranking criteria]

### OUTPUT SPECIFICATION:
- Result: [What the answer represents]
- Format: [Exact format required - be very specific]
- Final instruction: Return ONLY this value, no code or explanation

### CRITICAL PARAMETERS:
- [Only parameters essential for correct results]
- [Key thresholds that affect the outcome]

### REQUIREMENTS:
- Libraries: [Main R packages needed]
- Execution: Code must run without errors
- Output: Print only the final answer

## Reformatting Guidelines:

1. Focus on WHAT to do, not HOW to code it
2. Make the expected output format crystal clear
3. Include enough context for correct analysis without over-specifying
4. For file paths: local files use '{filepath}', URLs and package data keep original format
5. IMPORTANT: Do not expose the expert code details - include only data-loading patterns when necessary
6. DO NOT include the actual answer value anywhere in the prompt
7. Let the LLM demonstrate knowledge of standard bioinformatics practices
8. Ensure the task can be understood without seeing the expert code

Please reformat into the standardized format."

  # Format the prompt with actual values
  formatted_prompt <- glue::glue(meta_prompt_template,
    current_prompt = current_prompt,
    gold_code = gold_code,
    answer = answer
  )

  # Get LLM response
  standardized_prompt <- llmhelper::get_llm_response(
    prompt = formatted_prompt,
    llm_client = llm_client,
    max_retries = 3,
    verbose = FALSE,
    max_words = NULL,
    stream = FALSE
  )

  # Basic validation of the response
  if (is.null(standardized_prompt) || nchar(trimws(standardized_prompt)) == 0) {
    warning("LLM returned empty response. Please check inputs and try again.")
    return(NULL)
  }

  # Check if the response contains the expected sections
  required_sections <- c(
    "### TASK:", "### DATA SPECIFICATION:",
    "### ANALYSIS PROCEDURE:", "### OUTPUT SPECIFICATION:",
    "### CRITICAL PARAMETERS:", "### REQUIREMENTS:"
  )

  missing_sections <- required_sections[!sapply(
    required_sections,
    function(x) grepl(x, standardized_prompt, fixed = TRUE)
  )]

  if (length(missing_sections) > 0) {
    warning(paste(
      "The following sections might be missing:",
      paste(missing_sections, collapse = ", ")
    ))
  }

  # Check if answer was accidentally included
  if (grepl(answer, standardized_prompt, fixed = TRUE)) {
    warning("The standardized prompt may contain the actual answer. Please review and remove it.")
  }

  return(standardized_prompt)
}
