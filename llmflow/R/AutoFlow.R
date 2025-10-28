#' AutoFlow - Automated R Analysis Workflow with LLM
#'
#' @param react_llm Chat object for ReAct task execution (required)
#' @param task_prompt Task description (required)
#' @param rag_llm Chat object for RAG documentation retrieval (default: NULL, uses react_llm)
#' @param max_turns Maximum ReAct turns (default: 15)
#' @param pkgs_to_use Packages to load in R session
#' @param objects_to_use Named list of objects to load
#' @param existing_session Existing callr R session
#' @param verbose Verbose output (default: TRUE)
#' @param r_session_options Options for callr R session
#' @param context_window_size Context window size for history
#' @param max_observation_length Maximum observation length
#' @param error_escalation_threshold Error count threshold
#' @return ReAct result object
#' @export
#'
#' @details
#' **Dual-LLM Architecture:**
#'
#' AutoFlow supports using different models for different purposes:
#' - `rag_llm`: Retrieval-Augmented Generation - retrieves relevant function documentation
#' - `react_llm`: ReAct execution - performs reasoning and action loops
#'
#' **Why separate models?**
#' - RAG tasks are simple (extract function names) → use fast/cheap models
#' - ReAct tasks are complex (coding, reasoning) → use powerful models
#' - Cost savings: ~70% cheaper by using GPT-3.5 for RAG, GPT-4o for ReAct
#'
#' If `rag_llm` is NULL, both operations use `react_llm`.
#'
#' @examples
#' \dontrun{
#' # Simple: same model for both
#' llm <- llm_openai(model = "gpt-4o")
#' result <- AutoFlow(llm, "Load mtcars and plot mpg vs hp")
#'
#' # Optimized: lightweight RAG, powerful ReAct
#' rag <- llm_openai(model = "gpt-3.5-turbo") # Fast & cheap
#' react <- llm_openai(model = "gpt-4o") # Powerful
#' result <- AutoFlow(
#'   react_llm = react,
#'   task_prompt = "Perform PCA on iris dataset",
#'   rag_llm = rag
#' )
#'
#' # Cross-provider: DeepSeek RAG + Claude ReAct
#' rag <- chat_deepseek(model = "deepseek-chat")
#' react <- chat_anthropic(model = "claude-sonnet-4-20250514")
#' result <- AutoFlow(react, "Complex analysis", rag_llm = rag)
#'
#' # Batch evaluation with shared RAG
#' rag <- chat_deepseek(model = "deepseek-chat")
#' react <- chat_openai(model = "gpt-4o")
#'
#' for (task in tasks) {
#'   result <- AutoFlow(react, task, rag_llm = rag, verbose = FALSE)
#' }
#' }
AutoFlow <- function(
    react_llm,
    task_prompt,
    rag_llm = NULL,
    max_turns = 15,
    pkgs_to_use = c(),
    objects_to_use = list(),
    existing_session = NULL,
    verbose = TRUE,
    r_session_options = list(),
    context_window_size = 3000,
    max_observation_length = 800,
    error_escalation_threshold = 3) {
  # Use react_llm for RAG if rag_llm not provided
  if (is.null(rag_llm)) {
    rag_llm <- react_llm
  }

  # RAG: Retrieve documentation
  if (verbose) cat("\n📚 Retrieving function documentation...\n")

  ret_docs_prompt <- tryCatch(
    retrieve_docs(
      rag_llm, # Use RAG model for documentation retrieval
      prompt = task_prompt,
      skip_undocumented = FALSE,
      use_llm_fallback = TRUE,
      warn_skipped = FALSE
    ),
    error = function(e) {
      if (verbose) warning("⚠️  ", e$message, call. = FALSE)
      ""
    }
  )

  # Combine task with documentation
  final_prompt <- if (nchar(ret_docs_prompt) > 0) {
    paste0(task_prompt, "\n\n\n", ret_docs_prompt)
  } else {
    task_prompt
  }

  if (verbose) cat("\n🚀 Starting ReAct workflow...\n")

  # ReAct: Execute task
  react_r(
    react_llm, # Use ReAct model for task execution
    task = final_prompt,
    max_turns = max_turns,
    pkgs_to_use = pkgs_to_use,
    objects_to_use = objects_to_use,
    existing_session = existing_session,
    verbose = verbose,
    r_session_options = r_session_options,
    context_window_size = context_window_size,
    max_observation_length = max_observation_length,
    error_escalation_threshold = error_escalation_threshold
  )
}
