#' Execute Multiple Bioinformatics Tasks with Resume Capability (Dual Engine Support)
#'
#' @param tasks Named list where each element is a task prompt string
#' @param llm_client LLM client object. Type depends on execution_engine:
#'   - For "tidyprompt": LLM provider object
#'   - For "llmflow": Chat object from ellmer
#' @param execution_engine Execution engine to use: "tidyprompt" or "llmflow" (default: "tidyprompt")
#' @param timeout_sec Timeout in seconds (default: 1200)
#' @param max_interactions Maximum LLM interactions (default: 10)
#' @param return_mode Return mode specification:
#'   - For "tidyprompt": "full" or "only_response"
#'   - For "llmflow": "full", "code", "console", "object", "formatted_output", "llm_answer", "session"
#'   (default: "full" for both engines)
#' @param save_file Custom save file path (optional, auto-generated if NULL)
#' @param force_rerun Force re-execution of all tasks (default: FALSE)
#' @param pkgs_to_use (llmflow only) Packages to load in R session (default: c())
#' @param objects_to_use (llmflow only) Named list of objects to load in R session (default: list())
#' @param existing_session (llmflow only) Existing callr session to continue from (optional)
#' @param list_packages (llmflow only) Whether to list available packages in prompt (default: TRUE)
#' @param list_objects (llmflow only) Whether to list available objects in prompt (default: TRUE)
#' @param return_session_info (llmflow only) Whether to return session state information (default: TRUE)
#' @param evaluate_code (llmflow only) Whether to evaluate the generated code (default: TRUE)
#' @param r_session_options (llmflow only) Options for callr R session (default: list())
#' @param persist_session (llmflow only) Whether to persist R session across tasks (default: FALSE)
#'
#' @return Named list of completed task results
#' @export
#'
#' @examples
#' \dontrun{
#' # === tidyprompt mode ===
#' tasks <- list(
#'   task1 = "Calculate mean of iris$Sepal.Length",
#'   task2 = "Create a boxplot of iris data",
#'   task3 = "Run t-test on iris species"
#' )
#'
#' results <- Execute_Tasks(
#'   tasks = tasks,
#'   llm_client = my_llm_provider,
#'   execution_engine = "tidyprompt"
#' )
#'
#' # === llmflow mode ===
#' results <- Execute_Tasks(
#'   tasks = tasks,
#'   llm_client = my_chat_obj,
#'   execution_engine = "llmflow",
#'   pkgs_to_use = c("dplyr", "ggplot2"),
#'   persist_session = TRUE # Reuse session across tasks
#' )
#' }
Execute_Tasks <- function(tasks,
                          llm_client,
                          execution_engine = c("tidyprompt", "llmflow"),
                          timeout_sec = 1200,
                          max_interactions = 10,
                          return_mode = NULL,
                          save_file = NULL,
                          force_rerun = FALSE,
                          # llmflow specific parameters
                          pkgs_to_use = c(),
                          objects_to_use = list(),
                          existing_session = NULL,
                          list_packages = TRUE,
                          list_objects = TRUE,
                          return_session_info = TRUE,
                          evaluate_code = TRUE,
                          r_session_options = list(),
                          persist_session = FALSE) {
  # Input validation
  if (missing(llm_client)) stop("llm_client is required")
  if (!is.list(tasks) || is.null(names(tasks))) stop("tasks must be a named list")
  if (!all(sapply(tasks, is.character))) stop("each task must be a character string")

  # Validate and match execution_engine
  execution_engine <- match.arg(execution_engine)

  # Set default return_mode if not specified
  if (is.null(return_mode)) {
    return_mode <- "full"
  }

  # Determine save file
  if (is.null(save_file)) {
    if (execution_engine == "tidyprompt") {
      model_name <- llm_client$parameters$model %||% "default_model"
    } else {
      # For llmflow, try to extract model name from chat object
      model_name <- tryCatch(
        {
          llm_client$provider$model %||% "chat_model"
        },
        error = function(e) "chat_model"
      )
    }
    save_file <- paste0(model_name, "_results.rds")
  }

  # Load existing results
  results <- if (file.exists(save_file) && !force_rerun) {
    cat("📂 Loading existing results from:", save_file, "\n")
    readRDS(save_file)
  } else {
    cat("🆕 Starting fresh\n")
    list()
  }

  # Find remaining tasks
  remaining <- if (force_rerun) {
    names(tasks)
  } else {
    setdiff(names(tasks), names(results))
  }

  if (length(remaining) == 0) {
    cat("✅ All", length(tasks), "tasks already completed!\n")
    return(results)
  }

  cat("🔄 Processing", length(remaining), "remaining tasks\n")
  cat("🔧 Execution engine:", execution_engine, "\n")

  # Initialize session for llmflow if persist_session is TRUE
  current_session <- NULL
  if (execution_engine == "llmflow" && persist_session) {
    current_session <- existing_session
    cat("🔗 Session persistence enabled\n")
  }

  # Process remaining tasks
  for (i in seq_along(remaining)) {
    task_name <- remaining[i]
    cat(sprintf("\n⚡ [%d/%d] %s\n", i, length(remaining), task_name))

    # Execute task using Execute_Task
    task_result <- Execute_Task(
      task_prompt = tasks[[task_name]],
      task_name = task_name,
      llm_client = llm_client,
      execution_engine = execution_engine,
      timeout_sec = timeout_sec,
      max_interactions = max_interactions,
      return_mode = return_mode,
      save_file = NULL, # We handle saving at the batch level
      force_rerun = FALSE, # Already handled at batch level
      # llmflow specific parameters
      pkgs_to_use = pkgs_to_use,
      objects_to_use = objects_to_use,
      existing_session = current_session,
      list_packages = list_packages,
      list_objects = list_objects,
      return_session_info = return_session_info,
      evaluate_code = evaluate_code,
      r_session_options = r_session_options
    )

    # Update current_session if persisting (llmflow only)
    if (execution_engine == "llmflow" && persist_session &&
      !is.null(task_result) && is.list(task_result) && !is.null(task_result$session)) {
      current_session <- task_result$session
      cat("🔗 Session updated for next task\n")
    }

    # Save result (even if NULL)
    results[[task_name]] <- task_result
    saveRDS(results, save_file)
    cat("💾 Progress saved\n")
  }

  # Summary
  completed_count <- sum(!sapply(results, is.null))
  failed_tasks <- names(results)[sapply(results, is.null)]

  cat(sprintf("\n🎯 Finished! %d/%d tasks completed\n", completed_count, length(tasks)))

  if (length(failed_tasks) > 0) {
    cat("⚠️  Failed tasks:", paste(failed_tasks, collapse = ", "), "\n")
  }

  # Clean up session if it was created
  if (execution_engine == "llmflow" && persist_session && !is.null(current_session)) {
    tryCatch(
      {
        if (!is.null(existing_session)) {
          # Don't close if we were given an external session
          cat("🔗 Session kept alive (externally managed)\n")
        } else {
          current_session$close()
          cat("🔗 Session closed\n")
        }
      },
      error = function(e) {
        # Session might already be closed
      }
    )
  }

  results
}
