#' Execute Single Bioinformatics Task with Resume Capability (Dual Engine Support)
#'
#' @param task_prompt Character string containing the task prompt
#' @param task_name Optional task name for logging (default: auto-generated)
#' @param llm_client LLM client object. Type depends on execution_engine:
#'   - For "tidyprompt": LLM provider object
#'   - For "llmflow": Chat object from ellmer
#' @param execution_engine Execution engine to use: "tidyprompt" or "llmflow" (default: "tidyprompt")
#' @param timeout_sec Timeout in seconds (default: 1200)
#' @param max_interactions Maximum LLM interactions (default: 10)
#' @param return_mode Return mode specification:
#'   - For "tidyprompt": "full" or "only_response"
#'   - For "llmflow": "full", "code", "console", "object", "formatted_output", "llm_answer", "session"
#' @param save_file Optional save file path for result
#' @param force_rerun Force re-execution even if result exists (default: FALSE)
#' @param pkgs_to_use (llmflow only) Packages to load in R session (default: c())
#' @param objects_to_use (llmflow only) Named list of objects to load in R session (default: list())
#' @param existing_session (llmflow only) Existing callr session to continue from (optional)
#' @param list_packages (llmflow only) Whether to list available packages in prompt (default: TRUE)
#' @param list_objects (llmflow only) Whether to list available objects in prompt (default: TRUE)
#' @param return_session_info (llmflow only) Whether to return session state information (default: TRUE)
#' @param evaluate_code (llmflow only) Whether to evaluate the generated code (default: TRUE)
#' @param r_session_options (llmflow only) Options for callr R session (default: list())
#'
#' @return Task result object. Return structure depends on execution_engine and return_mode.
#' @export
#'
#' @examples
#' \dontrun{
#' # === tidyprompt mode ===
#' result <- Execute_Task(
#'   task_prompt = "Analyze this data",
#'   llm_client = my_llm_provider,
#'   execution_engine = "tidyprompt",
#'   return_mode = "full"
#' )
#'
#' # === llmflow mode ===
#' result <- Execute_Task(
#'   task_prompt = "Calculate mean of iris$Sepal.Length",
#'   llm_client = my_chat_obj,
#'   execution_engine = "llmflow",
#'   return_mode = "full",
#'   pkgs_to_use = c("dplyr"),
#'   objects_to_use = list(mydata = iris)
#' )
#' }
Execute_Task <- function(task_prompt,
                         task_name = NULL,
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
                         r_session_options = list()) {
  # Input validation
  if (missing(llm_client)) stop("llm_client is required")
  if (!is.character(task_prompt)) stop("task_prompt must be a character string")

  # Validate and match execution_engine
  execution_engine <- match.arg(execution_engine)

  # Set default return_mode based on execution_engine
  if (is.null(return_mode)) {
    return_mode <- if (execution_engine == "tidyprompt") "full" else "full"
  }

  # Validate return_mode based on execution_engine
  if (execution_engine == "tidyprompt") {
    return_mode <- match.arg(return_mode, choices = c("full", "only_response"))
  } else {
    return_mode <- match.arg(return_mode, choices = c(
      "full", "code", "console", "object",
      "formatted_output", "llm_answer", "session"
    ))
  }

  # Determine task name
  task_name <- task_name %||% paste0("task_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  # Check for existing result
  if (!is.null(save_file) && file.exists(save_file) && !force_rerun) {
    cat("📂 Loading existing result from:", save_file, "\n")
    result <- readRDS(save_file)
    cat("✅ Task already completed!\n")
    return(result)
  }

  cat("⚡ Executing task:", task_name, "\n")
  cat("🔧 Execution engine:", execution_engine, "\n")

  # Execute task with timeout protection based on execution_engine
  task_result <- tryCatch(
    {
      R.utils::withTimeout(
        {
          if (execution_engine == "tidyprompt") {
            # === tidyprompt mode ===
            task_prompt |>
              tidyprompt::add_text("IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).") |>
              tidyprompt::answer_using_r(evaluate_code = TRUE, output_as_tool = FALSE) |>
              tidyprompt::send_prompt(
                llm_client,
                max_interactions = max_interactions,
                return_mode = return_mode
              )
          } else {
            # === llmflow mode ===
            llmflow::response_to_r(
              chat_obj = llm_client,
              prompt = task_prompt,
              add_text = "IMPORTANT: Please use R language to solve this task. When calling any R functions, please use the package::function notation (for example, limma::normalizeBetweenArrays rather than just normalizeBetweenArrays).",
              pkgs_to_use = pkgs_to_use,
              objects_to_use = objects_to_use,
              existing_session = existing_session,
              list_packages = list_packages,
              list_objects = list_objects,
              return_session_info = return_session_info,
              evaluate_code = evaluate_code,
              r_session_options = r_session_options,
              return_mode = return_mode,
              max_iterations = max_interactions
            )
          }
        },
        timeout = timeout_sec
      )
    },
    TimeoutException = function(e) {
      cat("⏱️  Timeout after", timeout_sec, "seconds\n")
      NULL
    },
    error = function(e) {
      cat("❌ Error:", e$message, "\n")
      NULL
    }
  )

  # Save result if save_file is specified
  if (!is.null(save_file)) {
    saveRDS(task_result, save_file)
    cat("💾 Result saved to:", save_file, "\n")
  }

  # Status message
  if (!is.null(task_result)) {
    cat("✅ Task completed successfully!\n")

    # Show additional info for full mode (tidyprompt)
    if (execution_engine == "tidyprompt" && return_mode == "full" && is.list(task_result)) {
      if (!is.null(task_result$duration_seconds)) {
        cat("⏱️  Duration:", round(task_result$duration_seconds, 2), "seconds\n")
      }
      if (!is.null(task_result$interactions)) {
        cat("🔄 Interactions:", task_result$interactions, "\n")
      }
    }

    # Show additional info for llmflow mode
    if (execution_engine == "llmflow" && return_mode == "full" && is.list(task_result)) {
      if (!is.null(task_result$session_info)) {
        cat(
          "📦 Loaded packages:",
          if (length(task_result$session_info$loaded_packages) > 0) {
            paste(task_result$session_info$loaded_packages, collapse = ", ")
          } else {
            "none"
          }, "\n"
        )
        cat(
          "📊 Defined objects:",
          if (length(task_result$session_info$defined_objects) > 0) {
            length(task_result$session_info$defined_objects)
          } else {
            0
          }, "\n"
        )
      }
    }
  } else {
    cat("⚠️  Task failed or timed out\n")
  }

  task_result
}
