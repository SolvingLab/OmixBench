rm(list = ls())
library(OmixBenchR)
library(llmhelper)
library(jsonlite)

# Set working directory and get task directories
current_dir <- getwd()
task_dirs <- list.files('OmixQA_Tasks', full.names = TRUE)
task_dirs <- paste0(current_dir, '/', task_dirs)

# Read LLM configuration
# Input your LLM configuration file with three columns: model, api_key, and base_url.
llm_config <- read.csv('llm_config.csv', stringsAsFactors = FALSE)

# Main loop: iterate through each model
for (model_idx in 1:nrow(llm_config)) {

  # Extract configuration for current model
  current_model <- llm_config$model[model_idx]
  current_api_key <- llm_config$api_key[model_idx]
  current_base_url <- llm_config$base_url[model_idx]

  cat("\n========================================\n")
  cat("Processing Model:", current_model, "\n")
  cat("========================================\n")

  # Create LLM client for current model
  llm_client <- llm_provider(
    base_url = current_base_url,
    api_key = current_api_key,
    max_tokens = 15000,
    model = current_model
  )

  # Initialize result list for current model
  resultL <- list()

  # Loop through each task
  for (task_path in task_dirs) {

    # Change to task directory
    setwd(task_path)

    # Load task metadata
    task_meta <- fromJSON('task_prompt.json')
    task_prompt <- task_meta[['Task_prompt']]
    qid <- task_meta[['QID']]
    expected_answer <- task_meta[['Answer']]

    cat("\nProcessing Task:", qid, "\n")

    # Execute task with error handling
    response <- tryCatch({
      Execute_Task(
        task_prompt,
        llm_client = llm_client,
        timeout_sec = 1200,
        max_interactions = 10
      )
    }, error = function(e) {
      cat("Error in task", qid, ":", e$message, "\n")
      return(list(
        QID = qid,
        error = e$message,
        response = NULL,
        interactions = NA,
        duration_seconds = NA
      ))
    })

    # Store response
    resultL[[qid]] <- response

    # Extract and evaluate results
    if (!is.null(response) && !is.null(response[[qid]])) {

      # Extract result
      actual_result <- tryCatch({
        response[[qid]][["response"]][["output"]][["result"]]
      }, error = function(e) NA)

      # Extract interactions
      interactions <- tryCatch({
        length(response[[qid]][["interactions"]])
      }, error = function(e) NA)

      # Extract duration
      duration <- tryCatch({
        response[[qid]][["duration_seconds"]]
      }, error = function(e) NA)

      # Evaluate correctness (allowing small numeric tolerance)
      is_correct <- FALSE
      if (!is.na(actual_result) && !is.na(expected_answer)) {
        if (is.numeric(actual_result) && is.numeric(as.numeric(expected_answer))) {
          is_correct <- abs(as.numeric(actual_result) - as.numeric(expected_answer)) < 0.01
        } else {
          is_correct <- as.character(actual_result) == as.character(expected_answer)
        }
      }

      # Add evaluation metadata to response
      resultL[[qid]][["evaluation"]] <- list(
        expected_answer = expected_answer,
        actual_result = actual_result,
        is_correct = is_correct,
        num_interactions = interactions,
        duration_seconds = duration
      )

      # Print summary
      cat("  Expected:", expected_answer, "\n")
      cat("  Actual:", actual_result, "\n")
      cat("  Correct:", is_correct, "\n")
      cat("  Interactions:", interactions, "\n")
      cat("  Duration:", duration, "seconds\n")
    }
  }

  # Return to original directory
  setwd(current_dir)

  # Clean model name for filename (remove special characters)
  model_name_clean <- gsub("[^a-zA-Z0-9_-]", "_", current_model)

  # Save results for current model
  output_file <- paste0(model_name_clean, "_results.rds")
  saveRDS(resultL, file = output_file)

  cat("\n========================================\n")
  cat("Results saved to:", output_file, "\n")
  cat("========================================\n\n")

  # Generate summary statistics
  correct_count <- sum(sapply(resultL, function(x) {
    if (!is.null(x[["evaluation"]])) {
      return(x[["evaluation"]][["is_correct"]])
    }
    return(FALSE)
  }))

  total_tasks <- length(resultL)
  accuracy <- correct_count / total_tasks * 100

  cat("Summary for", current_model, ":\n")
  cat("  Total tasks:", total_tasks, "\n")
  cat("  Correct:", correct_count, "\n")
  cat("  Accuracy:", round(accuracy, 2), "%\n\n")
}
