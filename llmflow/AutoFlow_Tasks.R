library(llmflow)
library(jsonlite)

# Setup
current_dir <- getwd()
task_dirs <- list.files('OmixQA_Tasks', full.names = TRUE)

# ==================== LLM Configuration ====================

rag_llm <- chat_anthropic(model = 'claude-opus-4-20250514-thinking')  # For documentation
react_llm <- chat_deepseek(model = 'deepseek-chat')                   # For execution

# Results
resultL <- list()

# Main loop
for (task_path in task_dirs) {

  # Load task metadata
  task_meta <- tryCatch({
    setwd(task_path)
    fromJSON('task_prompt.json')
  }, error = function(e) {
    cat("\n⚠️  Skipping", basename(task_path), "-", e$message, "\n")
    return(NULL)
  })

  if (is.null(task_meta)) next

  task_prompt <- task_meta$Task_prompt
  qid <- task_meta$QID
  expected <- task_meta$Answer

  if (is.null(task_prompt) || is.null(qid)) next

  cat("\n========================================\n")
  cat("Task:", qid, "\n")
  cat("Expected:", expected, "\n")

  # Run AutoFlow with separate RAG and ReAct models
  response <- tryCatch({
    AutoFlow(
      react_llm = react_llm,
      task_prompt = task_prompt,
      rag_llm = rag_llm,  # Use specified RAG model
      verbose = FALSE
    )
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    list(
      final_answer = NA,
      error = e$message,
      total_turns = 0,
      performance = list(total_duration = 0)
    )
  })

  # Return to original dir
  setwd(current_dir)

  # Extract result
  actual <- response$final_answer

  # Smart comparison
  is_correct <- tryCatch({
    if (is.na(actual) || is.na(expected)) {
      FALSE
    } else {
      # Convert to character for comparison
      actual_str <- as.character(actual)
      expected_str <- as.character(expected)

      # Try exact match first
      if (actual_str == expected_str) {
        TRUE
      } else {
        # Try numeric comparison (with tolerance)
        actual_num <- suppressWarnings(as.numeric(actual_str))
        expected_num <- suppressWarnings(as.numeric(expected_str))

        if (!is.na(actual_num) && !is.na(expected_num)) {
          abs(actual_num - expected_num) < 0.01
        } else {
          # Try pattern matching for "Gene:Value" format
          if (grepl(":", expected_str) && grepl(":", actual_str)) {
            # Extract gene names
            expected_gene <- sub(":.*", "", expected_str)
            actual_gene <- sub(":.*", "", actual_str)

            if (expected_gene == actual_gene) {
              # Compare values
              expected_val <- suppressWarnings(as.numeric(sub(".*:", "", expected_str)))
              actual_val <- suppressWarnings(as.numeric(sub(".*:", "", actual_str)))

              if (!is.na(expected_val) && !is.na(actual_val)) {
                abs(expected_val - actual_val) < 0.01
              } else {
                FALSE
              }
            } else {
              FALSE
            }
          } else {
            FALSE
          }
        }
      }
    }
  }, error = function(e) FALSE)

  # Print result
  cat("Actual:", actual, "\n")
  cat("Correct:", is_correct, "\n")
  cat("Turns:", response$total_turns, "\n")
  cat("Duration:", round(response$performance$total_duration, 1), "s\n")

  # Store
  resultL[[qid]] <- list(
    qid = qid,
    expected = expected,
    actual = actual,
    is_correct = is_correct,
    turns = response$total_turns,
    duration = response$performance$total_duration,
    full_response = response
  )
}
