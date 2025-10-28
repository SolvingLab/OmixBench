#' ReAct (Reasoning and Acting) using R code execution - Optimized Version
#'
#' @param chat_obj Chat object from ellmer
#' @param task Character string. The task description to be solved
#' @param max_turns Integer. Maximum number of ReAct turns (default: 15)
#' @param pkgs_to_use Character vector. R packages to load in session
#' @param objects_to_use Named list. Objects to load in R session
#' @param existing_session Existing callr session to continue from (optional)
#' @param verbose Logical. Whether to print progress information
#' @param r_session_options List. Options for callr R session
#' @param context_window_size Integer. Maximum characters before history summary (default: 3000)
#' @param max_observation_length Integer. Maximum observation length (default: 800)
#' @param error_escalation_threshold Integer. Error count threshold for escalation (default: 3)
#' @return List with complete ReAct results
#' @export
react_using_r <- function(chat_obj,
                          task,
                          max_turns = 15,
                          pkgs_to_use = c(),
                          objects_to_use = list(),
                          existing_session = NULL,
                          verbose = TRUE,
                          r_session_options = list(),
                          context_window_size = 3000,
                          max_observation_length = 800,
                          error_escalation_threshold = 3) {
  # Record start time
  start_time <- Sys.time()

  # System configuration
  system_config <- list(
    max_turns = max_turns,
    context_window_size = context_window_size,
    max_observation_length = max_observation_length,
    error_escalation_threshold = error_escalation_threshold,
    packages_used = pkgs_to_use,
    r_session_options = r_session_options
  )

  # Validation
  stopifnot(
    is.character(task) && length(task) == 1,
    is.numeric(max_turns) && max_turns > 0,
    is.logical(verbose)
  )

  if (verbose) {
    if (requireNamespace("cli", quietly = TRUE)) {
      cli::cli_rule("ReAct Process Starting")
      cli::cli_alert_info("User Task: {.field {substr(task, 1, 80)}...}")
      cli::cli_alert_info("Max turns: {.val {max_turns}}")
      cli::cli_alert_info("Setting system prompt...")
      cli::cli_text("")
    } else {
      cat("=== Starting ReAct Process ===\n")
      cat("Task:", substr(task, 1, 80), "...\n")
      cat("Max turns:", max_turns, "\n\n")
    }
  }

  # Setup R session
  if (!requireNamespace("callr", quietly = TRUE)) install.packages("callr")

  session_result <- setup_evaluation_session(
    evaluate_code = TRUE,
    r_session_options,
    pkgs_to_use,
    objects_to_use,
    existing_session
  )

  r_session <- session_result$evaluation_session
  history <- list()

  # Set system prompt once at the beginning
  system_prompt <- build_system_prompt(task)
  chat_obj$set_system_prompt(system_prompt)

  if (verbose && requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_success("System prompt set ({nchar(system_prompt)} characters)")
    cli::cli_text("")
  }

  # Error tracking
  error_history <- list()
  guidance_queue <- list()

  # Main ReAct loop
  for (turn in 1:max_turns) {
    if (verbose) {
      cli::cli_rule("Turn {turn}")
      cli::cli_alert_info("Preparing prompt for LLM...")
    }

    # 1. Smart history management
    managed_history <- manage_history_context(
      history,
      context_window_size,
      chat_obj,
      verbose
    )

    # 2. Add any pending guidance to context
    if (length(guidance_queue) > 0) {
      guidance_content <- paste(guidance_queue, collapse = "\n\n")
      guidance_queue <- list() # Clear after use

      if (verbose) cli::cli_alert_info("Adding guidance to context")
      managed_history <- add_guidance_to_history(managed_history, guidance_content)
    }

    # 3. Get JSON from LLM
    json_response <- get_react_json(chat_obj, managed_history, verbose)

    # 4. Parse response
    thought <- json_response$thought
    action <- json_response$action

    if (verbose) {
      cli::cli_h2("LLM Response")
      cli::cli_text("{.strong Thought:} {thought}")
      if (nchar(action) > 0) {
        cli::cli_text("{.strong Action:} R code provided")
        cli::cli_code(action)
      } else {
        cli::cli_text("{.strong Action:} {.emph No R code}")
      }
      cli::cli_text("")
    }

    # 5. Check if finished
    finished <- grepl("Finished!", thought, ignore.case = TRUE)
    if (verbose && finished) cli::cli_alert_success("LLM indicated completion")

    # 6. Execute R code
    execution <- execute_r_action(action, r_session, verbose, max_observation_length)
    r_session <- execution$session

    # 7. Record turn
    current_turn <- list(
      turn = turn,
      thought = thought,
      action = action,
      observation = execution$observation,
      success = execution$success,
      result = execution$result,
      timestamp = Sys.time()
    )
    history[[turn]] <- current_turn

    # 8. Error handling and escalation
    if (!execution$success) {
      error_history <- append(error_history, list(list(
        turn = turn,
        error = execution$observation,
        action = action
      )))

      # Check if error escalation is needed
      if (length(error_history) >= error_escalation_threshold) {
        recent_errors <- tail(error_history, error_escalation_threshold)
        if (detect_error_pattern(recent_errors)) {
          if (verbose) cli::cli_alert_warning("Error escalation triggered")

          # Generate error guidance
          error_guidance <- generate_error_guidance(recent_errors)
          guidance_queue <- append(guidance_queue, error_guidance)
        }
      }

      # Check if task degradation is needed (6 consecutive errors)
      if (length(error_history) >= 6) {
        recent_6_errors <- tail(error_history, 6)
        if (all(sapply(recent_6_errors, function(e) !is.null(e$error)))) {
          if (verbose) cli::cli_alert_warning("Task degradation suggested")

          # Generate task degradation guidance
          degradation_guidance <- generate_task_degradation_guidance(task)
          guidance_queue <- append(guidance_queue, degradation_guidance)
        }
      }
    }

    # 9. Exit if finished AND execution was successful AND result exists
    if (finished && execution$success && !is.null(execution$result)) {
      if (verbose && requireNamespace("cli", quietly = TRUE)) {
        cli::cli_rule("🎯 Task Completed Successfully")
        cli::cli_alert_success("LLM indicated completion with valid result")

        cli::cli_text("")
        cli::cli_h2("Final Answer")
        if (is.character(execution$result) && length(execution$result) == 1) {
          cli::cli_alert_success("Result: {.val {execution$result}}")
        } else {
          cli::cli_alert_success("Result: {.val {format_result_summary(execution$result)}}")
        }
        cli::cli_text("")
      }
      return(build_final_result(task, history, r_session,
        completed = TRUE,
        start_time = start_time,
        system_config = system_config
      ))
    }

    # 10. Approaching turn limit - activate result protection
    if (turn >= max_turns - 2) {
      protection_guidance <- handle_intermediate_result_protection(
        history, turn, max_turns, task
      )

      if (!is.null(protection_guidance)) {
        guidance_queue <- append(guidance_queue, protection_guidance)
        if (verbose) cli::cli_alert_warning("Result protection activated")
      }
    }

    # 11. Check natural stopping
    if (should_stop_naturally(history, turn)) {
      if (verbose) {
        cli::cli_rule("ReAct Process Stopped Naturally")
        cli::cli_alert_info("Natural stopping condition detected")
      }
      return(build_final_result(task, history, r_session,
        completed = FALSE,
        start_time = start_time,
        system_config = system_config
      ))
    }

    if (verbose) cli::cli_text("")
  }

  # Max turns reached
  if (verbose) {
    cli::cli_rule("Max Turns Reached")
    cli::cli_alert_warning("Reached maximum {max_turns} turns")
  }
  return(build_final_result(task, history, r_session,
    completed = FALSE,
    start_time = start_time,
    system_config = system_config
  ))
}

#' Smart history management
#' @param history Current history
#' @param context_window_size Maximum context size
#' @param chat_obj Chat object for summarization
#' @param verbose Verbose flag
#' @return Managed history
#' @noRd
manage_history_context <- function(history, context_window_size, chat_obj, verbose) {
  if (length(history) == 0) {
    return(history)
  }

  # Calculate current history length
  current_prompt <- build_history_prompt(history)
  current_length <- nchar(current_prompt)

  if (current_length <= context_window_size) {
    return(history)
  }

  if (verbose) {
    cli::cli_alert_info("History length ({current_length}) exceeds limit ({context_window_size})")
    cli::cli_alert_info("Generating history summary...")
  }

  # Generate history summary
  history_summary <- generate_history_summary(history, chat_obj, verbose)

  # Keep recent 2-3 turns + summary
  recent_turns <- tail(history, 3)

  # Create summary turn
  summary_turn <- list(
    turn = "SUMMARY",
    thought = "History Summary",
    action = "",
    observation = history_summary,
    success = TRUE,
    result = NULL,
    timestamp = Sys.time(),
    is_summary = TRUE
  )

  managed_history <- c(list(summary_turn), recent_turns)

  if (verbose) {
    cli::cli_alert_success("History summarized and compressed")
  }

  return(managed_history)
}

#' Generate history summary
#' @param history Full history
#' @param chat_obj Chat object
#' @param verbose Verbose flag
#' @return Summary string
#' @noRd
generate_history_summary <- function(history, chat_obj, verbose) {
  # Build summary request
  summary_prompt <- paste0(
    "Summarize the following ReAct conversation history, focusing on:\n",
    "1. Data status: what variables/objects exist, what they represent\n",
    "2. Analysis progress: which steps completed, current position in analysis workflow\n",
    "3. Key findings: important intermediate results or patterns\n",
    "4. Outstanding issues: errors encountered or next steps needed\n",
    "5. Method parameters: statistical methods used, thresholds, key settings\n\n",
    "Summary should be concise but complete, providing necessary context for continued analysis.\n\n",
    "History records:\n",
    build_full_history_text(history)
  )

  tryCatch(
    {
      summary_response <- chat_obj$chat(summary_prompt)
      return(summary_response)
    },
    error = function(e) {
      if (verbose) cli::cli_alert_warning("History summarization failed: {e$message}")
      return("History summary generation failed, continuing with recent conversation context.")
    }
  )
}

#' Build full history text for summarization
#' @param history History list
#' @return Full history text
#' @noRd
build_full_history_text <- function(history) {
  text_parts <- c()

  for (h in history) {
    if (isTRUE(h$is_summary)) next # Skip existing summaries

    turn_text <- paste0(
      "Turn ", h$turn, ":\n",
      "Thought: ", h$thought, "\n",
      "Action: ", if (nchar(h$action) > 0) "R code executed" else "No code", "\n",
      "Observation: ", h$observation, "\n"
    )
    text_parts <- append(text_parts, turn_text)
  }

  return(paste(text_parts, collapse = "\n"))
}

#' Add guidance to history context
#' @param history Current history
#' @param guidance_content Guidance text
#' @return History with guidance added
#' @noRd
add_guidance_to_history <- function(history, guidance_content) {
  guidance_turn <- list(
    turn = "GUIDANCE",
    thought = "System Guidance",
    action = "",
    observation = guidance_content,
    success = TRUE,
    result = NULL,
    timestamp = Sys.time(),
    is_guidance = TRUE
  )

  return(c(history, list(guidance_turn)))
}

#' Detect error patterns
#' @param error_history Recent error history
#' @return Boolean indicating if pattern detected
#' @noRd
detect_error_pattern <- function(error_history) {
  if (length(error_history) < 3) {
    return(FALSE)
  }

  # Check for similar error types
  error_types <- sapply(error_history, function(e) {
    if (grepl("syntax|Invalid", e$error)) {
      return("syntax")
    }
    if (grepl("object.*not found|could not find", e$error)) {
      return("missing_object")
    }
    if (grepl("subscript|index", e$error)) {
      return("indexing")
    }
    if (grepl("argument|function", e$error)) {
      return("function_call")
    }
    return("other")
  })

  # If recent 3 errors have 2+ of same type, consider it a pattern
  return(length(unique(error_types)) <= 2)
}

#' Generate error guidance
#' @param error_history Recent error history
#' @return Error guidance string
#' @noRd
generate_error_guidance <- function(error_history) {
  # Analyze error types
  error_types <- sapply(error_history, function(e) {
    if (grepl("syntax|Invalid", e$error)) {
      return("syntax")
    }
    if (grepl("object.*not found|could not find", e$error)) {
      return("missing_object")
    }
    if (grepl("subscript|index", e$error)) {
      return("indexing")
    }
    if (grepl("argument|function", e$error)) {
      return("function_call")
    }
    return("other")
  })

  dominant_type <- names(sort(table(error_types), decreasing = TRUE))[1]

  guidance <- switch(dominant_type,
    "syntax" = paste0(
      "=== ERROR ESCALATION GUIDANCE ===\n",
      "Repeated syntax errors detected. Suggestions:\n",
      "1. Break complex statements into simpler steps\n",
      "2. Check parentheses and quote matching\n",
      "3. Test core logic first, then add complexity\n",
      "4. Use basic R functions instead of complex wrappers"
    ),
    "missing_object" = paste0(
      "=== ERROR ESCALATION GUIDANCE ===\n",
      "Repeated object not found errors. Suggestions:\n",
      "1. Check current environment: ls()\n",
      "2. Verify variable name spelling\n",
      "3. Confirm data was loaded successfully\n",
      "4. Validate each variable's existence stepwise"
    ),
    "indexing" = paste0(
      "=== ERROR ESCALATION GUIDANCE ===\n",
      "Repeated indexing errors detected. Suggestions:\n",
      "1. Check data structure: str(data)\n",
      "2. Verify dimensions: dim(), length()\n",
      "3. Use safe indexing methods\n",
      "4. Test with small data subset first"
    ),
    "function_call" = paste0(
      "=== ERROR ESCALATION GUIDANCE ===\n",
      "Repeated function call errors. Suggestions:\n",
      "1. Check function arguments carefully\n",
      "2. Confirm required packages are loaded\n",
      "3. Review function help: ?function_name\n",
      "4. Try simpler alternative approaches"
    ),
    paste0(
      "=== ERROR ESCALATION GUIDANCE ===\n",
      "Repeated errors detected. General suggestions:\n",
      "1. Simplify current approach\n",
      "2. Debug step by step\n",
      "3. Check data and variable states\n",
      "4. Consider alternative implementation strategies"
    )
  )

  return(guidance)
}

#' Generate task degradation guidance
#' @param original_task Original task description
#' @return Task degradation guidance
#' @noRd
generate_task_degradation_guidance <- function(original_task) {
  guidance <- paste0(
    "=== TASK DEGRADATION SUGGESTION ===\n",
    "Current approach encountering repeated difficulties. Consider:\n\n",
    "1. **Simplify implementation**: Can you use more basic methods to achieve the same goal?\n",
    "   - Avoid complex wrapper functions, use basic R operations\n",
    "   - Break into steps: get partial results first, then refine\n\n",
    "2. **Reduce precision requirements**: While maintaining core objectives:\n",
    "   - Reduce data processing complexity\n",
    "   - Use approximate methods instead of exact algorithms\n\n",
    "3. **Partial completion strategy**:\n",
    "   - Provide best answer based on current results\n",
    "   - Clearly state which steps were completed\n",
    "   - Explain what's missing if result is incomplete\n\n",
    "Original task: ", original_task, "\n",
    "Please reconsider the problem with a more direct, stable solution approach."
  )

  return(guidance)
}

#' Handle intermediate result protection
#' @param history Current history
#' @param turn Current turn number
#' @param max_turns Maximum turns
#' @param task Original task
#' @return Protection guidance
#' @noRd
handle_intermediate_result_protection <- function(history, turn, max_turns, task) {
  if (turn < max_turns - 2) {
    return(NULL)
  }

  # Find recent successful results
  successful_results <- list()
  for (i in length(history):1) {
    if (history[[i]]$success && !is.null(history[[i]]$result)) {
      successful_results <- append(successful_results, list(list(
        turn = history[[i]]$turn,
        result = history[[i]]$result,
        thought = history[[i]]$thought
      )), after = 0)
      if (length(successful_results) >= 3) break
    }
  }

  # Build summary without using %>%
  result_summaries <- sapply(successful_results, function(r) {
    paste0("   Turn ", r$turn, ": ", substr(r$thought, 1, 50), "...")
  })

  protection_guidance <- paste0(
    "=== INTERMEDIATE RESULT PROTECTION ===\n",
    "Approaching analysis turn limit. Please:\n\n",
    "1. **Current progress summary**:\n",
    paste(result_summaries, collapse = "\n"),
    "\n\n2. **Requirements**:\n",
    "   - Provide best answer based on current results\n",
    "   - Clearly state which analysis steps were completed\n",
    "   - If result is incomplete, explain what's missing\n",
    "   - Must include specific result in 'Finished!' statement\n\n",
    "Please immediately compile current results and complete the task."
  )

  return(protection_guidance)
}

#' Enhanced system prompt builder
#' @param task Original task description
#' @return Enhanced system prompt string
#' @noRd
build_system_prompt <- function(task) {
  paste0(
    "You are an expert data analyst using the ReAct (Reasoning and Acting) method to solve problems.\n\n",
    "CURRENT TASK:\n", task, "\n\n",
    "WORKING METHOD:\n",
    "1. Analyze the current situation based on previous observations\n",
    "2. Plan your next action step by step\n",
    "3. Execute R code to gather information or perform analysis\n",
    "4. Use the results to inform your next reasoning step\n\n",
    "OUTPUT FORMAT:\n",
    "Always respond with JSON in this exact format:\n",
    "{\n",
    "  \"thought\": \"Your detailed analysis and planning\",\n",
    "  \"action\": \"R code to execute (empty string if no code needed)\"\n",
    "}\n\n",
    "COMPLETION RULE:\n",
    "- Include 'Finished!' in your thought ONLY when you have the exact final answer requested\n",
    "- Do NOT say 'Finished!' during data loading, preprocessing, or intermediate steps\n",
    "- Do NOT say 'Finished!' if there were any errors in the previous observation\n",
    "- Only say 'Finished!' when you can provide the precise result format requested\n\n",
    "ERROR HANDLING:\n",
    "- If previous observation shows ERROR, focus entirely on fixing it\n",
    "- If you receive system guidance (marked with ===), follow the suggestions carefully\n",
    "- When stuck on repeated errors, try simpler approaches to achieve the same goal\n",
    "- Consider breaking complex operations into smaller, testable steps\n\n",
    "IMPORTANT GUIDELINES:\n",
    "- If previous observation shows success, decide what to do next\n",
    "- Build upon the current R session state (variables persist)\n",
    "- Keep R code focused, executable, and error-free\n",
    "- Think step by step toward the final goal\n",
    "- Return results in the exact format requested by the user\n",
    "- When approaching turn limits, prioritize giving the best available answer"
  )
}

#' Modified build_history_prompt to handle summaries and guidance
#' @param history Execution history (may include summaries and guidance)
#' @return History prompt string
#' @noRd
build_history_prompt <- function(history) {
  if (length(history) == 0) {
    return("Begin your analysis.")
  }

  history_text <- "=== PREVIOUS TURNS ===\n"

  for (h in history) {
    if (isTRUE(h$is_summary)) {
      # Handle summary
      history_text <- paste0(
        history_text,
        "=== PREVIOUS ANALYSIS SUMMARY ===\n",
        h$observation, "\n\n"
      )
    } else if (isTRUE(h$is_guidance)) {
      # Handle guidance
      history_text <- paste0(
        history_text,
        h$observation, "\n\n"
      )
    } else {
      # Handle regular turn
      history_text <- paste0(
        history_text,
        sprintf("Turn %s:\n", as.character(h$turn)),
        "Thought: ", substr(h$thought, 1, 100), "...\n"
      )

      if (nchar(h$action) > 0) {
        history_text <- paste0(history_text, "Action: R code executed\n")
      } else {
        history_text <- paste0(history_text, "Action: No code\n")
      }

      # Show observation results
      if (h$success && !is.null(h$result)) {
        result_summary <- format_result_summary(h$result)
        history_text <- paste0(
          history_text,
          "Observation: SUCCESS - ", result_summary, "\n"
        )
      } else if (!h$success) {
        history_text <- paste0(
          history_text,
          "Observation: ERROR - ", substr(h$observation, 1, 150), "\n"
        )
      } else {
        history_text <- paste0(
          history_text,
          "Observation: ", substr(h$observation, 1, 100), "\n"
        )
      }

      history_text <- paste0(history_text, "\n")
    }
  }

  history_text <- paste0(
    history_text,
    "=== CONTINUE ANALYSIS ===\n",
    "Based on the above results, provide your next thought and action:"
  )

  return(history_text)
}

#' Modified execute_r_action with observation length control
#' @param action R code string
#' @param r_session Current R session
#' @param verbose Verbose flag
#' @param max_observation_length Maximum observation length
#' @return List with execution results
#' @noRd
execute_r_action <- function(action, r_session, verbose, max_observation_length = 800) {
  if (is.null(action) || nchar(trimws(action)) == 0) {
    if (verbose) {
      cli::cli_h3("R Execution")
      cli::cli_alert_info("No R code to execute")
    }
    return(list(
      observation = "No R code to execute",
      success = TRUE,
      result = NULL,
      session = r_session
    ))
  }

  if (verbose) {
    cli::cli_h3("R Execution")
    cli::cli_alert_info("Executing R code...")
  }

  # Clean R code
  clean_code <- clean_r_code(action)
  if (is.null(clean_code)) {
    if (verbose) cli::cli_alert_danger("Invalid R code syntax")
    return(list(
      observation = "Invalid R code syntax",
      success = FALSE,
      result = NULL,
      session = r_session
    ))
  }

  # Execute
  tryCatch(
    {
      session_clone <- r_session$clone()
      output <- session_clone$run_with_output(
        function(code) eval(parse(text = code), envir = .GlobalEnv),
        args = list(code = clean_code)
      )

      observation <- format_observation(output, max_observation_length)
      success <- is.null(output$error)

      if (verbose) {
        if (success) {
          cli::cli_alert_success("R code executed successfully")
        } else {
          cli::cli_alert_danger("R code execution failed")
        }
        cli::cli_h3("Observation")
        cli::cli_text("{observation}")

        if (!is.null(output$result) && success) {
          cli::cli_text("{.emph Result details:} {.val {class(output$result)[1]}} object")
        }
        if (!is.null(output$stderr) && nchar(trimws(output$stderr)) > 0) {
          cli::cli_text("{.emph Warnings:} {output$stderr}")
        }
      }

      return(list(
        observation = observation,
        success = success,
        result = output$result,
        session = session_clone
      ))
    },
    error = function(e) {
      obs <- paste("Execution error:", e$message)
      if (verbose) {
        cli::cli_alert_danger("Critical execution error")
        cli::cli_h3("Observation")
        cli::cli_text("{obs}")
      }
      return(list(
        observation = obs,
        success = FALSE,
        result = NULL,
        session = r_session
      ))
    }
  )
}

#' Format observation with length control
#' @param output R session output
#' @param max_length Maximum observation length
#' @return Formatted observation string
#' @noRd
format_observation <- function(output, max_length = 800) {
  if (!is.null(output$error)) {
    return(paste("ERROR:", output$error))
  }

  parts <- c()

  # Result (return value)
  if (!is.null(output$result)) {
    result_str <- format_result(output$result)
    if (nchar(result_str) > 0) {
      parts <- c(parts, paste("RESULT:", result_str))
    }
  }

  # Standard output
  if (!is.null(output$stdout) && nchar(trimws(output$stdout)) > 0) {
    stdout_str <- trimws(output$stdout)
    # Limit stdout length
    if (nchar(stdout_str) > 300) {
      stdout_str <- paste0(substr(stdout_str, 1, 297), "...")
    }
    parts <- c(parts, paste("OUTPUT:", stdout_str))
  }

  # Warnings
  if (!is.null(output$stderr) && nchar(trimws(output$stderr)) > 0) {
    stderr_str <- trimws(output$stderr)
    if (nchar(stderr_str) > 200) {
      stderr_str <- paste0(substr(stderr_str, 1, 197), "...")
    }
    parts <- c(parts, paste("WARNING:", stderr_str))
  }

  if (length(parts) == 0) {
    final_obs <- "Code executed successfully"
  } else {
    final_obs <- paste(parts, collapse = " | ")
  }

  # Overall length control
  if (nchar(final_obs) > max_length) {
    final_obs <- paste0(substr(final_obs, 1, max_length - 15), "... [truncated]")
  }

  return(final_obs)
}

#' Enhanced result formatting
#' @param result R object
#' @return Formatted string with length control
#' @noRd
format_result <- function(result) {
  if (is.null(result)) {
    return("")
  }

  tryCatch(
    {
      if (is.data.frame(result)) {
        if (nrow(result) <= 5 && ncol(result) <= 5) {
          # Small data frames can show more info
          return(paste0(
            "data.frame(", nrow(result), "×", ncol(result), ") with columns: ",
            paste(head(colnames(result), 3), collapse = ", "),
            if (ncol(result) > 3) "..." else ""
          ))
        } else {
          return(paste0("data.frame(", nrow(result), "×", ncol(result), ")"))
        }
      } else if (is.matrix(result)) {
        return(paste0("matrix(", nrow(result), "×", ncol(result), ")"))
      } else if (is.vector(result)) {
        if (length(result) <= 5) {
          return(paste("vector:", paste(result, collapse = ", ")))
        } else {
          return(paste0(
            "vector(length=", length(result), "): ",
            paste(head(result, 3), collapse = ", "), ", ..."
          ))
        }
      } else if (is.list(result)) {
        return(paste0("list(length=", length(result), ")"))
      } else {
        result_str <- as.character(result)
        if (nchar(result_str) > 100) {
          return(paste0(substr(result_str, 1, 97), "..."))
        } else {
          return(result_str)
        }
      }
    },
    error = function(e) "complex object"
  )
}

#' Format result summary
#' @param result R result object
#' @return Summary string
#' @noRd
format_result_summary <- function(result) {
  if (is.null(result)) {
    return("No result")
  }

  tryCatch(
    {
      if (is.numeric(result) && length(result) == 1) {
        return(paste("Value:", result))
      } else if (is.character(result) && length(result) == 1) {
        return(paste("String:", result))
      } else if (is.data.frame(result)) {
        return(paste0("data.frame(", nrow(result), "×", ncol(result), ")"))
      } else if (is.vector(result) && length(result) <= 5) {
        return(paste("Vector:", paste(result, collapse = ", ")))
      } else if (is.vector(result)) {
        return(paste0("Vector(", length(result), "): ", paste(head(result, 3), collapse = ", "), "..."))
      } else {
        return(paste("Object:", class(result)[1]))
      }
    },
    error = function(e) "Complex object"
  )
}

#' Check if should stop naturally
#' @param history Execution history
#' @param current_turn Current turn
#' @return Boolean
#' @noRd
should_stop_naturally <- function(history, current_turn) {
  if (length(history) < 3) {
    return(FALSE)
  }

  # Filter out system-generated content
  user_history <- history[sapply(history, function(h) !isTRUE(h$is_summary) && !isTRUE(h$is_guidance))]
  if (length(user_history) < 3) {
    return(FALSE)
  }

  recent <- tail(user_history, 3)

  # Multiple turns with no action
  no_action_count <- sum(sapply(recent, function(h) nchar(h$action) == 0))
  if (no_action_count >= 2) {
    return(TRUE)
  }

  # Multiple failures
  failure_count <- sum(sapply(recent, function(h) !h$success))
  if (failure_count >= 3) {
    return(TRUE)
  }

  # Repetitive thoughts
  thoughts <- sapply(tail(recent, 2), function(h) substr(h$thought, 1, 50))
  if (length(unique(thoughts)) == 1) {
    return(TRUE)
  }

  FALSE
}

#' Clean R code
#' @param code Raw code text
#' @return Cleaned code or NULL
#' @noRd
clean_r_code <- function(code) {
  if (is.null(code)) {
    return(NULL)
  }

  # Extract R code blocks if present
  r_blocks <- extract_r_code(code)
  if (length(r_blocks) > 0) {
    code <- paste(r_blocks, collapse = "\n")
  }

  code <- trimws(code)
  if (nchar(code) == 0) {
    return(NULL)
  }

  # Basic syntax check
  tryCatch(
    {
      parse(text = code)
      return(code)
    },
    error = function(e) NULL
  )
}

#' Build final result structure - Comprehensive version
#' @param task Original task
#' @param history Execution history
#' @param r_session Final R session
#' @param completed Whether completed normally
#' @param start_time Process start time
#' @param system_config System configuration
#' @return Comprehensive final result list
#' @noRd
build_final_result <- function(task, history, r_session, completed,
                               start_time = NULL, system_config = list()) {
  end_time <- Sys.time()

  # Get final answer from most recent successful result
  final_answer <- NULL
  for (i in length(history):1) {
    if (history[[i]]$success && !is.null(history[[i]]$result) &&
      !isTRUE(history[[i]]$is_summary) && !isTRUE(history[[i]]$is_guidance)) {
      final_answer <- history[[i]]$result
      break
    }
  }

  if (is.null(final_answer)) {
    final_answer <- "No conclusive result"
  }

  # Filter user turns (exclude system-generated content)
  user_turns <- history[sapply(history, function(h) !isTRUE(h$is_summary) && !isTRUE(h$is_guidance))]

  # Session info
  session_info <- tryCatch(
    {
      get_session_info(r_session)
    },
    error = function(e) {
      list(loaded_packages = character(0), defined_objects = character(0))
    }
  )

  # Performance metrics
  performance_metrics <- list(
    total_duration = if (!is.null(start_time)) as.numeric(difftime(end_time, start_time, units = "secs")) else NA,
    total_turns = length(user_turns),
    successful_turns = sum(sapply(user_turns, function(h) h$success)),
    failed_turns = sum(sapply(user_turns, function(h) !h$success)),
    total_tokens = sum(sapply(history, function(h) if (!is.null(h$tokens_total)) h$tokens_total else 0), na.rm = TRUE)
  )

  # Error analysis
  error_analysis <- analyze_errors(user_turns)

  # Code summary
  code_summary <- summarize_code_execution(user_turns)

  # Determine stop reason
  stop_reason <- determine_stop_reason(completed, user_turns, system_config$max_turns)

  # Collect optimization statistics
  optimization_info <- list(
    summary_used = any(sapply(history, function(h) isTRUE(h$is_summary))),
    summary_count = sum(sapply(history, function(h) isTRUE(h$is_summary))),
    guidance_provided = sum(sapply(history, function(h) isTRUE(h$is_guidance))),
    error_escalations = sum(sapply(history, function(h) isTRUE(h$is_guidance) && grepl("ERROR ESCALATION", h$observation, ignore.case = TRUE))),
    task_degradations = sum(sapply(history, function(h) isTRUE(h$is_guidance) && grepl("TASK DEGRADATION", h$observation, ignore.case = TRUE))),
    result_protection = any(sapply(history, function(h) isTRUE(h$is_guidance) && grepl("RESULT PROTECTION", h$observation, ignore.case = TRUE)))
  )

  # Build comprehensive result
  result <- list(
    # Core results
    task = task,
    final_answer = final_answer,
    completed = completed,
    success = completed || length(user_turns) > 0,
    stop_reason = stop_reason,

    # Execution details
    execution_history = history,
    total_turns = length(user_turns),

    # Performance metrics
    performance = performance_metrics,

    # Analysis summaries
    error_analysis = error_analysis,
    code_summary = code_summary,

    # System information
    final_session = r_session,
    session_info = session_info,
    system_config = system_config,

    # Enhancement statistics
    optimization_info = optimization_info,

    # Timestamps
    start_time = start_time,
    end_time = end_time
  )

  return(result)
}

#' Analyze errors from execution history
#' @param user_turns User turns only
#' @return Error analysis summary
#' @noRd
analyze_errors <- function(user_turns) {
  failed_turns <- user_turns[sapply(user_turns, function(h) !h$success)]

  if (length(failed_turns) == 0) {
    return(list(
      total_errors = 0,
      error_types = character(0),
      error_recovery_rate = 1.0,
      consecutive_errors = 0
    ))
  }

  # Classify error types
  error_types <- sapply(failed_turns, function(h) {
    obs <- h$observation
    if (grepl("syntax|Invalid", obs, ignore.case = TRUE)) {
      return("syntax")
    }
    if (grepl("object.*not found|could not find", obs, ignore.case = TRUE)) {
      return("missing_object")
    }
    if (grepl("subscript|index", obs, ignore.case = TRUE)) {
      return("indexing")
    }
    if (grepl("argument|function", obs, ignore.case = TRUE)) {
      return("function_call")
    }
    if (grepl("undefined columns|Can't subset", obs, ignore.case = TRUE)) {
      return("data_structure")
    }
    return("other")
  })

  # Calculate recovery rate (errors followed by success)
  recovery_count <- 0
  for (i in 1:(length(user_turns) - 1)) {
    if (!user_turns[[i]]$success && user_turns[[i + 1]]$success) {
      recovery_count <- recovery_count + 1
    }
  }
  recovery_rate <- if (length(failed_turns) > 0) recovery_count / length(failed_turns) else 1.0

  # Find longest consecutive error sequence
  consecutive_errors <- 0
  current_consecutive <- 0
  for (turn in user_turns) {
    if (!turn$success) {
      current_consecutive <- current_consecutive + 1
      consecutive_errors <- max(consecutive_errors, current_consecutive)
    } else {
      current_consecutive <- 0
    }
  }

  return(list(
    total_errors = length(failed_turns),
    error_types = table(error_types),
    error_recovery_rate = recovery_rate,
    consecutive_errors = consecutive_errors,
    failed_turns = sapply(failed_turns, function(h) h$turn)
  ))
}

#' Summarize code execution
#' @param user_turns User turns only
#' @return Code execution summary
#' @noRd
summarize_code_execution <- function(user_turns) {
  # Extract all executed code
  all_code <- c()
  successful_code <- c()
  failed_code <- c()

  for (turn in user_turns) {
    if (nchar(turn$action) > 0) {
      clean_code <- clean_r_code(turn$action)
      if (!is.null(clean_code)) {
        all_code <- c(all_code, clean_code)
        if (turn$success) {
          successful_code <- c(successful_code, clean_code)
        } else {
          failed_code <- c(failed_code, clean_code)
        }
      }
    }
  }

  return(list(
    total_code_blocks = length(all_code),
    successful_code_blocks = length(successful_code),
    failed_code_blocks = length(failed_code),
    complete_script = paste(successful_code, collapse = "\n\n"),
    all_attempted_code = paste(all_code, collapse = "\n\n# --- Next attempt ---\n"),
    execution_success_rate = if (length(all_code) > 0) length(successful_code) / length(all_code) else 0
  ))
}

#' Determine stop reason
#' @param completed Whether task completed successfully
#' @param user_turns User execution history
#' @param max_turns Maximum allowed turns
#' @return Stop reason string
#' @noRd
determine_stop_reason <- function(completed, user_turns, max_turns) {
  if (completed) {
    return("task_completed")
  } else if (length(user_turns) >= max_turns) {
    return("max_turns_reached")
  } else {
    # Check for natural stopping conditions
    if (length(user_turns) >= 3) {
      recent <- tail(user_turns, 3)

      # Multiple turns with no action
      no_action_count <- sum(sapply(recent, function(h) nchar(h$action) == 0))
      if (no_action_count >= 2) {
        return("no_progress_detected")
      }

      # Multiple failures
      failure_count <- sum(sapply(recent, function(h) !h$success))
      if (failure_count >= 3) {
        return("repeated_failures")
      }

      # Repetitive thoughts
      thoughts <- sapply(tail(recent, 2), function(h) substr(h$thought, 1, 50))
      if (length(unique(thoughts)) == 1) {
        return("repetitive_behavior")
      }
    }

    return("natural_stopping")
  }
}

#' Get session information helper function
#' @param r_session R session object
#' @return Session information list
#' @noRd
get_session_info <- function(r_session) {
  tryCatch(
    {
      session_data <- r_session$run(function() {
        list(
          loaded_packages = if (is.null(sessionInfo()$otherPkgs)) character(0) else names(sessionInfo()$otherPkgs),
          defined_objects = ls(envir = .GlobalEnv)
        )
      })
      return(session_data)
    },
    error = function(e) {
      return(list(loaded_packages = character(0), defined_objects = character(0)))
    }
  )
}

#' Modified get_react_json function with guidance support
#' @param chat_obj Chat object (with system prompt already set)
#' @param history Execution history (may include guidance)
#' @param verbose Verbose flag
#' @return Parsed JSON response
#' @noRd
get_react_json <- function(chat_obj, history, verbose) {
  # Build history prompt with guidance support
  history_prompt <- build_history_prompt(history)

  if (verbose && requireNamespace("cli", quietly = TRUE)) {
    cli::cli_h3("Prompt to LLM")
    cli::cli_text("{.emph History prompt ({nchar(history_prompt)} chars):}")
    if (nchar(history_prompt) > 200) {
      cli::cli_code(paste0(substr(history_prompt, 1, 200), "..."))
    } else {
      cli::cli_code(history_prompt)
    }
    cli::cli_text("")
    cli::cli_alert_info("Waiting for LLM response...")
  }

  tryCatch(
    {
      response_as_json(
        chat_obj = chat_obj,
        prompt = history_prompt,
        schema = list(
          type = "object",
          properties = list(
            thought = list(type = "string"),
            action = list(type = "string")
          ),
          required = c("thought", "action"),
          additionalProperties = FALSE
        ),
        schema_strict = TRUE,
        max_iterations = 2
      )
    },
    error = function(e) {
      if (verbose && requireNamespace("cli", quietly = TRUE)) {
        cli::cli_alert_danger("JSON parsing failed: {e$message}")
      }
      stop("Failed to get valid JSON from LLM")
    }
  )
}

#' Simplified interface - Enhanced react_r
#' @param chat_obj Chat object
#' @param task Task description
#' @param ... Additional arguments passed to react_using_r
#' @return Formatted result display
#' @export
react_r <- function(chat_obj, task, ...) {
  result <- react_using_r(chat_obj = chat_obj, task = task, ...)

  # Load cli for display
  if (!requireNamespace("cli", quietly = TRUE)) {
    install.packages("cli")
  }

  cli::cli_rule("ReAct Analysis Summary")
  cli::cli_text("")

  cli::cli_alert_info("Task: {.field {substr(result$task, 1, 80)}...}")
  cli::cli_alert_info("Status: {.val {ifelse(result$completed, 'Completed', 'Stopped')}}")
  cli::cli_alert_info("Total turns: {.val {result$total_turns}}")
  cli::cli_alert_info("Stop reason: {.val {result$stop_reason}}")

  # Show performance metrics
  if (!is.null(result$performance)) {
    cli::cli_text("")
    cli::cli_h3("Performance Metrics")
    cli::cli_alert_info("Duration: {.val {round(result$performance$total_duration, 1)}} seconds")
    cli::cli_alert_info("Success rate: {.val {round(result$performance$successful_turns / result$performance$total_turns * 100, 1)}}%")
    if (result$performance$total_tokens > 0) {
      cli::cli_alert_info("Total tokens: {.val {result$performance$total_tokens}}")
    }
  }

  # Show enhancement statistics
  if (!is.null(result$optimization_info)) {
    cli::cli_text("")
    cli::cli_h3("Enhancement Statistics")
    if (result$optimization_info$summary_used) {
      cli::cli_alert_info("Context summaries generated: {result$optimization_info$summary_count}")
    }
    if (result$optimization_info$guidance_provided > 0) {
      cli::cli_alert_info("Guidance interventions: {result$optimization_info$guidance_provided}")
    }
    if (result$optimization_info$error_escalations > 0) {
      cli::cli_alert_info("Error escalations: {result$optimization_info$error_escalations}")
    }
  }

  # Show error analysis if there were errors
  if (!is.null(result$error_analysis) && result$error_analysis$total_errors > 0) {
    cli::cli_text("")
    cli::cli_h3("Error Analysis")
    cli::cli_alert_info("Total errors: {result$error_analysis$total_errors}")
    cli::cli_alert_info("Recovery rate: {round(result$error_analysis$error_recovery_rate * 100, 1)}%")
    if (result$error_analysis$consecutive_errors > 0) {
      cli::cli_alert_info("Max consecutive errors: {result$error_analysis$consecutive_errors}")
    }
  }

  cli::cli_text("")

  if (!is.null(result$final_answer)) {
    cli::cli_h2("🎯 Final Answer")
    if (is.numeric(result$final_answer) && length(result$final_answer) == 1) {
      cli::cli_alert_success("{.val {result$final_answer}}")
    } else if (is.character(result$final_answer) && length(result$final_answer) == 1) {
      cli::cli_alert_success("{.val {result$final_answer}}")
    } else {
      cli::cli_text("{as.character(result$final_answer)}")
    }
    cli::cli_text("")
  }

  cli::cli_h2("Execution Timeline")
  for (turn in result$execution_history) {
    if (isTRUE(turn$is_summary) || isTRUE(turn$is_guidance)) next # Skip system-generated content

    status <- if (turn$success) cli::col_green("[OK]") else cli::col_red("[X]")
    cli::cli_text("{status} Turn {turn$turn}: {substr(turn$thought, 1, 60)}...")
    if (!turn$success) {
      cli::cli_text("   {cli::col_red('Error:')} {substr(turn$observation, 1, 80)}")
    } else if (!is.null(turn$result)) {
      cli::cli_text("   {cli::col_blue('Result:')} {substr(format_result_summary(turn$result), 1, 80)}")
    }
  }

  cli::cli_text("")
  cli::cli_rule()

  # Return invisible result for programmatic access
  invisible(result)
}
