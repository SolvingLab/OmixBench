#' Response to R code generation and execution with session continuity
#'
#' @param chat_obj Chat object from ellmer
#' @param prompt User prompt for R code generation
#' @param add_text Additional instruction text
#' @param pkgs_to_use Packages to load in R session
#' @param objects_to_use Named list of objects to load in R session
#' @param existing_session Existing callr session to continue from (optional)
#' @param list_packages Whether to list available packages in prompt
#' @param list_objects Whether to list available objects in prompt
#' @param return_session_info Whether to return session state information
#' @param evaluate_code Whether to evaluate the generated code
#' @param r_session_options Options for callr R session
#' @param return_mode Return mode specification
#' @param max_iterations Maximum retry attempts
#' @return Result based on return_mode
#' @export
response_to_r <- function(chat_obj,
                          prompt,
                          add_text = NULL,
                          pkgs_to_use = c(),
                          objects_to_use = list(),
                          existing_session = NULL,
                          list_packages = TRUE,
                          list_objects = TRUE,
                          return_session_info = TRUE,
                          evaluate_code = TRUE,
                          r_session_options = list(),
                          return_mode = c("full", "code", "console", "object", "formatted_output", "llm_answer", "session"),
                          max_iterations = 3) {
  return_mode <- match.arg(return_mode)
  stopifnot(
    is.null(add_text) || (is.character(add_text) && length(add_text) == 1),
    length(pkgs_to_use) == 0 ||
      is.vector(pkgs_to_use) & all(sapply(pkgs_to_use, is.character)),
    is.list(objects_to_use),
    length(objects_to_use) == 0 || !is.null(names(objects_to_use)),
    is.null(existing_session) || inherits(existing_session, "r_session"),
    is.logical(list_packages),
    is.logical(list_objects),
    is.logical(return_session_info),
    is.logical(evaluate_code),
    is.list(r_session_options)
  )

  if (evaluate_code & !requireNamespace("callr", quietly = TRUE)) install.packages("callr")
  if (!evaluate_code & return_mode %in% c("console", "object", "formatted_output", "session")) {
    stop("The return mode must be 'full', 'code', or 'llm_answer' if 'evaluate_code' is FALSE")
  }

  ## Setup evaluation session and load packages/objects
  session_result <- setup_evaluation_session(evaluate_code, r_session_options, pkgs_to_use, objects_to_use, existing_session)
  evaluation_session <- session_result$evaluation_session
  loaded_pkgs <- session_result$loaded_pkgs
  loaded_objects <- session_result$loaded_objects
  existing_session_info <- session_result$existing_session_info

  # Record starting turn for retry context
  start_turn <- length(chat_obj$get_turns()) + 1

  # Modify prompt with additional context including session state
  prompt_text <- modify_prompt(
    prompt, add_text, list_packages, list_objects,
    loaded_pkgs, objects_to_use, evaluate_code, return_mode, existing_session_info
  )

  # Get LLM response
  res <- chat_obj$chat(prompt_text)

  return_list <- list()
  return_list$llm_answer <- res

  # Extract and retry extraction
  extracted_code <- extract_r_code(res)
  extracted_code <- retry_r_extraction(chat_obj, extracted_code, start_turn, max_iterations)

  # Parse and retry parsing
  parse_result <- retry_r_parsing(chat_obj, extracted_code, start_turn, max_iterations)
  extracted_code <- parse_result$extracted_code
  parsed_code <- parse_result$parsed_code

  if (!evaluate_code) {
    return_list$code <- parsed_code
    if (return_mode == "code") {
      return(parsed_code)
    }
    if (return_mode == "llm_answer") {
      return(res)
    }
    return(return_list)
  }

  # Run and retry execution
  run_result <- retry_r_run(
    chat_obj, extracted_code, parsed_code, start_turn, max_iterations,
    evaluation_session, return_mode
  )
  extracted_code <- run_result$extracted_code
  parsed_code <- run_result$parsed_code
  output <- run_result$output
  final_session <- run_result$final_session

  return_list$code <- parsed_code
  return_list$output <- output
  return_list$formatted_output <- format_output(extracted_code, output)
  return_list$session <- final_session

  # Return session state information
  if (return_session_info) {
    return_list$session_info <- get_session_info(final_session)
  }

  # Return based on mode
  if (return_mode == "full") {
    return(return_list)
  }
  if (return_mode == "code") {
    return(return_list$code)
  }
  if (return_mode == "console") {
    return(return_list$output$stdout)
  }
  if (return_mode == "object") {
    return(return_list$output$result)
  }
  if (return_mode == "formatted_output") {
    return(return_list$formatted_output)
  }
  if (return_mode == "llm_answer") {
    return(res)
  }
  if (return_mode == "session") {
    return(return_list$session)
  }

  return(output$stdout)
}

#' Format output for display
#'
#' @param extracted_code Extracted R code
#' @param output Execution output
#' @return Formatted output string
#' @noRd
format_output <- function(extracted_code, output) {
  # Helper function to truncate text
  str_trunc <- function(text, max_chars) {
    if (nchar(text) <= max_chars) {
      return(text)
    } else {
      return(paste0(substr(text, 1, max_chars - 3), "..."))
    }
  }

  # Helper function to safely convert any R object to string
  safe_object_to_string <- function(obj) {
    if (is.null(obj)) {
      return("NULL")
    }

    tryCatch(
      {
        # Handle different object types
        if (is.function(obj)) {
          # For functions, show the function signature and first few lines
          func_def <- deparse(obj)
          if (length(func_def) > 5) {
            # Show first few lines and indicate there's more
            paste(c(func_def[1:5], "# ... (function continues)"), collapse = "\n")
          } else {
            paste(func_def, collapse = "\n")
          }
        } else if (is.data.frame(obj) || is.matrix(obj)) {
          # For data frames and matrices, use print output
          paste(capture.output(print(obj)), collapse = "\n")
        } else if (is.list(obj) && !is.data.frame(obj)) {
          # For lists, use str() to show structure
          paste(capture.output(str(obj, max.level = 2)), collapse = "\n")
        } else if (is.vector(obj) && length(obj) > 0) {
          # For vectors, use paste safely
          if (length(obj) > 10) {
            # For long vectors, show first few elements
            paste(
              c(
                paste(obj[1:10], collapse = " "),
                paste("... (", length(obj) - 10, "more elements)")
              ),
              collapse = "\n"
            )
          } else {
            paste(obj, collapse = " ")
          }
        } else {
          # For other objects, use print output
          paste(capture.output(print(obj)), collapse = "\n")
        }
      },
      error = function(e) {
        # Fallback: show object class and length/dim info
        obj_info <- paste("Object of class:", paste(class(obj), collapse = ", "))
        if (has_length <- !is.null(try(length(obj), silent = TRUE))) {
          obj_info <- paste(obj_info, "| Length:", length(obj))
        }
        if (has_dim <- !is.null(try(dim(obj), silent = TRUE))) {
          obj_info <- paste(obj_info, "| Dimensions:", paste(dim(obj), collapse = "x"))
        }
        paste(obj_info, "| (Cannot display content)")
      }
    )
  }

  glue::glue(
    "--- R code: ---\n",
    "{paste(extracted_code, collapse = '\\n')}\n\n",
    "--- Console output: ---\n",
    "{
      if (is.null(output$stdout) || output$stdout == '') {
        'No console output produced.'
      } else {
        str_trunc(paste(output$stdout, collapse = '\\n'), 1000)
      }
    }\n\n",
    "--- Last object: ---\n",
    "{
      if (is.null(output$result)) {
        'No object produced.'
      } else {
        str_trunc(safe_object_to_string(output$result), 500)
      }
    }"
  )
}

#' Setup evaluation session and load packages/objects
#'
#' @param evaluate_code Whether to evaluate code
#' @param r_session_options Session options
#' @param pkgs_to_use Packages to load
#' @param objects_to_use Objects to load
#' @param existing_session Existing session to reuse
#' @return List with evaluation_session, loaded_pkgs, loaded_objects, existing_session_info
#' @noRd
setup_evaluation_session <- function(evaluate_code, r_session_options, pkgs_to_use, objects_to_use, existing_session = NULL) {
  existing_session_info <- NULL

  if (evaluate_code) {
    if (!is.null(existing_session)) {
      # Clone existing session safely
      evaluation_session <- existing_session$clone()

      # Get existing session information
      existing_session_info <- get_session_info(existing_session)
      loaded_pkgs <- existing_session_info$loaded_packages
      loaded_objects <- existing_session_info$defined_objects

      # Check and load new packages (avoid duplicates)
      new_pkgs <- setdiff(pkgs_to_use, loaded_pkgs)
      if (length(new_pkgs) > 0) {
        # Check if new packages are installed
        installed_pkgs <- evaluation_session$run(
          function(pkgs_to_check) {
            sapply(pkgs_to_check, function(pkg) {
              requireNamespace(pkg, quietly = TRUE)
            }, simplify = TRUE, USE.NAMES = TRUE)
          },
          args = list(pkgs_to_check = new_pkgs)
        )

        if (any(installed_pkgs == FALSE)) {
          stop(paste0(
            "The following packages are not installed: ",
            paste(names(installed_pkgs)[installed_pkgs == FALSE], collapse = ", ")
          ))
        }

        # Load new packages
        new_loaded_pkgs <- evaluation_session$run(
          function(pkgs_to_load) {
            for (pkg_name in pkgs_to_load) {
              library(pkg_name, character.only = TRUE)
            }
            session_info <- utils::sessionInfo()
            if (is.null(session_info$otherPkgs)) {
              character(0)
            } else {
              names(session_info$otherPkgs)
            }
          },
          args = list(pkgs_to_load = new_pkgs)
        )
        loaded_pkgs <- new_loaded_pkgs
      }

      # Load new objects (avoid duplicates)
      new_objects <- setdiff(names(objects_to_use), loaded_objects)
      if (length(new_objects) > 0) {
        new_objects_to_load <- objects_to_use[new_objects]
        new_loaded_objects <- evaluation_session$run(
          function(objects_to_load) {
            for (i in seq_along(objects_to_load)) {
              obj <- objects_to_load[[i]]
              obj_name <- names(objects_to_load)[i]
              assign(obj_name, obj, envir = .GlobalEnv)
            }
            ls(envir = .GlobalEnv)
          },
          args = list(objects_to_load = new_objects_to_load)
        )
        loaded_objects <- new_loaded_objects
      }
    } else {
      # Create new session (original logic)
      if (length(r_session_options) == 0) {
        r_session_options <- callr::r_session_options()
        r_session_options$system_profile <- FALSE
        r_session_options$user_profile <- FALSE
      }
      evaluation_session <- callr::r_session$new(options = r_session_options)

      # Check and load packages
      if (length(pkgs_to_use) > 0) {
        installed_pkgs <- evaluation_session$run(
          function(pkgs_to_use) {
            sapply(pkgs_to_use, function(pkg) {
              requireNamespace(pkg, quietly = TRUE)
            }, simplify = TRUE, USE.NAMES = TRUE)
          },
          args = list(pkgs_to_use = pkgs_to_use)
        )

        if (any(installed_pkgs == FALSE)) {
          stop(paste0(
            "The following packages are not installed: ",
            paste(names(installed_pkgs)[installed_pkgs == FALSE], collapse = ", ")
          ))
        }

        loaded_pkgs <- evaluation_session$run(
          function(pkgs_to_use) {
            for (pkg_name in pkgs_to_use) {
              library(pkg_name, character.only = TRUE)
            }
            session_info <- utils::sessionInfo()
            if (is.null(session_info$otherPkgs)) {
              character(0)
            } else {
              names(session_info$otherPkgs)
            }
          },
          args = list(pkgs_to_use = pkgs_to_use)
        )
      } else {
        loaded_pkgs <- character(0)
      }

      # Load objects
      if (length(objects_to_use) > 0) {
        loaded_objects <- evaluation_session$run(
          function(objects_to_use) {
            for (i in seq_along(objects_to_use)) {
              obj <- objects_to_use[[i]]
              obj_name <- names(objects_to_use)[i]
              assign(obj_name, obj, envir = .GlobalEnv)
            }
            ls(envir = .GlobalEnv)
          },
          args = list(objects_to_use = objects_to_use)
        )

        if (!all(names(objects_to_use) %in% loaded_objects)) {
          stop(paste0(
            "The following objects could not be loaded: ",
            paste(names(objects_to_use)[!(names(objects_to_use) %in% loaded_objects)], collapse = ", ")
          ))
        }
      } else {
        loaded_objects <- character(0)
      }
    }
  } else {
    evaluation_session <- NULL
    loaded_pkgs <- if (length(pkgs_to_use) == 0) character(0) else pkgs_to_use
    loaded_objects <- if (length(objects_to_use) == 0) character(0) else names(objects_to_use)
  }

  return(list(
    evaluation_session = evaluation_session,
    loaded_pkgs = loaded_pkgs,
    loaded_objects = loaded_objects,
    existing_session_info = existing_session_info
  ))
}

#' Get session state information with detailed object summaries
#'
#' @param session R session object
#' @return List with session state info
#' @noRd
get_session_info <- function(session) {
  if (is.null(session)) {
    return(NULL)
  }

  session_info <- tryCatch(
    {
      session$run(
        function() {
          # Get loaded packages
          session_info <- utils::sessionInfo()
          loaded_packages <- if (is.null(session_info$otherPkgs)) {
            character(0)
          } else {
            names(session_info$otherPkgs)
          }

          # Get objects in global environment
          all_objects <- ls(envir = .GlobalEnv)
          object_info <- if (length(all_objects) > 0) {
            sapply(all_objects, function(obj_name) {
              tryCatch(
                {
                  obj <- get(obj_name, envir = .GlobalEnv)
                  obj_class <- class(obj)[1]

                  # Add dimension information for different object types
                  if (is.matrix(obj) || is.data.frame(obj)) {
                    dims <- dim(obj)
                    if (prod(dims) > 1e6) { # Large objects (>1M elements)
                      dim_str <- paste0(" ", format(dims[1], big.mark = ","), "×", format(dims[2], big.mark = ","))
                    } else {
                      dim_str <- paste0(" ", dims[1], "×", dims[2])
                    }
                    paste0(obj_name, " (", obj_class, ":", dim_str, ")")
                  } else if (is.vector(obj) || is.list(obj)) {
                    len <- length(obj)
                    if (len > 1000) {
                      len_str <- format(len, big.mark = ",")
                    } else {
                      len_str <- as.character(len)
                    }
                    paste0(obj_name, " (", obj_class, ": length ", len_str, ")")
                  } else if (is.array(obj)) {
                    dims <- dim(obj)
                    dim_str <- paste(dims, collapse = "×")
                    paste0(obj_name, " (", obj_class, ": ", dim_str, ")")
                  } else {
                    paste0(obj_name, " (", obj_class, ")")
                  }
                },
                error = function(e) {
                  paste0(obj_name, " (unknown)")
                }
              )
            }, simplify = TRUE, USE.NAMES = FALSE)
          } else {
            character(0)
          }

          list(
            loaded_packages = loaded_packages,
            defined_objects = all_objects,
            object_summary = object_info
          )
        }
      )
    },
    error = function(e) {
      warning("Failed to get session info: ", e$message)
      list(
        loaded_packages = character(0),
        defined_objects = character(0),
        object_summary = character(0)
      )
    }
  )

  return(list(
    loaded_packages = session_info$loaded_packages,
    defined_objects = session_info$defined_objects,
    object_summary = session_info$object_summary
  ))
}

#' Modify prompt with additional context
#'
#' @param original_prompt_text Original prompt
#' @param add_text Additional instruction text (can be NULL)
#' @param list_packages Whether to list packages
#' @param list_objects Whether to list objects
#' @param loaded_pkgs Loaded packages
#' @param objects_to_use Objects to use
#' @param evaluate_code Whether evaluating code
#' @param return_mode Return mode
#' @param existing_session_info Existing session information
#' @return Modified prompt text
#' @noRd
modify_prompt <- function(original_prompt_text, add_text, list_packages, list_objects,
                          loaded_pkgs, objects_to_use, evaluate_code, return_mode, existing_session_info = NULL) {
  # Always include the default R code instruction
  default_instruction <- "You must code in the programming language 'R' to answer this prompt. Provide R code between ```r and ```."

  # Start with original prompt and default instruction
  new_text <- paste(original_prompt_text, "\n\n", default_instruction, sep = "")

  # Add user's additional text if provided
  if (!is.null(add_text)) {
    new_text <- paste(new_text, "\n\n", add_text, sep = "")
  }

  # Add existing session state information
  if (!is.null(existing_session_info)) {
    session_context <- "\n**IMPORTANT: You are continuing from an existing R session with the following state:**\n"

    if (length(existing_session_info$loaded_packages) > 0) {
      session_context <- paste0(
        session_context,
        "- Already loaded packages: ", paste(existing_session_info$loaded_packages, collapse = ", "), "\n"
      )
    }

    if (length(existing_session_info$defined_objects) > 0) {
      session_context <- paste0(
        session_context,
        "- Already defined objects: ", paste(existing_session_info$object_summary, collapse = ", "), "\n"
      )
    }

    session_context <- paste0(
      session_context,
      "\n**Please build upon the existing session state. Avoid redundant operations (re-loading packages, re-defining existing variables unless explicitly requested).**\n"
    )

    new_text <- paste(new_text, "\n", session_context, sep = "")
  }

  if (list_packages & length(loaded_pkgs) > 0) {
    new_text <- paste(new_text, "\n",
      "You can use functions from these packages: ",
      paste(loaded_pkgs, collapse = ", "), ".",
      sep = ""
    )
  }

  if (list_objects & length(objects_to_use) > 0) {
    objects_summary <- paste(names(objects_to_use),
      paste0("(", sapply(objects_to_use, function(obj) class(obj)[1]), ")"),
      sep = " ", collapse = ", "
    )

    new_text <- paste(new_text, "\n",
      "These objects already exist in the R session: ", objects_summary, ".\n\n",
      "Do not define these objects in your R code.",
      sep = ""
    )
  }

  if (evaluate_code & return_mode == "console") {
    new_text <- paste(new_text, "\n",
      "The R code should produce console output that answers the prompt.",
      sep = ""
    )
  }
  if (evaluate_code & return_mode == "object") {
    new_text <- paste(new_text, "\n",
      "The R code should produce an object that answers the prompt.",
      sep = ""
    )
  }

  return(new_text)
}

#' Retry R codes extraction with intelligent retry prompts
#'
#' @param chat_obj Chat object
#' @param extracted_code Current extracted codes
#' @param start_turn Starting turn index
#' @param max_iterations Maximum attempts
#' @return Extracted codes
#' @noRd
retry_r_extraction <- function(chat_obj, extracted_code, start_turn, max_iterations) {
  if (length(extracted_code) > 0) {
    return(extracted_code)
  }

  attempt <- 0

  while (length(extracted_code) == 0 && attempt < max_iterations) {
    attempt <- attempt + 1

    if (attempt == 1) {
      # First retry: use full history with context
      prompt_retry <- prompt_from_history(
        chat_obj,
        add_text = "No R code detected. You must provide R code between ```r and ```.",
        add_role = "user",
        start_turn_index = start_turn
      )
    } else {
      # Subsequent retries: use simple retry prompt without full history
      prompt_retry <- retry_code_prompt(
        chat_obj,
        error_message = "No R code detected in previous response",
        correction_instruction = "You must provide R code between ```r and ```."
      )
    }

    llm_response <- chat_obj$chat(prompt_retry)
    extracted_code <- extract_r_code(llm_response)

    if (length(extracted_code) > 0) {
      cat(sprintf(
        "R codes extracted successfully on attempt %d/%d\n",
        attempt, max_iterations
      ))
    }
  }

  if (length(extracted_code) == 0) {
    stop(sprintf(
      "LLM failed to provide a valid R code response after %d attempts",
      max_iterations
    ))
  }

  return(extracted_code)
}

#' Retry R codes parsing with intelligent retry prompts
#'
#' @param chat_obj Chat object
#' @param extracted_code Current extracted codes
#' @param start_turn Starting turn index
#' @param max_iterations Maximum attempts
#' @return List with extracted_code and parsed_code
#' @noRd
retry_r_parsing <- function(chat_obj, extracted_code, start_turn, max_iterations) {
  parsed_code <- tryCatch(parse(text = extracted_code), error = function(e) e)

  if (inherits(parsed_code, "error")) {
    attempt <- 0

    while (inherits(parsed_code, "error") && attempt < max_iterations) {
      attempt <- attempt + 1

      if (attempt == 1) {
        # First retry: use full history with context
        prompt_retry <- prompt_from_history(
          chat_obj,
          add_text = paste("Invalid R code detected:\n",
            "    ", parsed_code$message, "\n",
            "Please provide syntactically correct R code.",
            sep = ""
          ),
          add_role = "user",
          start_turn_index = start_turn
        )
      } else {
        # Subsequent retries: use simple retry prompt without full history
        prompt_retry <- retry_code_prompt(
          chat_obj,
          error_message = paste("R syntax error:", parsed_code$message),
          correction_instruction = "Please provide syntactically correct R code."
        )
      }

      llm_response <- chat_obj$chat(prompt_retry)
      extracted_code <- extract_r_code(llm_response)
      extracted_code <- retry_r_extraction(chat_obj, extracted_code, start_turn, max_iterations)

      parsed_code <- tryCatch(parse(text = extracted_code), error = function(e) e)
    }

    if (inherits(parsed_code, "error")) {
      stop(sprintf(
        "LLM failed to provide a valid R code parsing after %d attempts",
        max_iterations
      ))
    }
  }

  return(list(extracted_code = extracted_code, parsed_code = parsed_code))
}

#' Retry R codes running with intelligent retry prompts
#'
#' @param chat_obj Chat object
#' @param extracted_code Current extracted code
#' @param parsed_code Current parsed code
#' @param start_turn Starting turn index
#' @param max_iterations Maximum attempts
#' @param evaluation_session R session for evaluation
#' @param return_mode Return mode specification
#' @return List with extracted_code, parsed_code, output, and final_session
#' @noRd
retry_r_run <- function(chat_obj, extracted_code, parsed_code, start_turn, max_iterations, evaluation_session, return_mode) {
  clone_session <- evaluation_session$clone()
  output <- clone_session$run_with_output(
    function(r_code) {
      eval(parse(text = r_code), envir = .GlobalEnv)
    },
    args = list(r_code = extracted_code)
  )

  # Check if errors occurred during execution
  if (!is.null(output$error)) {
    attempt <- 0

    while (!is.null(output$error) && attempt < max_iterations) {
      attempt <- attempt + 1

      if (attempt == 1) {
        # First retry: use full history with context
        prompt_retry <- prompt_from_history(
          chat_obj,
          add_text = paste("An error occurred while executing the R code:\n",
            "    ", output$error, "\n",
            "Please provide corrected R code.",
            sep = ""
          ),
          add_role = "user",
          start_turn_index = start_turn
        )
      } else {
        # Subsequent retries: use simple retry prompt without full history
        prompt_retry <- retry_code_prompt(
          chat_obj,
          error_message = paste("R execution error:", output$error),
          correction_instruction = "Please provide corrected R code."
        )
      }

      llm_response <- chat_obj$chat(prompt_retry)
      new_extracted_code <- extract_r_code(llm_response)
      new_extracted_code <- retry_r_extraction(chat_obj, new_extracted_code, start_turn, max_iterations)

      parse_result <- retry_r_parsing(chat_obj, new_extracted_code, start_turn, max_iterations)
      extracted_code <- parse_result$extracted_code
      parsed_code <- parse_result$parsed_code

      clone_session <- evaluation_session$clone()
      output <- clone_session$run_with_output(
        function(r_code) {
          eval(parse(text = r_code), envir = .GlobalEnv)
        },
        args = list(r_code = extracted_code)
      )
    }

    if (!is.null(output$error)) {
      stop(sprintf(
        "LLM failed to provide working R code after %d attempts. Last error: %s",
        max_iterations, output$error
      ))
    }
  }

  # Check if the code produced any relevant output
  if (output$stdout == "" & return_mode == "console") {
    attempt <- 0

    while (output$stdout == "" && attempt < max_iterations) {
      attempt <- attempt + 1

      if (attempt == 1) {
        # First retry: use full history with context
        prompt_retry <- prompt_from_history(
          chat_obj,
          add_text = "The R code did not produce any console output. Please provide R code that produces console output.",
          add_role = "user",
          start_turn_index = start_turn
        )
      } else {
        # Subsequent retries: use simple retry prompt without full history
        prompt_retry <- retry_code_prompt(
          chat_obj,
          error_message = "No console output produced",
          correction_instruction = "Please provide R code that produces console output."
        )
      }

      llm_response <- chat_obj$chat(prompt_retry)
      new_extracted_code <- extract_r_code(llm_response)
      new_extracted_code <- retry_r_extraction(chat_obj, new_extracted_code, start_turn, max_iterations)

      parse_result <- retry_r_parsing(chat_obj, new_extracted_code, start_turn, max_iterations)
      extracted_code <- parse_result$extracted_code
      parsed_code <- parse_result$parsed_code

      clone_session <- evaluation_session$clone()
      output <- clone_session$run_with_output(
        function(r_code) {
          eval(parse(text = r_code), envir = .GlobalEnv)
        },
        args = list(r_code = extracted_code)
      )
    }

    if (output$stdout == "") {
      stop(sprintf(
        "LLM failed to provide R code with console output after %d attempts",
        max_iterations
      ))
    }
  }

  if (is.null(output$result) & return_mode == "object") {
    attempt <- 0

    while (is.null(output$result) && attempt < max_iterations) {
      attempt <- attempt + 1

      if (attempt == 1) {
        # First retry: use full history with context
        prompt_retry <- prompt_from_history(
          chat_obj,
          add_text = "The R code did not produce an object. Please provide R code that produces an object.",
          add_role = "user",
          start_turn_index = start_turn
        )
      } else {
        # Subsequent retries: use simple retry prompt without full history
        prompt_retry <- retry_code_prompt(
          chat_obj,
          error_message = "No object result produced",
          correction_instruction = "Please provide R code that produces an object."
        )
      }

      llm_response <- chat_obj$chat(prompt_retry)
      new_extracted_code <- extract_r_code(llm_response)
      new_extracted_code <- retry_r_extraction(chat_obj, new_extracted_code, start_turn, max_iterations)

      parse_result <- retry_r_parsing(chat_obj, new_extracted_code, start_turn, max_iterations)
      extracted_code <- parse_result$extracted_code
      parsed_code <- parse_result$parsed_code

      clone_session <- evaluation_session$clone()
      output <- clone_session$run_with_output(
        function(r_code) {
          eval(parse(text = r_code), envir = .GlobalEnv)
        },
        args = list(r_code = extracted_code)
      )
    }

    if (is.null(output$result)) {
      stop(sprintf(
        "LLM failed to provide R code that produces an object after %d attempts",
        max_iterations
      ))
    }
  }

  if (is.null(output$stdout) & is.null(output$result)) {
    stop("The R code did not produce any output.")
  }

  return(list(
    extracted_code = extracted_code,
    parsed_code = parsed_code,
    output = output,
    final_session = clone_session
  ))
}

#' Build retry prompt for code correction
#' @param chat_obj Chat object
#' @param error_message Error message from execution
#' @param correction_instruction Specific correction instruction
#' @return Concise retry prompt
#' @noRd
retry_code_prompt <- function(chat_obj, error_message, correction_instruction) {
  turns <- chat_obj$get_turns()
  if (length(turns) > 0) {
    last_assistant <- tail(turns[sapply(turns, function(x) x@role == "assistant")], 1)
    if (length(last_assistant) > 0) {
      # Find the corresponding user question
      assistant_turns <- which(sapply(turns, function(x) x@role == "assistant"))
      user_turns <- which(sapply(turns, function(x) x@role == "user"))
      last_assistant_idx <- tail(assistant_turns, 1)
      corresponding_user_idx <- tail(user_turns[user_turns < last_assistant_idx], 1)

      # Build context with user question if available
      context <- ""
      if (length(corresponding_user_idx) > 0) {
        context <- paste("Original request:\n", turns[[corresponding_user_idx]]@text, "\n\n")
      }

      last_code <- extract_r_code(last_assistant[[1]]@text)
      if (length(last_code) > 0) {
        context <- paste0(
          context, "Previous code attempt:\n```r\n",
          paste(last_code, collapse = "\n"),
          "\n```\n\n"
        )
      }
    } else {
      context <- ""
    }
  } else {
    context <- ""
  }

  # 构建简洁的重试prompt
  retry_prompt <- paste0(
    context,
    "Error occurred:\n", error_message, "\n\n",
    correction_instruction
  )

  return(retry_prompt)
}
