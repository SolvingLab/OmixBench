#' Create OpenAI-compatible LLM provider with enhanced error handling
#'
#' This function creates an OpenAI-compatible LLM provider with comprehensive
#' error handling and testing capabilities. It automatically handles max_tokens
#' limits by falling back to the model's maximum when exceeded.
#'
#' @param base_url The base URL for the OpenAI-compatible API
#' @param api_key The API key for authentication. If NULL, will use LLM_API_KEY env var
#' @param model The model name to use
#' @param temperature The temperature parameter for response randomness
#' @param max_tokens Maximum number of tokens in response (will auto-adjust if exceeds model limit)
#' @param timeout Request timeout in seconds
#' @param stream Whether to use streaming responses
#' @param verbose Whether to show verbose output
#' @param skip_test Whether to skip the availability test (useful for problematic providers)
#' @param test_mode The testing mode: "full", "http_only", "skip"
#' @param ... Additional parameters to pass to the model
#' @return A configured LLM provider object
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
llm_provider <- function(
    base_url = "https://api.openai.com/v1/chat/completions",
    api_key = NULL,
    model = "gpt-4o-mini",
    temperature = 0.2,
    max_tokens = 5000,
    timeout = 100,
    stream = FALSE,
    verbose = TRUE,
    skip_test = FALSE,
    test_mode = c("full", "http_only", "skip"),
    ...) {
  test_mode <- match.arg(test_mode)

  # Handle API key
  if (is.null(api_key)) {
    api_key <- Sys.getenv("LLM_API_KEY")
    if (api_key == "") {
      stop("API key not provided. Please set api_key parameter or LLM_API_KEY environment variable.")
    }
  }

  # Store original max_tokens for comparison
  original_max_tokens <- max_tokens
  adjusted_max_tokens <- max_tokens

  # Create the model provider using tidyprompt's architecture
  model_provider <- create_provider_internal(
    base_url, api_key, model, temperature, stream,
    timeout, adjusted_max_tokens, verbose, ...
  )

  # Skip all tests if requested
  if (skip_test || test_mode == "skip") {
    cli::cli_alert_info(paste("Skipping availability test for model:", model))
    cli::cli_alert_success(paste("Model provider created for", model))
    return(model_provider)
  }

  # Perform testing based on test_mode with max_tokens fallback
  test_result <- NULL
  max_tokens_adjusted <- FALSE

  if (test_mode == "http_only") {
    test_result <- test_http_connection_with_fallback(
      base_url, api_key, model, adjusted_max_tokens, verbose
    )

    # Check if max_tokens was adjusted
    if (!is.null(test_result$adjusted_max_tokens)) {
      adjusted_max_tokens <- test_result$adjusted_max_tokens
      max_tokens_adjusted <- TRUE

      # Recreate provider with adjusted max_tokens
      model_provider <- create_provider_internal(
        base_url, api_key, model, temperature, stream,
        timeout, adjusted_max_tokens, verbose, ...
      )
    }

    success <- test_result$success
  } else {
    # Full tidyprompt compatibility test
    test_result <- test_full_compatibility_with_fallback(
      model_provider, base_url, api_key, model, temperature,
      stream, timeout, adjusted_max_tokens, verbose, ...
    )

    if (!is.null(test_result$adjusted_max_tokens)) {
      adjusted_max_tokens <- test_result$adjusted_max_tokens
      max_tokens_adjusted <- TRUE
      model_provider <- test_result$provider
    }

    success <- test_result$success
  }

  # Show warning if max_tokens was adjusted
  if (max_tokens_adjusted) {
    cli::cli_alert_warning(
      paste0(
        "max_tokens adjusted from ", original_max_tokens,
        " to ", adjusted_max_tokens, " (model limit)"
      )
    )
  }

  if (success) {
    cli::cli_alert_success(paste("Model", model, "is ready to use"))
  } else {
    cli::cli_alert_warning(paste("Model", model, "may have compatibility issues"))
    cli::cli_alert_info("You can still use this provider, but consider using skip_test = TRUE for future calls")
  }

  return(model_provider)
}

#' Internal helper to create provider
#' @keywords internal
create_provider_internal <- function(base_url, api_key, model, temperature,
                                     stream, timeout, max_tokens, verbose, ...) {
  tidyprompt::llm_provider_openai(
    parameters = list(
      model = model,
      temperature = temperature,
      stream = stream,
      timeout = timeout,
      max_tokens = max_tokens,
      ...
    ),
    url = base_url,
    api_key = api_key,
    verbose = verbose
  )
}

#' Parse max_tokens from error message
#' @keywords internal
parse_max_tokens_from_error <- function(error_message) {
  # Try to extract max value from patterns like:
  # "the valid range of max_tokens is [1, 8192]"
  # "max_tokens must be between 1 and 8192"
  # "maximum tokens: 8192"

  patterns <- list(
    "\\[\\d+,\\s*(\\d+)\\]", # [1, 8192]
    "between \\d+ and (\\d+)", # between 1 and 8192
    "maximum.*?(\\d+)", # maximum tokens: 8192
    "max.*?is (\\d+)", # max is 8192
    "up to (\\d+)", # up to 8192
    "range.*?(\\d+)$" # range ... 8192
  )

  for (pattern in patterns) {
    match <- regmatches(error_message, regexec(pattern, error_message, ignore.case = TRUE))
    if (length(match[[1]]) > 1) {
      max_val <- as.integer(match[[1]][2])
      if (!is.na(max_val) && max_val > 0) {
        return(max_val)
      }
    }
  }

  return(NULL)
}

#' Test HTTP connection with max_tokens fallback
#' @keywords internal
test_http_connection_with_fallback <- function(base_url, api_key, model,
                                               max_tokens, verbose = TRUE) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    cli::cli_alert_danger("httr package required for HTTP testing")
    return(list(success = FALSE, adjusted_max_tokens = NULL))
  }

  if (verbose) cli::cli_alert_info("Testing HTTP connection...")

  # Prepare minimal test request
  body <- list(
    model = model,
    messages = list(list(role = "user", content = "Hi")),
    max_tokens = max_tokens,
    temperature = 0
  )

  # Make HTTP request
  response <- tryCatch(
    {
      httr::POST(
        url = base_url,
        httr::add_headers(
          "Content-Type" = "application/json",
          "Authorization" = paste("Bearer", api_key)
        ),
        body = jsonlite::toJSON(body, auto_unbox = TRUE),
        httr::timeout(30)
      )
    },
    error = function(e) {
      if (verbose) cli::cli_alert_danger(paste("HTTP request error:", e$message))
      return(NULL)
    }
  )

  if (is.null(response)) {
    return(list(success = FALSE, adjusted_max_tokens = NULL))
  }

  status_code <- httr::status_code(response)

  if (status_code == 200) {
    if (verbose) cli::cli_alert_success("HTTP test successful")
    return(list(success = TRUE, adjusted_max_tokens = NULL))
  } else if (status_code == 400) {
    # Check if it's a max_tokens error
    content <- tryCatch(
      httr::content(response, "text", encoding = "UTF-8"),
      error = function(e) NULL
    )

    if (!is.null(content) && grepl("max_tokens", content, ignore.case = TRUE)) {
      if (verbose) cli::cli_alert_info("Detected max_tokens limit error, attempting to parse...")

      # Try to parse the model's max_tokens limit
      model_max_tokens <- parse_max_tokens_from_error(content)

      if (!is.null(model_max_tokens)) {
        if (verbose) {
          cli::cli_alert_info(paste("Detected model max_tokens limit:", model_max_tokens))
        }

        # Retry with adjusted max_tokens
        body$max_tokens <- model_max_tokens

        response_retry <- tryCatch(
          {
            httr::POST(
              url = base_url,
              httr::add_headers(
                "Content-Type" = "application/json",
                "Authorization" = paste("Bearer", api_key)
              ),
              body = jsonlite::toJSON(body, auto_unbox = TRUE),
              httr::timeout(30)
            )
          },
          error = function(e) NULL
        )

        if (!is.null(response_retry) && httr::status_code(response_retry) == 200) {
          if (verbose) cli::cli_alert_success("HTTP test successful with adjusted max_tokens")
          return(list(success = TRUE, adjusted_max_tokens = model_max_tokens))
        }
      }
    }

    if (verbose) {
      cli::cli_alert_warning(paste("HTTP test returned status 400"))
      if (!is.null(content)) {
        cli::cli_alert_info(paste("Response:", substr(content, 1, 200)))
      }
    }
    return(list(success = FALSE, adjusted_max_tokens = NULL))
  } else if (status_code == 401) {
    if (verbose) cli::cli_alert_danger("Authentication failed - check API key")
    return(list(success = FALSE, adjusted_max_tokens = NULL))
  } else {
    if (verbose) {
      cli::cli_alert_warning(paste("HTTP test returned status:", status_code))
      content <- tryCatch(httr::content(response, "text"), error = function(e) NULL)
      if (!is.null(content)) {
        cli::cli_alert_info(paste("Response:", substr(content, 1, 200)))
      }
    }
    return(list(success = FALSE, adjusted_max_tokens = NULL))
  }
}

#' Test full compatibility with tidyprompt and max_tokens fallback
#' @keywords internal
test_full_compatibility_with_fallback <- function(model_provider, base_url, api_key,
                                                  model, temperature, stream, timeout,
                                                  max_tokens, verbose, ...) {
  if (verbose) cli::cli_alert_info("Testing tidyprompt compatibility...")

  result <- tryCatch(
    {
      # Test with a very simple prompt to minimize data parsing issues
      result <- tidyprompt::send_prompt("Hi", model_provider, verbose = FALSE)

      # Check if we got a valid result
      if (!is.null(result) && length(result) > 0) {
        if (verbose) cli::cli_alert_success("tidyprompt compatibility test passed")
        return(list(success = TRUE, adjusted_max_tokens = NULL, provider = model_provider))
      } else {
        if (verbose) cli::cli_alert_warning("tidyprompt test returned empty result")
        return(list(success = FALSE, adjusted_max_tokens = NULL, provider = model_provider))
      }
    },
    error = function(e) {
      error_msg <- e$message

      # Check if error is related to max_tokens
      if (grepl("max_tokens", error_msg, ignore.case = TRUE)) {
        if (verbose) cli::cli_alert_info("Detected max_tokens error in tidyprompt test")

        # Try to parse max_tokens from error
        model_max_tokens <- parse_max_tokens_from_error(error_msg)

        if (!is.null(model_max_tokens)) {
          if (verbose) {
            cli::cli_alert_info(paste("Detected model max_tokens limit:", model_max_tokens))
            cli::cli_alert_info("Recreating provider with adjusted max_tokens...")
          }

          # Create new provider with adjusted max_tokens
          new_provider <- create_provider_internal(
            base_url, api_key, model, temperature, stream,
            timeout, model_max_tokens, verbose, ...
          )

          # Retry the test
          retry_result <- tryCatch(
            {
              result <- tidyprompt::send_prompt("Hi", new_provider, verbose = FALSE)
              if (!is.null(result) && length(result) > 0) {
                if (verbose) cli::cli_alert_success("tidyprompt test passed with adjusted max_tokens")
                return(list(success = TRUE, adjusted_max_tokens = model_max_tokens, provider = new_provider))
              } else {
                return(list(success = FALSE, adjusted_max_tokens = model_max_tokens, provider = new_provider))
              }
            },
            error = function(e2) {
              if (verbose) cli::cli_alert_warning("tidyprompt test still failed after adjustment")
              return(list(success = FALSE, adjusted_max_tokens = model_max_tokens, provider = new_provider))
            }
          )

          return(retry_result)
        }
      }

      # Other errors (not max_tokens related)
      if (verbose) {
        cli::cli_alert_warning("tidyprompt compatibility test failed")
        cli::cli_alert_info(paste("Error:", error_msg))

        if (grepl("differing number of rows", error_msg)) {
          cli::cli_alert_info("This may be due to data parsing issues in tidyprompt")
          cli::cli_alert_info("The API itself may work fine - consider using skip_test = TRUE")
        }
      }

      return(list(success = FALSE, adjusted_max_tokens = NULL, provider = model_provider))
    }
  )

  return(result)
}



#' Create Ollama LLM provider with enhanced availability check and auto-download
#'
#' This function creates an Ollama LLM provider with better error handling
#' and follows tidyprompt best practices.
#'
#' @param base_url The base URL for the Ollama API
#' @param model The model name to use
#' @param temperature The temperature parameter for response randomness
#' @param max_tokens Maximum number of tokens in response
#' @param timeout Request timeout in seconds
#' @param stream Whether to use streaming responses
#' @param verbose Whether to show verbose output
#' @param skip_test Whether to skip the availability test
#' @param auto_download Whether to automatically download missing models
#' @param ... Additional parameters to pass to the model
#' @return A configured LLM provider object
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
llm_ollama <- function(
    base_url = "http://localhost:11434/api/chat",
    model = "qwen2.5:1.5b-instruct",
    temperature = 0.2,
    max_tokens = 5000,
    timeout = 100,
    stream = TRUE,
    verbose = TRUE,
    skip_test = FALSE,
    auto_download = TRUE,
    ...) {
  # Create the model provider using tidyprompt's architecture
  model_provider <- tidyprompt::llm_provider_ollama(
    parameters = list(
      model = model,
      temperature = temperature,
      timeout = timeout,
      max_tokens = max_tokens,
      stream = stream,
      ...
    ),
    url = base_url,
    verbose = verbose
  )

  # Skip test if requested
  if (skip_test) {
    cli::cli_alert_info(paste("Skipping availability test for model:", model))
    cli::cli_alert_success(paste("Model provider created for", model))
    return(model_provider)
  }

  # Extract server URL for model management
  server_url <- gsub("/api/chat$", "", base_url)

  # Test model availability with enhanced error handling
  test_result <- tryCatch(
    {
      # Test with tidyprompt
      result <- tidyprompt::send_prompt("Hi", model_provider, verbose = FALSE)
      if (!is.null(result) && length(result) > 0) {
        TRUE
      } else {
        FALSE
      }
    },
    error = function(e) {
      if (verbose) {
        cli::cli_alert_warning(paste("Model", model, "test failed:", e$message))
      }

      # If auto_download is enabled, try to download the model
      if (auto_download) {
        if (verbose) cli::cli_alert_info(paste("Attempting to download model:", model))

        download_result <- tryCatch(
          {
            ollama_download_model(model, server_url)
            TRUE
          },
          error = function(download_error) {
            if (verbose) cli::cli_alert_danger(paste("Failed to download model:", download_error$message))
            FALSE
          }
        )

        if (download_result) {
          # Test again after download
          tryCatch(
            {
              result <- tidyprompt::send_prompt("Hi", model_provider, verbose = FALSE)
              if (!is.null(result) && length(result) > 0) {
                if (verbose) cli::cli_alert_success(paste("Model", model, "downloaded and ready"))
                return(TRUE)
              } else {
                return(FALSE)
              }
            },
            error = function(e2) {
              if (verbose) cli::cli_alert_danger(paste("Model still not working after download:", e2$message))
              return(FALSE)
            }
          )
        } else {
          return(FALSE)
        }
      } else {
        return(FALSE)
      }
    }
  )

  if (test_result) {
    cli::cli_alert_success(paste("Model", model, "is ready to use"))
  } else {
    cli::cli_alert_warning(paste("Model", model, "may not be available"))
    cli::cli_alert_info("You can still try to use this provider")
  }

  return(model_provider)
}

#' Comprehensive LLM connection diagnostics
#'
#' This function provides detailed diagnostics for LLM connection issues,
#' helping identify problems at different levels of the stack.
#'
#' @param base_url The API base URL
#' @param api_key The API key
#' @param model The model name
#' @param test_tidyprompt Whether to test tidyprompt compatibility
#' @export
diagnose_llm_connection <- function(base_url, api_key, model, test_tidyprompt = TRUE) {
  cli::cli_h1("LLM Connection Comprehensive Diagnosis")
  cli::cli_alert_info(paste("Testing connection to:", base_url))
  cli::cli_alert_info(paste("Model:", model))
  cli::cli_alert_info(paste("API key preview:", paste0(substr(api_key, 1, 8), "...")))

  results <- list()

  # Test 1: Basic connectivity
  cli::cli_h2("Test 1: Basic Network Connectivity")
  results$connectivity <- tryCatch(
    {
      # Try to reach the base domain
      base_domain <- gsub("^https?://([^/]+).*", "\\1", base_url)
      response <- httr::GET(paste0("https://", base_domain), httr::timeout(10))
      status <- httr::status_code(response)
      cli::cli_alert_success(paste("Network connectivity OK, status:", status))
      TRUE
    },
    error = function(e) {
      cli::cli_alert_danger(paste("Network connectivity failed:", e$message))
      FALSE
    }
  )

  # Test 2: API endpoint accessibility
  cli::cli_h2("Test 2: API Endpoint Accessibility")
  results$endpoint <- tryCatch(
    {
      # Try to reach the API endpoint
      response <- httr::GET(base_url, httr::timeout(10))
      status <- httr::status_code(response)
      if (status == 405) { # Method not allowed is expected for GET on POST endpoint
        cli::cli_alert_success("API endpoint accessible (405 Method Not Allowed is expected)")
        TRUE
      } else {
        cli::cli_alert_info(paste("API endpoint returned status:", status))
        TRUE
      }
    },
    error = function(e) {
      cli::cli_alert_danger(paste("API endpoint not accessible:", e$message))
      FALSE
    }
  )

  # Test 3: Authentication and model availability
  cli::cli_h2("Test 3: Authentication and Model Test")
  results$auth <- test_http_connection(base_url, api_key, model, verbose = TRUE)

  # Test 4: tidyprompt compatibility (if requested)
  if (test_tidyprompt) {
    cli::cli_h2("Test 4: tidyprompt Compatibility")
    results$tidyprompt <- tryCatch(
      {
        provider <- tidyprompt::llm_provider_openai(
          parameters = list(model = model, temperature = 0.1, max_tokens = 5),
          url = base_url,
          api_key = api_key,
          verbose = FALSE
        )

        result <- tidyprompt::send_prompt("Hi", provider, verbose = FALSE)
        if (!is.null(result) && length(result) > 0) {
          cli::cli_alert_success("tidyprompt compatibility: PASSED")
          True
        } else {
          cli::cli_alert_warning("tidyprompt compatibility: FAILED (empty result)")
          FALSE
        }
      },
      error = function(e) {
        cli::cli_alert_danger(paste("tidyprompt compatibility: FAILED -", e$message))
        if (grepl("differing number of rows", e$message)) {
          cli::cli_alert_info("This is a known data parsing issue in tidyprompt")
          cli::cli_alert_info("The API likely works fine - use skip_test = TRUE")
        }
        FALSE
      }
    )
  }

  # Summary and recommendations
  cli::cli_h2("Summary and Recommendations")

  if (all(unlist(results[1:3]))) {
    if (!test_tidyprompt || results$tidyprompt) {
      cli::cli_alert_success("All tests passed! Your connection should work perfectly.")
    } else {
      cli::cli_alert_warning("API works but tidyprompt has parsing issues")
      cli::cli_alert_info("Recommended solution:")
      cli::cli_code("llm_client <- llm_openai(..., skip_test = TRUE)")
    }
  } else {
    cli::cli_alert_danger("Some tests failed. Issues found:")
    if (!results$connectivity) cli::cli_alert("Network connectivity problems")
    if (!results$endpoint) cli::cli_alert("API endpoint not accessible")
    if (!results$auth) cli::cli_alert("Authentication or model issues")

    cli::cli_alert_info("Check your network, API key, and model name")
  }

  invisible(results)
}

#' List available models from Ollama server
#'
#' @param server_url The Ollama server URL (without /api/chat)
#' @return A data frame of available models or NULL if failed
ollama_list_models <- function(server_url = "http://localhost:11434") {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("httr package required for Ollama operations")
  }

  tryCatch(
    {
      response <- httr::GET(paste0(server_url, "/api/tags"), httr::timeout(10))

      if (httr::status_code(response) == 200) {
        content <- httr::content(response, "parsed")
        if (!is.null(content$models) && length(content$models) > 0) {
          models_df <- data.frame(
            name = sapply(content$models, function(x) x$name),
            size = sapply(content$models, function(x) x$size),
            modified = sapply(content$models, function(x) x$modified_at),
            stringsAsFactors = FALSE
          )
          return(models_df)
        } else {
          return(data.frame(name = character(0), size = numeric(0), modified = character(0)))
        }
      } else {
        cli::cli_alert_danger(paste("Failed to list models. Status:", httr::status_code(response)))
        return(NULL)
      }
    },
    error = function(e) {
      cli::cli_alert_danger(paste("Error listing models:", e$message))
      return(NULL)
    }
  )
}

#' Download a model from Ollama
#'
#' @param model The model name to download
#' @param server_url The Ollama server URL (without /api/chat)
#' @return TRUE if successful, FALSE otherwise
ollama_download_model <- function(model, server_url = "http://localhost:11434") {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("httr package required for Ollama operations")
  }

  cli::cli_alert_info(paste("Downloading model:", model))
  cli::cli_alert_info("This may take several minutes depending on model size...")

  tryCatch(
    {
      body <- list(name = model)

      response <- httr::POST(
        url = paste0(server_url, "/api/pull"),
        body = jsonlite::toJSON(body, auto_unbox = TRUE),
        httr::content_type_json(),
        httr::timeout(600) # 10 minutes timeout for download
      )

      if (httr::status_code(response) == 200) {
        cli::cli_alert_success(paste("Model", model, "downloaded successfully"))
        return(TRUE)
      } else {
        cli::cli_alert_danger(paste("Failed to download model. Status:", httr::status_code(response)))
        content <- tryCatch(httr::content(response, "text"), error = function(e) "Unknown error")
        cli::cli_alert_info(paste("Response:", content))
        return(FALSE)
      }
    },
    error = function(e) {
      cli::cli_alert_danger(paste("Error downloading model:", e$message))
      return(FALSE)
    }
  )
}

#' Check if Ollama server is running
#'
#' @param server_url The Ollama server URL (without /api/chat)
#' @return TRUE if server is running, FALSE otherwise
ollama_check_server <- function(server_url = "http://localhost:11434") {
  if (!requireNamespace("httr", quietly = TRUE)) {
    return(FALSE)
  }

  tryCatch(
    {
      response <- httr::GET(paste0(server_url, "/api/tags"), httr::timeout(5))
      return(httr::status_code(response) == 200)
    },
    error = function(e) {
      return(FALSE)
    }
  )
}

#' Get Ollama server version information
#'
#' @param server_url The Ollama server URL (without /api/chat)
#' @return Version information or NULL if failed
ollama_get_version <- function(server_url = "http://localhost:11434") {
  if (!requireNamespace("httr", quietly = TRUE)) {
    return(NULL)
  }

  tryCatch(
    {
      response <- httr::GET(paste0(server_url, "/api/version"), httr::timeout(5))

      if (httr::status_code(response) == 200) {
        content <- httr::content(response, "parsed")
        return(content)
      } else {
        return(NULL)
      }
    },
    error = function(e) {
      return(NULL)
    }
  )
}
