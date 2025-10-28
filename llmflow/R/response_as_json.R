#' Get JSON response from LLM with validation and retry
#'
#' @param chat_obj Chat object (LLM client). Must be a properly initialized LLM client instance.
#' @param prompt Character string. The user prompt to send to the LLM requesting JSON output.
#' @param schema List or NULL. Optional JSON schema for response validation (as R list structure).
#'   When provided, validates the LLM response against this schema. Default is NULL.
#' @param schema_strict Logical. Whether to use strict schema validation (no additional properties allowed).
#'   Only applies when schema is provided. Default is FALSE.
#' @param max_iterations Integer. Maximum number of retry attempts for invalid JSON or schema validation failures.
#'   Must be positive. Default is 3.
#' @return List. Parsed JSON response from the LLM.
#' @examples
#' \dontrun{
#' # Basic usage without schema
#' result <- response_as_json(
#'   chat_obj = llm_client,
#'   prompt = "List three colors"
#' )
#'
#' # With schema validation
#' schema <- list(
#'   type = "object",
#'   properties = list(
#'     equation = list(type = "string"),
#'     solution = list(type = "number")
#'   ),
#'   required = c("equation", "solution")
#' )
#' result <- response_as_json(
#'   chat_obj = llm_client,
#'   prompt = "How can I solve 8x + 7 = -23?",
#'   schema = schema,
#'   schema_strict = TRUE,
#'   max_iterations = 3
#' )
#' }
#'
#' @export
response_as_json <- function(chat_obj,
                             prompt,
                             schema = NULL,
                             schema_strict = FALSE,
                             max_iterations = 3) {
  # Check jsonvalidate package
  if (!requireNamespace("jsonvalidate", quietly = TRUE)) {
    utils::install.packages("jsonvalidate")
  }

  # Validate schema argument
  if (!is.null(schema)) {
    if (!is.list(schema)) {
      stop("The 'schema' argument must be a list (R object) representing a JSON schema")
    }
    if (length(schema) == 0) {
      stop("The 'schema' argument must be a non-empty list")
    }
  }

  # Record starting turn for retry context
  start_turn <- length(chat_obj$get_turns()) + 1

  # Build initial prompt with JSON instruction
  prompt_text <- build_json_prompt(prompt, schema)

  # Get initial response
  llm_response <- chat_obj$chat(prompt_text)
  jsons <- extract_json(llm_response)

  # Retry if no JSON found
  jsons <- retry_json_extraction(
    chat_obj,
    jsons,
    start_turn,
    max_iterations
  )

  # If schema provided, validate and retry if needed
  if (!is.null(schema)) {
    jsons <- validate_and_retry(
      chat_obj,
      jsons,
      schema,
      schema_strict,
      start_turn,
      max_iterations
    )
  }

  # Return single JSON if only one extracted
  if (length(jsons) == 1) {
    return(jsons[[1]])
  }

  return(jsons)
}

#' Build JSON prompt with schema instruction
#'
#' @param prompt Original prompt
#' @param schema Optional JSON schema
#' @return Formatted prompt string
#' @noRd
build_json_prompt <- function(prompt, schema = NULL) {
  base_prompt <- paste0(
    prompt, "\n\n",
    "You must format your response as a JSON object."
  )

  if (is.null(schema)) {
    schema_instruction <- "Your JSON object should match this example JSON object:\n{}"
  } else {
    schema_instruction <- paste0(
      "Your JSON object should match this JSON schema:\n",
      jsonlite::toJSON(schema, auto_unbox = TRUE, pretty = TRUE)
    )
  }

  return(paste0(base_prompt, "\n\n", schema_instruction))
}

#' Retry JSON extraction with iteration limit
#'
#' @param chat_obj Chat object
#' @param jsons Current extracted JSONs
#' @param start_turn Starting turn index
#' @param max_iterations Maximum attempts
#' @return Extracted JSONs
#' @noRd
retry_json_extraction <- function(chat_obj, jsons, start_turn, max_iterations) {
  if (length(jsons) > 0) {
    return(jsons)
  }

  attempt <- 0

  while (length(jsons) == 0 && attempt < max_iterations) {
    attempt <- attempt + 1
    if (attempt != 1) {
      start_turn <- length(chat_obj$get_turns()) - 1
    }

    prompt_retry <- prompt_from_history(
      chat_obj,
      add_text = "You must respond as a valid JSON object.",
      add_role = "user",
      start_turn_index = start_turn
    )

    llm_response <- chat_obj$chat(prompt_retry)
    jsons <- extract_json(llm_response)

    if (length(jsons) > 0) {
      cat(sprintf(
        "JSON extracted successfully on attempt %d/%d\n",
        attempt, max_iterations
      ))
    }
  }

  if (length(jsons) == 0) {
    stop(sprintf(
      "LLM failed to provide a valid JSON response after %d attempts",
      max_iterations
    ))
  }

  return(jsons)
}

#' Validate JSON against schema and retry if needed
#'
#' @param chat_obj Chat object
#' @param jsons Extracted JSONs
#' @param schema JSON schema
#' @param schema_strict Strict validation flag
#' @param start_turn Starting turn index
#' @param max_iterations Maximum attempts
#' @return Validated JSONs
#' @noRd
validate_and_retry <- function(chat_obj,
                               jsons,
                               schema,
                               schema_strict,
                               start_turn,
                               max_iterations) {
  # Prepare for validation
  if (length(jsons) == 1) {
    json_to_validate <- jsons[[1]]
  } else {
    json_to_validate <- jsons
  }

  answer_json <- jsonlite::toJSON(json_to_validate, auto_unbox = TRUE, pretty = TRUE)
  schema_json <- jsonlite::toJSON(schema, auto_unbox = TRUE, pretty = TRUE)

  # Initial validation
  validation_result <- jsonvalidate::json_validate(
    answer_json,
    schema_json,
    strict = schema_strict,
    verbose = TRUE
  )

  if (validation_result) {
    return(jsons)
  }

  # Retry with error feedback
  attempt <- 0

  while (!validation_result && attempt < max_iterations) {
    attempt <- attempt + 1
    if (attempt != 1) {
      start_turn <- length(chat_obj$get_turns()) - 1
    }

    # Get error details
    error_details <- attr(validation_result, "errors")

    # Build retry prompt with error feedback
    schema_instruction <- paste0(
      "Your JSON object should match this JSON schema:\n",
      jsonlite::toJSON(schema, auto_unbox = TRUE, pretty = TRUE)
    )

    prompt_retry <- prompt_from_history(
      chat_obj,
      add_text = paste0(
        "Your response did not match the expected JSON schema.\n\n",
        "Validation errors:\n",
        df_to_string(error_details),
        "\n\n",
        schema_instruction
      ),
      add_role = "user",
      start_turn_index = start_turn
    )

    llm_response <- chat_obj$chat(prompt_retry)
    jsons <- extract_json(llm_response)

    # Ensure we have JSON
    jsons <- retry_json_extraction(chat_obj, jsons, start_turn, max_iterations)

    # Prepare for revalidation
    if (length(jsons) == 1) {
      json_to_validate <- jsons[[1]]
    } else {
      json_to_validate <- jsons
    }

    answer_json <- jsonlite::toJSON(json_to_validate, auto_unbox = TRUE, pretty = TRUE)

    # Revalidate
    validation_result <- jsonvalidate::json_validate(
      answer_json,
      schema_json,
      strict = schema_strict,
      verbose = TRUE
    )

    if (validation_result) {
      cat(sprintf(
        "Schema validation passed on attempt %d/%d\n",
        attempt, max_iterations
      ))
    }
  }

  if (!validation_result) {
    stop(sprintf(
      "LLM failed to provide a schema-compliant JSON response after %d attempts",
      max_iterations
    ))
  }

  return(jsons)
}

#' Helper function to convert data frame to string
#'
#' @param df Data frame
#' @return String representation
#' @noRd
df_to_string <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return("No specific errors available")
  }

  paste(capture.output(print(df)), collapse = "\n")
}
