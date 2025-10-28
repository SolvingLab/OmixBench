#' Interactive JSON Schema Generator using tidyprompt
#'
#' This function creates an interactive system to generate JSON schemas based on user
#' descriptions. It supports multi-turn conversations until the user is satisfied
#' with the generated schema.
#'
#' @param description Initial description of the desired JSON structure
#' @param llm_client The LLM provider object (from llm_openai or llm_ollama)
#' @param max_iterations Maximum number of refinement iterations (default: 5)
#' @param interactive Whether to run in interactive mode (default: TRUE)
#' @param verbose Whether to show detailed conversation logs (default: TRUE)
#' @return A list containing the final JSON schema and conversation history
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
generate_json_schema <- function(description,
                                 llm_client,
                                 max_iterations = 5,
                                 interactive = TRUE,
                                 verbose = TRUE) {
  # Initialize conversation state
  conversation_state <- list(
    iteration = 1,
    description = description,
    current_schema = NULL,
    user_feedback = NULL,
    history = list(),
    satisfied = FALSE
  )

  # Initial schema generation prompt
  initial_prompt <- create_schema_prompt(conversation_state)

  if (verbose) {
    cli::cli_h2("Starting JSON Schema Generation")
    cli::cli_alert_info(paste("Description:", description))
    cli::cli_rule()
  }

  # Main conversation loop
  while (!conversation_state$satisfied && conversation_state$iteration <= max_iterations) {
    if (verbose) {
      cli::cli_h3(paste("Iteration", conversation_state$iteration))
    }

    # Generate or refine schema
    result <- generate_schema_step(initial_prompt, llm_client, conversation_state, verbose)

    # Update conversation state
    conversation_state$current_schema <- result$schema
    conversation_state$history <- append(
      conversation_state$history,
      list(list(
        iteration = conversation_state$iteration,
        prompt = result$prompt_used,
        schema = result$schema,
        timestamp = Sys.time()
      ))
    )

    if (verbose) {
      cli::cli_alert_success("Schema generated successfully!")
      print_schema_preview(result$schema)
    }

    # Get user feedback
    if (interactive) {
      conversation_state <- get_user_feedback(conversation_state, verbose)
    } else {
      conversation_state$satisfied <- TRUE
    }

    # Prepare next iteration
    if (!conversation_state$satisfied) {
      initial_prompt <- create_refinement_prompt(conversation_state)
      conversation_state$iteration <- conversation_state$iteration + 1
    }
  }

  # Final result
  if (conversation_state$satisfied) {
    if (verbose) {
      cli::cli_alert_success("Schema generation completed successfully!")
    }
  } else {
    if (verbose) {
      cli::cli_alert_warning("Maximum iterations reached. You can continue with the current schema.")
    }
  }

  return(list(
    schema = conversation_state$current_schema,
    description = conversation_state$description,
    iterations = conversation_state$iteration - 1,
    history = conversation_state$history,
    satisfied = conversation_state$satisfied
  ))
}

#' Create initial schema generation prompt
#' @param state Current conversation state
#' @return tidyprompt object
#' @noRd
create_schema_prompt <- function(state) {
  system_prompt <- paste(
    "You are a JSON schema expert. Your task is to generate schemas in a specific R list format",
    "with 'name', 'description', and 'schema' components.",
    "Focus on being precise and comprehensive."
  )

  user_prompt <- paste(
    "Generate a JSON schema structure for the following description:",
    paste0("Description: ", state$description),
    "",
    "Return a JSON object with exactly this structure:",
    "{",
    '  "name": "descriptive_name_for_schema",',
    '  "description": "Clear description of what this schema represents",',
    '  "schema": {',
    '    "type": "object",',
    '    "properties": {',
    '      "PropertyName1": {',
    '        "type": "string",',
    '        "description": "Description of this property"',
    "      },",
    '      "PropertyName2": {',
    '        "type": "string",',
    '        "description": "Description of this property"',
    "      }",
    "    },",
    '    "required": ["PropertyName1", "PropertyName2"],',
    '    "additionalProperties": false',
    "  }",
    "}",
    "",
    "Requirements:",
    "1. 'name' should be a concise identifier (snake_case)",
    "2. 'description' should explain the schema's purpose",
    "3. 'schema' should contain proper JSON schema with type, properties, required fields",
    "4. Use meaningful property names relevant to the description",
    "5. Set 'additionalProperties' to false for strict validation",
    "",
    "Return only the JSON object."
  )

  set_prompt(system = system_prompt, user = user_prompt)
}

#' Create refinement prompt based on user feedback
#' @param state Current conversation state
#' @return tidyprompt object
#' @noRd
create_refinement_prompt <- function(state) {
  system_prompt <- paste(
    "You are a JSON schema expert helping to refine an existing schema structure.",
    "Maintain the required format with 'name', 'description', and 'schema' components.",
    "Carefully consider the user's feedback and modify accordingly."
  )

  current_schema_text <- jsonlite::toJSON(state$current_schema, pretty = TRUE, auto_unbox = TRUE)

  user_prompt <- paste(
    "Here is the current schema structure:",
    "```json",
    current_schema_text,
    "```",
    "",
    paste0("User feedback: ", state$user_feedback),
    "",
    "Please modify the schema based on this feedback while maintaining the exact format:",
    "- Keep the 'name', 'description', and 'schema' structure",
    "- Update relevant parts based on feedback",
    "- Ensure 'schema' contains proper JSON schema format",
    "",
    "Return the updated complete JSON object."
  )

  set_prompt(system = system_prompt, user = user_prompt)
}

#' Generate schema step with error handling
#' @param prompt tidyprompt object
#' @param llm_client LLM provider
#' @param state Current state
#' @param verbose Show logs
#' @return List with schema and prompt used
#' @noRd
generate_schema_step <- function(prompt, llm_client, state, verbose) {
  # Create custom validation function for the specific schema format
  schema_validator <- function(response) {
    # Check if it's a valid list (R object from JSON)
    if (!is.list(response)) {
      return(tidyprompt::llm_feedback(
        "The response must be a valid JSON object. Please provide a properly formatted schema structure."
      ))
    }

    # Check for required top-level components
    required_components <- c("name", "description", "schema")
    missing_components <- setdiff(required_components, names(response))
    if (length(missing_components) > 0) {
      return(tidyprompt::llm_feedback(
        paste(
          "Missing required components:", paste(missing_components, collapse = ", "),
          ". Please include 'name', 'description', and 'schema' in your response."
        )
      ))
    }

    # Check schema structure
    schema_obj <- response$schema
    if (!is.list(schema_obj)) {
      return(tidyprompt::llm_feedback(
        "The 'schema' component must be a valid object. Please provide a proper schema structure."
      ))
    }

    # Check for basic schema properties
    if (!"type" %in% names(schema_obj)) {
      return(tidyprompt::llm_feedback(
        "The schema must include a 'type' property. Please add the appropriate type."
      ))
    }

    # For object types, check for properties
    if (schema_obj$type == "object" && (!"properties" %in% names(schema_obj) || length(schema_obj$properties) == 0)) {
      return(tidyprompt::llm_feedback(
        "Object type schemas should include 'properties'. Please add relevant properties to the schema."
      ))
    }

    # Check that name is a string
    if (!is.character(response$name) || length(response$name) != 1 || nchar(response$name) == 0) {
      return(tidyprompt::llm_feedback(
        "The 'name' field must be a non-empty string. Please provide a descriptive name."
      ))
    }

    # Check that description is a string
    if (!is.character(response$description) || length(response$description) != 1 || nchar(response$description) == 0) {
      return(tidyprompt::llm_feedback(
        "The 'description' field must be a non-empty string. Please provide a clear description."
      ))
    }

    return(response)
  }

  # Generate schema with validation
  result <- tryCatch(
    {
      schema <- prompt |>
        tidyprompt::answer_as_json(type = "auto") |>
        tidyprompt::prompt_wrap(validation_fn = schema_validator) |>
        tidyprompt::send_prompt(llm_client, verbose = verbose)

      list(schema = schema, prompt_used = prompt$construct_prompt_text())
    },
    error = function(e) {
      cli::cli_alert_danger(paste("Error generating schema:", e$message))
      stop("Failed to generate schema. Please check your LLM provider and try again.")
    }
  )

  return(result)
}

#' Get user feedback interactively
#' @param state Current conversation state
#' @param verbose Show logs
#' @return Updated conversation state
get_user_feedback <- function(state, verbose) {
  if (verbose) {
    cli::cli_rule()
  }

  cat("\nOptions:\n")
  cat("1. I'm satisfied with this schema\n")
  cat("2. I want to modify this schema\n")
  cat("3. Show me the full schema in detail\n")
  cat("4. Show me the R code format\n")
  cat("5. Start over with a new description\n")

  choice <- readline(prompt = "Enter your choice (1-5): ")

  switch(choice,
    "1" = {
      state$satisfied <- TRUE
      if (verbose) cli::cli_alert_success("Great! Schema finalized.")
    },
    "2" = {
      feedback <- readline(prompt = "What changes would you like? Describe specifically: ")
      if (nchar(trimws(feedback)) > 0) {
        state$user_feedback <- feedback
        if (verbose) cli::cli_alert_info("Feedback recorded. Generating refined schema...")
      } else {
        if (verbose) cli::cli_alert_warning("No feedback provided. Keeping current schema.")
        state$satisfied <- TRUE
      }
    },
    "3" = {
      if (verbose) {
        cli::cli_h3("Full Schema Details")
        cat(jsonlite::toJSON(state$current_schema, pretty = TRUE, auto_unbox = TRUE))
        cat("\n")
      }
      # Recursively ask for feedback
      state <- get_user_feedback(state, verbose)
    },
    "4" = {
      if (verbose) {
        cli::cli_h3("R Code Format")
        print_r_code(state$current_schema)
      }
      # Recursively ask for feedback
      state <- get_user_feedback(state, verbose)
    },
    "5" = {
      new_desc <- readline(prompt = "Enter new description: ")
      if (nchar(trimws(new_desc)) > 0) {
        state$description <- new_desc
        state$user_feedback <- "Start over with new description"
        state$iteration <- 1 # Reset iteration
        if (verbose) cli::cli_alert_info("Starting over with new description...")
      } else {
        if (verbose) cli::cli_alert_warning("No new description provided.")
        state <- get_user_feedback(state, verbose)
      }
    },
    {
      if (verbose) cli::cli_alert_warning("Invalid choice. Please try again.")
      state <- get_user_feedback(state, verbose)
    }
  )

  return(state)
}

#' Print schema preview
#' @param schema_obj Complete schema object with name, description, schema
#' @noRd
print_schema_preview <- function(schema_obj) {
  cli::cli_h3("Generated Schema Preview")

  # Show top-level info
  cat(paste("Name:", schema_obj$name %||% "not specified", "\n"))
  cat(paste("Description:", schema_obj$description %||% "not specified", "\n"))
  cat("\n")

  # Show schema details
  schema <- schema_obj$schema
  if (!is.null(schema)) {
    cat(paste("Schema Type:", schema$type %||% "not specified", "\n"))

    if ("properties" %in% names(schema) && length(schema$properties) > 0) {
      cat("Properties:\n")
      for (prop_name in names(schema$properties)) {
        prop <- schema$properties[[prop_name]]
        prop_type <- prop$type %||% "unspecified"
        prop_desc <- if ("description" %in% names(prop)) paste(" -", prop$description) else ""
        cat(paste("  -", prop_name, "(", prop_type, ")", prop_desc, "\n"))
      }
    }

    if ("required" %in% names(schema) && length(schema$required) > 0) {
      cat(paste("Required fields:", paste(schema$required, collapse = ", "), "\n"))
    }

    if ("additionalProperties" %in% names(schema)) {
      cat(paste("Additional properties allowed:", schema$additionalProperties, "\n"))
    }
  }

  cat("\n")
}

#' Extract only the schema part from generated result
#' @param schema_result Result from generate_json_schema
#' @return Just the schema portion for use with tidyprompt
#' @export
extract_schema_only <- function(schema_result) {
  if ("schema" %in% names(schema_result)) {
    return(schema_result$schema)
  } else {
    return(schema_result) # In case it's already just the schema
  }
}

#' Format schema result as R list assignment code
#' @param schema_result Complete schema object
#' @param variable_name Variable name for the assignment (default: "json_schema")
#' @return Character string with R code
#' @noRd
format_as_r_code <- function(schema_result, variable_name = "json_schema") {
  # Convert to properly formatted R code
  r_code <- paste0(variable_name, " <- ")
  r_code <- paste0(r_code, deparse(schema_result, width.cutoff = 80))
  r_code <- paste(r_code, collapse = "\n")
  return(r_code)
}

#' Print formatted R code to console
#' @param schema_result Complete schema object
#' @param variable_name Variable name for the assignment
#' @noRd
print_r_code <- function(schema_result, variable_name = "json_schema") {
  cat("# Generated R code:\n")
  cat(format_as_r_code(schema_result, variable_name))
  cat("\n")
}
