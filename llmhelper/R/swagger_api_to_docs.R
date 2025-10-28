# Swagger API to LangChain Tools Converter
# Generate LangChain tools from Swagger/OpenAPI documentation

# Using explicit package references for clarity
# Required packages: httr, jsonlite, stringr
# Install with: install.packages(c("httr", "jsonlite", "stringr"))

#' Convert Swagger API documentation to LangChain tools
#'
#' @param swagger_base_url Character. The base URL of the Swagger/OpenAPI server
#' @param output_dir Character. Directory where tool files will be saved
#' @param verbose Logical. Whether to print progress information
#'
#' @return List with two components:
#'   - tool_details: Complete tool definitions for LangChain usage
#'   - tool_summary: Data frame with overview of all tools (name, method, params, etc.)
#'
#' @export
swagger_api_to_docs <- function(
    swagger_base_url = "http://solvinglab.top:5002",
    output_dir = "langchain_tools",
    verbose = TRUE) {
  if (verbose) cat("🌐 Fetching API specification...\n")

  # Get OpenAPI specification
  api_spec <- get_openapi_spec(swagger_base_url)
  if (is.null(api_spec)) {
    stop("Failed to retrieve API specification")
  }

  if (verbose) cat("📋 Parsing API endpoints...\n")

  # Parse to tools
  tools <- parse_openapi_to_tools(api_spec, swagger_base_url)

  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Save tools
  if (verbose) cat("💾 Saving tool files...\n")
  for (tool in tools) {
    filename <- file.path(output_dir, paste0(tool$name, ".json"))
    jsonlite::write_json(tool, filename, pretty = TRUE, auto_unbox = TRUE)
  }

  # Create summary table
  tool_summary <- create_summary_table(tools)

  # Save summary files
  jsonlite::write_json(tools, file.path(output_dir, "all_tools.json"), pretty = TRUE, auto_unbox = TRUE)
  write.csv(tool_summary, file.path(output_dir, "tools_summary.csv"), row.names = FALSE)

  if (verbose) {
    cat("✅ Generated", length(tools), "tools\n")
    cat("📊 Summary columns: tool_name, display_name, method, endpoint, description, total_params, required_params, tags\n")
    cat("📄 Files saved: tools_summary.csv, all_tools.json\n")
  }

  # Return structured list
  return(list(
    tool_details = tools,
    tool_summary = tool_summary
  ))
}

#' Get OpenAPI specification from common endpoints
#' @noRd
get_openapi_spec <- function(base_url) {
  endpoints <- c("/openapi.json", "/swagger.json", "/__docs__/openapi.json")

  for (endpoint in endpoints) {
    url <- paste0(base_url, endpoint)

    tryCatch(
      {
        response <- httr::GET(url)
        if (httr::status_code(response) == 200) {
          spec <- httr::content(response, as = "parsed")
          return(spec)
        }
      },
      error = function(e) NULL
    )
  }

  return(NULL)
}

#' Parse OpenAPI specification to LangChain tools
#' @noRd
parse_openapi_to_tools <- function(spec, base_url) {
  tools <- list()
  api_title <- spec$info$title %||% "API"

  # Process each path and method
  for (path in names(spec$paths)) {
    path_info <- spec$paths[[path]]

    for (method in names(path_info)) {
      if (method %in% c("get", "post", "put", "delete", "patch")) {
        operation <- path_info[[method]]

        tool <- create_tool_from_operation(path, method, operation, base_url)

        if (!is.null(tool)) {
          tools[[length(tools) + 1]] <- tool
        }
      }
    }
  }

  return(tools)
}

#' Create tool from OpenAPI operation
#' @noRd
create_tool_from_operation <- function(path, method, operation, base_url) {
  # Generate tool name
  operation_id <- operation$operationId %||% paste0(method, "_", path)
  tool_name <- stringr::str_replace_all(operation_id, "[^a-zA-Z0-9_]", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace("^_+|_+$", "")

  if (stringr::str_detect(tool_name, "^[0-9]")) {
    tool_name <- paste0("api_", tool_name)
  }

  # Extract parameters
  parameters <- extract_parameters(operation)

  # Build tool
  tool <- list(
    name = tool_name,
    display_name = operation$summary %||% paste(toupper(method), path),
    description = operation$description %||% operation$summary %||% paste("API endpoint:", toupper(method), path),
    api_info = list(
      method = toupper(method),
      url = paste0(base_url, path),
      content_type = "application/json"
    ),
    parameters = parameters,
    tags = operation$tags %||% character(0),
    langchain_config = list(
      return_direct = FALSE,
      handle_tool_error = TRUE
    )
  )

  return(tool)
}

#' Extract parameters from operation
#' @noRd
extract_parameters <- function(operation) {
  properties <- list()
  required <- character()

  # URL parameters (path and query)
  if (!is.null(operation$parameters)) {
    for (param in operation$parameters) {
      properties[[param$name]] <- list(
        type = param$schema$type %||% "string",
        description = param$description %||% ""
      )

      if (param$required %||% FALSE) {
        required <- c(required, param$name)
      }
    }
  }

  # Request body parameters
  if (!is.null(operation$requestBody)) {
    json_content <- operation$requestBody$content$`application/json`

    if (!is.null(json_content$schema$properties)) {
      for (prop_name in names(json_content$schema$properties)) {
        prop <- json_content$schema$properties[[prop_name]]

        properties[[prop_name]] <- list(
          type = prop$type %||% "string",
          description = prop$description %||% ""
        )

        if (prop_name %in% (json_content$schema$required %||% character())) {
          required <- c(required, prop_name)
        }
      }
    }
  }

  return(list(
    type = "object",
    properties = properties,
    required = unique(required)
  ))
}

#' Create summary table of tools
#' @noRd
create_summary_table <- function(tools) {
  # Extract key information for each tool
  summary_data <- data.frame(
    tool_name = sapply(tools, function(x) x$name),
    display_name = sapply(tools, function(x) x$display_name),
    method = sapply(tools, function(x) x$api_info$method),
    endpoint = sapply(tools, function(x) {
      # Extract just the path part from the full URL
      url <- x$api_info$url
      path <- sub("^https?://[^/]+", "", url)
      return(path)
    }),
    description = sapply(tools, function(x) substr(x$description, 1, 80)),
    total_params = sapply(tools, function(x) length(x$parameters$properties)),
    required_params = sapply(tools, function(x) length(x$parameters$required)),
    tags = sapply(tools, function(x) paste(x$tags, collapse = ", ")),
    stringsAsFactors = FALSE
  )

  return(summary_data)
}

# NULL coalescing operator
# `%||%` <- function(x, y) if (is.null(x)) y else x

# Example usage:
# result <- swagger_api_to_docs("http://solvinglab.top:5002")
#
# # Access detailed tool definitions
# tools <- result$tool_details
#
# # View summary table
# summary <- result$tool_summary
# View(summary)
#
# # Check specific tool
# print(tools[[1]]$api_info)
