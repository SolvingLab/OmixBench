#' Extract Examples from a Package Function
#'
#' This function extracts and cleans the examples section from a specific
#' function's documentation in an R package. It uses the `tools` package
#' to access the Rd database and extracts examples using `tools::Rd2ex()`.
#' The output is cleaned to remove metadata headers and formatting artifacts.
#'
#' @param package_name A character string specifying the name of the package
#' @param function_name A character string specifying the name of the function
#'
#' @return A character string containing the cleaned examples code, or `NA`
#'   if no examples are found or an error occurs
#'
#' @examples
#' \dontrun{
#' # Extract examples from ggplot2's geom_point function
#' examples <- extract_function_examples("ggplot2", "geom_point")
#' cat(examples)
#' }
#'
#' @export
extract_function_examples <- function(package_name, function_name) {
  tryCatch(
    {
      # Locate the Rd database file for the package
      rdbfile <- file.path(find.package(package_name), "help", package_name)

      # Fetch the documentation for the specific function
      rdb <- tools:::fetchRdDB(rdbfile, key = function_name)
      if (is.null(rdb)) {
        return(NA)
      }

      # Extract examples using tools::Rd2ex
      examples_text <- capture.output(tools::Rd2ex(rdb))
      if (length(examples_text) == 0) {
        return(NA)
      }

      # Clean examples: remove header information, keep only actual code
      # Look for the "### ** Examples" line
      examples_start <- which(grepl("### \\*\\* Examples", examples_text))
      if (length(examples_start) > 0) {
        # Start from the line after "### ** Examples"
        clean_examples <- examples_text[(examples_start[1] + 1):length(examples_text)]

        # Remove leading empty lines
        while (length(clean_examples) > 0 && clean_examples[1] == "") {
          clean_examples <- clean_examples[-1]
        }
        return(paste(clean_examples, collapse = "\n"))
      } else {
        # If standard format not found, return all (remove obvious headers)
        # Remove lines starting with ###
        clean_examples <- examples_text[!grepl("^###", examples_text)]
        return(paste(clean_examples, collapse = "\n"))
      }
    },
    error = function(e) {
      return(NA)
    }
  )
}


#' Extract Documentation for a Single Function
#'
#' This function extracts comprehensive documentation for a specific function
#' in an R package. It retrieves the function's title, description, usage
#' pattern, parameters with defaults, return value, and examples. Unlike
#' \code{\link{extract_package_docs}}, this function targets a single function
#' and returns its documentation directly rather than in a list.
#'
#' @param package_name A character string specifying the name of the package
#' @param function_name A character string specifying the name of the function
#'
#' @return A list containing the function's documentation with the following elements:
#'   \itemize{
#'     \item \code{package}: Package name
#'     \item \code{function_name}: Function name
#'     \item \code{title}: Function title
#'     \item \code{description}: Detailed description
#'     \item \code{usage}: Usage syntax
#'     \item \code{parameters}: Data frame of parameters with defaults
#'     \item \code{return_value}: Description of return value
#'     \item \code{examples}: Example code
#'     \item \code{formatted_arguments}: Formatted parameter descriptions
#'     \item \code{simple_arguments}: Required parameters only
#'   }
#'   Returns \code{NULL} if the function is not found or an error occurs.
#'
#' @examples
#' \dontrun{
#' # Extract documentation for dplyr's filter function
#' filter_doc <- extract_function_docs("dplyr", "filter")
#'
#' # View the description
#' cat(filter_doc$description)
#'
#' # View the parameters
#' print(filter_doc$parameters)
#' }
#'
#' @export
extract_function_docs <- function(package_name, function_name) {
  tryCatch(
    {
      # Get help database for the package
      help_db <- tools::Rd_db(package_name)

      # Search through documentation files for the function
      for (rd_file in names(help_db)) {
        rd_content <- help_db[[rd_file]]
        tags <- tools:::RdTags(rd_content)

        # Get function aliases
        alias_indices <- which(tags == "\\alias")
        if (length(alias_indices) == 0) next

        # Check if this file contains our target function
        for (i in alias_indices) {
          func_name <- clean_text(extract_text(rd_content[[i]]))
          if (is.na(func_name) || func_name != function_name) next

          # Verify it's an actual function
          if (!exists(func_name, envir = asNamespace(package_name))) next
          obj <- get(func_name, envir = asNamespace(package_name))
          if (!is.function(obj)) next

          # Extract documentation sections
          title <- get_section(rd_content, tags, "\\title")
          description <- get_section(rd_content, tags, "\\description")
          usage <- get_section(rd_content, tags, "\\usage")
          arguments_text <- get_section(rd_content, tags, "\\arguments")
          return_value <- get_section(rd_content, tags, "\\value")

          # Extract examples using the proper method
          examples <- extract_function_examples(package_name, func_name)

          # Process parameters
          params <- get_parameters(usage, func_name)

          if (nrow(params) > 0) {
            param_descs <- get_param_descriptions(arguments_text, params$name)
            params$description <- param_descs

            # Format parameter descriptions (newline separated)
            all_params <- character()
            required_params <- character()

            for (j in 1:nrow(params)) {
              desc <- if (is.na(params$description[j])) "[No description]" else params$description[j]
              param_line <- paste0(params$name[j], ": ", desc)
              all_params <- c(all_params, param_line)

              if (!params$has_default[j]) {
                required_params <- c(required_params, param_line)
              }
            }

            formatted_arguments <- paste(all_params, collapse = "\n")
            simple_arguments <- paste(required_params, collapse = "\n")
          } else {
            formatted_arguments <- ""
            simple_arguments <- ""
          }

          # Return the result
          result <- list(
            package = package_name,
            function_name = func_name,
            title = title,
            description = description,
            usage = usage,
            parameters = params,
            return_value = return_value,
            examples = examples,
            formatted_arguments = formatted_arguments,
            simple_arguments = simple_arguments
          )

          return(result)
        }
      }

      # Function not found
      warning("Function '", function_name, "' not found in package '", package_name, "'")
      return(NULL)
    },
    error = function(e) {
      warning("Error extracting documentation: ", e$message)
      return(NULL)
    }
  )
}


#' Extract Documentation for All Functions in a Package
#'
#' This function extracts comprehensive documentation for all exported functions
#' in an R package. It retrieves function names, titles, descriptions, usage
#' patterns, parameters with defaults, return values, and examples. The function
#' processes the package's Rd database and returns a structured list containing
#' all documentation elements.
#'
#' @param package_name A character string specifying the name of the package
#'   to extract documentation from
#'
#' @return A named list where each element corresponds to a function in the
#'   package. Each function's documentation includes:
#'   \itemize{
#'     \item \code{package}: Package name
#'     \item \code{function_name}: Function name
#'     \item \code{title}: Function title
#'     \item \code{description}: Detailed description
#'     \item \code{usage}: Usage syntax
#'     \item \code{parameters}: Data frame of parameters with defaults
#'     \item \code{return_value}: Description of return value
#'     \item \code{examples}: Example code
#'     \item \code{formatted_arguments}: Formatted parameter descriptions
#'     \item \code{simple_arguments}: Required parameters only
#'   }
#'
#' @examples
#' \dontrun{
#' # Extract documentation from dplyr
#' docs <- extract_package_docs("dplyr")
#'
#' # View documentation for a specific function
#' docs$filter
#'
#' # List all extracted functions
#' names(docs)
#' }
#'
#' @export
extract_package_docs <- function(package_name) {
  cat("Extracting package:", package_name, "\n")

  # Get help database for the package
  help_db <- tools::Rd_db(package_name)
  cat("Documentation files:", length(help_db), "\n")

  results <- list()
  count <- 0

  # Process each documentation file
  for (rd_file in names(help_db)) {
    rd_content <- help_db[[rd_file]]
    tags <- tools:::RdTags(rd_content)

    # Get function aliases
    alias_indices <- which(tags == "\\alias")
    if (length(alias_indices) == 0) next

    for (i in alias_indices) {
      func_name <- clean_text(extract_text(rd_content[[i]]))
      if (is.na(func_name) || !grepl("^[a-zA-Z]", func_name)) next

      # Check if it's an actual function
      if (!exists(func_name, envir = asNamespace(package_name))) next
      obj <- get(func_name, envir = asNamespace(package_name))
      if (!is.function(obj)) next

      # Extract documentation sections
      title <- get_section(rd_content, tags, "\\title")
      description <- get_section(rd_content, tags, "\\description")
      usage <- get_section(rd_content, tags, "\\usage")
      arguments_text <- get_section(rd_content, tags, "\\arguments")
      return_value <- get_section(rd_content, tags, "\\value")

      # Extract examples using the proper method
      examples <- extract_function_examples(package_name, func_name)

      # Process parameters
      params <- get_parameters(usage, func_name)

      if (nrow(params) > 0) {
        param_descs <- get_param_descriptions(arguments_text, params$name)
        params$description <- param_descs

        # Format parameter descriptions (newline separated)
        all_params <- character()
        required_params <- character()

        for (j in 1:nrow(params)) {
          desc <- if (is.na(params$description[j])) "[No description]" else params$description[j]
          param_line <- paste0(params$name[j], ": ", desc)
          all_params <- c(all_params, param_line)

          if (!params$has_default[j]) {
            required_params <- c(required_params, param_line)
          }
        }

        formatted_arguments <- paste(all_params, collapse = "\n")
        simple_arguments <- paste(required_params, collapse = "\n")
      } else {
        formatted_arguments <- ""
        simple_arguments <- ""
      }

      # Save results
      results[[func_name]] <- list(
        package = package_name,
        function_name = func_name,
        title = title,
        description = description,
        usage = usage,
        parameters = params,
        return_value = return_value,
        examples = examples,
        formatted_arguments = formatted_arguments,
        simple_arguments = simple_arguments
      )

      count <- count + 1
      cat("[OK]", func_name, "\n")
    }
  }

  cat("Complete! Extracted", count, "functions\n")
  return(results)
}


# Helper Functions --------------------------------------------------------

#' @noRd
clean_text <- function(text) {
  if (is.null(text) || length(text) == 0 || text == "") {
    return(NA)
  }

  # Remove LaTeX-style commands
  text <- gsub("\\\\[a-zA-Z]+\\{[^}]*\\}", "", text)
  # Remove curly braces
  text <- gsub("\\{|\\}", "", text)
  # Collapse multiple spaces
  text <- gsub("\\s+", " ", text)
  # Trim whitespace
  text <- trimws(text)

  if (nchar(text) == 0) {
    return(NA)
  }
  return(text)
}


#' @noRd
extract_text <- function(rd_element) {
  if (is.character(rd_element)) {
    return(paste(rd_element, collapse = " "))
  } else if (is.list(rd_element)) {
    texts <- sapply(rd_element, extract_text)
    return(paste(texts[texts != ""], collapse = " "))
  }
  return("")
}


#' @noRd
get_section <- function(rd_content, tags, section_name) {
  indices <- which(tags == section_name)
  if (length(indices) == 0) {
    return(NA)
  }

  content <- ""
  for (idx in indices) {
    content <- paste(content, extract_text(rd_content[[idx]]))
  }
  return(clean_text(content))
}


#' @noRd
get_parameters <- function(usage_text, func_name) {
  if (is.na(usage_text)) {
    return(data.frame())
  }

  # Find the target function's signature (handles multiple functions)
  func_pattern <- paste0("\\b", func_name, "\\s*\\(")
  func_match <- regexpr(func_pattern, usage_text, perl = TRUE)
  if (func_match == -1) {
    return(data.frame())
  }

  # Extract from function name onwards
  start_pos <- func_match
  usage_from_func <- substr(usage_text, start_pos, nchar(usage_text))

  # Find the complete parentheses content of the function call (handle nested parentheses)
  paren_start <- regexpr("\\(", usage_from_func)
  if (paren_start == -1) {
    return(data.frame())
  }

  # From the opening parenthesis, find the matching closing parenthesis
  paren_count <- 0
  end_pos <- paren_start
  for (i in paren_start:nchar(usage_from_func)) {
    char <- substr(usage_from_func, i, i)
    if (char == "(") paren_count <- paren_count + 1
    if (char == ")") paren_count <- paren_count - 1
    if (paren_count == 0) {
      end_pos <- i
      break
    }
  }

  # Extract parameter text within parentheses
  if (paren_count != 0) {
    return(data.frame())
  } # Unmatched parentheses
  params_text <- substr(usage_from_func, paren_start + 1, end_pos - 1)
  params_text <- trimws(params_text)

  if (nchar(params_text) == 0) {
    return(data.frame())
  }

  # Smart parameter splitting - considering nested parentheses
  params <- split_parameters_smart(params_text)

  result <- data.frame(
    name = character(),
    has_default = logical(),
    default_value = character(),
    stringsAsFactors = FALSE
  )

  for (param in params) {
    param <- trimws(param)
    if (nchar(param) > 0) {
      if (grepl("=", param)) {
        # Find the first equals sign not within parentheses
        eq_pos <- find_main_equals(param)
        if (eq_pos > 0) {
          param_name <- trimws(substr(param, 1, eq_pos - 1))
          default_val <- trimws(substr(param, eq_pos + 1, nchar(param)))

          result <- rbind(result, data.frame(
            name = param_name,
            has_default = TRUE,
            default_value = default_val,
            stringsAsFactors = FALSE
          ))
        } else {
          # If no main equals found, treat as no default
          result <- rbind(result, data.frame(
            name = trimws(param),
            has_default = FALSE,
            default_value = NA,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        result <- rbind(result, data.frame(
          name = trimws(param),
          has_default = FALSE,
          default_value = NA,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(result)
}


#' @noRd
split_parameters_smart <- function(params_text) {
  if (nchar(params_text) == 0) {
    return(character(0))
  }

  params <- character()
  current_param <- ""
  paren_count <- 0

  for (i in 1:nchar(params_text)) {
    char <- substr(params_text, i, i)

    if (char == "(") {
      paren_count <- paren_count + 1
      current_param <- paste0(current_param, char)
    } else if (char == ")") {
      paren_count <- paren_count - 1
      current_param <- paste0(current_param, char)
    } else if (char == "," && paren_count == 0) {
      # Only split on commas that are not within parentheses
      if (nchar(trimws(current_param)) > 0) {
        params <- c(params, current_param)
      }
      current_param <- ""
    } else {
      current_param <- paste0(current_param, char)
    }
  }

  # Add the last parameter
  if (nchar(trimws(current_param)) > 0) {
    params <- c(params, current_param)
  }

  return(params)
}


#' @noRd
find_main_equals <- function(param_text) {
  paren_count <- 0

  for (i in 1:nchar(param_text)) {
    char <- substr(param_text, i, i)

    if (char == "(") {
      paren_count <- paren_count + 1
    } else if (char == ")") {
      paren_count <- paren_count - 1
    } else if (char == "=" && paren_count == 0) {
      return(i)
    }
  }

  return(-1) # Not found
}


#' @noRd
get_param_descriptions <- function(arguments_text, param_names) {
  descriptions <- rep(NA, length(param_names))

  if (is.na(arguments_text) || length(param_names) == 0) {
    return(descriptions)
  }

  # Find positions of each parameter name in the arguments text
  positions <- numeric(length(param_names))
  for (i in seq_along(param_names)) {
    param <- param_names[i]
    pos <- regexpr(param, arguments_text, fixed = TRUE)
    positions[i] <- if (pos == -1) 0 else pos
  }

  # Filter valid positions and sort
  valid_indices <- which(positions > 0)
  if (length(valid_indices) == 0) {
    return(descriptions)
  }

  sorted_order <- order(positions[valid_indices])
  sorted_indices <- valid_indices[sorted_order]

  # Extract description for each parameter
  for (i in seq_along(sorted_indices)) {
    idx <- sorted_indices[i]
    param <- param_names[idx]
    start_pos <- positions[idx] + nchar(param)

    # End position is the start of the next parameter, or end of text
    if (i < length(sorted_indices)) {
      next_idx <- sorted_indices[i + 1]
      end_pos <- positions[next_idx] - 1
    } else {
      end_pos <- nchar(arguments_text)
    }

    # Extract and clean description
    desc <- substr(arguments_text, start_pos, end_pos)
    desc <- trimws(desc)
    desc <- gsub("^[^a-zA-Z]*", "", desc)

    if (nchar(desc) > 3) {
      descriptions[idx] <- desc
    }
  }

  return(descriptions)
}


#' Extract Documentation from Package Tarball
#'
#' This function extracts comprehensive documentation from an uninstalled R package
#' in tar.gz format. It processes the .Rd files directly from the tarball without
#' requiring package installation. The return structure is identical to
#' \code{\link{extract_package_docs}}, allowing you to use \code{\link{generate_qa_from_docs}}
#' directly on the results.
#'
#' This is particularly useful for:
#' \itemize{
#'   \item Batch processing large collections of packages (e.g., Bioconductor)
#'   \item Extracting documentation from packages that cannot be installed
#'   \item Processing packages in environments without installation privileges
#' }
#'
#' @param tar_path A character string specifying the full path to the package
#'   tar.gz file
#'
#' @return A named list where each element corresponds to a function in the
#'   package. Each function's documentation includes:
#'   \itemize{
#'     \item \code{package}: Package name
#'     \item \code{function_name}: Function name
#'     \item \code{title}: Function title
#'     \item \code{description}: Detailed description
#'     \item \code{usage}: Usage syntax
#'     \item \code{parameters}: Data frame of parameters with defaults
#'     \item \code{return_value}: Description of return value
#'     \item \code{examples}: Example code
#'     \item \code{formatted_arguments}: Formatted parameter descriptions
#'     \item \code{simple_arguments}: Required parameters only
#'   }
#'   Returns \code{NULL} if extraction fails or no documentation is found.
#'
#' @examples
#' \dontrun{
#' # Extract documentation from a tarball
#' docs <- extract_docs_from_tarball("path/to/dplyr_1.0.0.tar.gz")
#'
#' # Use with existing QA generation function
#' qa_pairs <- generate_qa_from_docs(docs)
#'
#' # View documentation for a specific function
#' docs$filter
#'
#' # List all extracted functions
#' names(docs)
#' }
#'
#' @export
extract_docs_from_tarball <- function(tar_path) {
  if (!file.exists(tar_path)) {
    warning("Tarball file not found: ", tar_path)
    return(NULL)
  }

  package_name <- sub("_.*\\.tar\\.gz$", "", basename(tar_path))
  cat("Extracting package:", package_name, "\n")

  # Create package-specific temporary directory
  temp_dir <- tempdir()
  pkg_temp_dir <- file.path(temp_dir, "tarball_extraction", package_name)
  if (dir.exists(pkg_temp_dir)) unlink(pkg_temp_dir, recursive = TRUE)
  dir.create(pkg_temp_dir, recursive = TRUE)

  # Extract tar.gz file
  tryCatch(
    {
      untar(tar_path, exdir = pkg_temp_dir)
    },
    error = function(e) {
      cat("Extraction failed:", package_name, "-", e$message, "\n")
      unlink(pkg_temp_dir, recursive = TRUE)
      return(NULL)
    }
  )

  # Find extracted package directory
  extracted_dirs <- list.dirs(pkg_temp_dir, recursive = FALSE)
  if (length(extracted_dirs) == 0) {
    cat("No extracted directory found:", package_name, "\n")
    unlink(pkg_temp_dir, recursive = TRUE)
    return(NULL)
  }

  pkg_dir <- extracted_dirs[1]
  man_dir <- file.path(pkg_dir, "man")

  # Check if man directory exists
  if (!dir.exists(man_dir)) {
    cat("No man directory found:", package_name, "\n")
    unlink(pkg_temp_dir, recursive = TRUE)
    return(NULL)
  }

  # Get all .Rd files
  rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  if (length(rd_files) == 0) {
    cat("No .Rd files found:", package_name, "\n")
    unlink(pkg_temp_dir, recursive = TRUE)
    return(NULL)
  }

  cat("Documentation files:", length(rd_files), "\n")

  # Initialize results list
  results <- list()
  count <- 0

  # Process each .Rd file
  for (rd_file in rd_files) {
    tryCatch(
      {
        # Read and parse .Rd file
        rd_content <- tools::parse_Rd(rd_file)
        tags <- tools:::RdTags(rd_content)

        # Get function aliases
        alias_indices <- which(tags == "\\alias")
        if (length(alias_indices) == 0) next

        for (i in alias_indices) {
          func_name <- clean_text(extract_text(rd_content[[i]]))
          if (is.na(func_name) || !grepl("^[a-zA-Z]", func_name)) next

          # Extract documentation sections
          title <- get_section(rd_content, tags, "\\title")
          description <- get_section(rd_content, tags, "\\description")
          usage <- get_section(rd_content, tags, "\\usage")
          arguments_text <- get_section(rd_content, tags, "\\arguments")
          return_value <- get_section(rd_content, tags, "\\value")

          # Extract examples with proper formatting (preserving line breaks)
          examples <- extract_examples_from_tarball(rd_content, tags)

          # Process parameters
          params <- get_parameters(usage, func_name)

          if (nrow(params) > 0) {
            param_descs <- get_param_descriptions(arguments_text, params$name)
            params$description <- param_descs

            # Format parameter descriptions (newline separated)
            all_params <- character()
            required_params <- character()

            for (j in 1:nrow(params)) {
              desc <- if (is.na(params$description[j])) "[No description]" else params$description[j]
              param_line <- paste0(params$name[j], ": ", desc)
              all_params <- c(all_params, param_line)

              if (!params$has_default[j]) {
                required_params <- c(required_params, param_line)
              }
            }

            formatted_arguments <- paste(all_params, collapse = "\n")
            simple_arguments <- paste(required_params, collapse = "\n")
          } else {
            formatted_arguments <- ""
            simple_arguments <- ""
          }

          # Save results in the same structure as extract_package_docs
          results[[func_name]] <- list(
            package = package_name,
            function_name = func_name,
            title = title,
            description = description,
            usage = usage,
            parameters = params,
            return_value = return_value,
            examples = examples,
            formatted_arguments = formatted_arguments,
            simple_arguments = simple_arguments
          )

          count <- count + 1
          cat("[OK]", func_name, "\n")
        }
      },
      error = function(e) {
        cat("Failed to process .Rd file:", basename(rd_file), "-", e$message, "\n")
      }
    )
  }

  # Clean up temporary directory
  unlink(pkg_temp_dir, recursive = TRUE)

  cat("Complete! Extracted", count, "functions\n")
  return(results)
}


# Helper Functions --------------------------------------------------------

#' @noRd
extract_text_preserve_lines <- function(rd_element) {
  if (is.character(rd_element)) {
    return(paste(rd_element, collapse = "\n"))
  } else if (is.list(rd_element)) {
    text_parts <- character()
    for (elem in rd_element) {
      sub_text <- extract_text_preserve_lines(elem)
      if (nchar(sub_text) > 0) {
        text_parts <- c(text_parts, sub_text)
      }
    }
    return(paste(text_parts, collapse = "\n"))
  }
  return("")
}


#' @noRd
extract_examples_from_tarball <- function(rd_content, tags) {
  example_indices <- which(tags == "\\examples")
  if (length(example_indices) == 0) {
    return(NA)
  }

  # Extract examples content directly without clean_text
  examples_lines <- character()
  for (idx in example_indices) {
    # Recursively extract text, preserving structure
    example_text <- extract_text_preserve_lines(rd_content[[idx]])
    if (nchar(example_text) > 0) {
      examples_lines <- c(examples_lines, example_text)
    }
  }

  if (length(examples_lines) == 0) {
    return(NA)
  }

  # Merge all examples, preserving line breaks
  full_examples <- paste(examples_lines, collapse = "\n")

  # Minimal cleaning, preserve line breaks and code structure
  # Preserve content within LaTeX commands
  full_examples <- gsub("\\\\[a-zA-Z]+\\{([^}]*)\\}", "\\1", full_examples)
  # Remove other LaTeX commands
  full_examples <- gsub("\\\\[a-zA-Z]+", "", full_examples)
  full_examples <- trimws(full_examples)

  if (nchar(full_examples) < 10) {
    return(NA)
  }

  return(full_examples)
}
