#' Generate Question-Answer Pairs from Package Documentation
#'
#' This function generates comprehensive question-answer (QA) pairs optimized
#' for R code generation and learning from package documentation. It processes
#' the documentation structure returned by \code{\link{extract_package_docs}}
#' and creates multiple types of QA pairs for each function, including package
#' context, purpose, syntax, parameters, examples, and return values.
#'
#' @param docs A named list of function documentation, typically returned by
#'   \code{\link{extract_package_docs}}. Each element should contain function
#'   metadata including package name, function name, description, usage,
#'   parameters, examples, and return value.
#'
#' @return A list of QA pairs where each element contains:
#'   \itemize{
#'     \item \code{question}: A natural language question about the function
#'     \item \code{answer}: The corresponding answer, often containing R code
#'     \item \code{package}: The source package name
#'     \item \code{function_name}: The function name
#'     \item \code{qa_type}: The type of QA pair (e.g., "package_context",
#'       "purpose", "syntax", "examples", etc.)
#'   }
#'   The function generates up to 8 different types of QA pairs per function:
#'   \enumerate{
#'     \item \strong{package_context}: Which package to load
#'     \item \strong{purpose}: When and why to use the function
#'     \item \strong{syntax}: Correct function call syntax
#'     \item \strong{essential_parameters}: Required arguments
#'     \item \strong{all_parameters}: All available parameters
#'     \item \strong{examples}: Working code examples
#'     \item \strong{code_completion}: Parameter completion patterns
#'     \item \strong{return_value}: Expected output
#'   }
#'
#' @examples
#' \dontrun{
#' # Extract documentation from a package
#' docs <- extract_package_docs("dplyr")
#'
#' # Generate QA pairs
#' qa_pairs <- generate_qa_from_docs(docs)
#'
#' # View a specific QA pair
#' qa_pairs[[1]]
#'
#' # Filter by QA type
#' examples_qa <- Filter(function(x) x$qa_type == "examples", qa_pairs)
#' }
#'
#' @export
generate_qa_from_docs <- function(docs) {
  all_qa_pairs <- list()
  skipped_functions <- character()

  for (func_name in names(docs)) {
    func_data <- docs[[func_name]]

    # Get description (title as fallback)
    description <- get_description(func_data)
    if (is.na(description)) {
      skipped_functions <- c(skipped_functions, func_name)
      next
    }

    # Generate QA pairs for this function
    func_qa_pairs <- generate_function_qa_pairs(func_data, description)

    # Add to main list
    all_qa_pairs <- c(all_qa_pairs, func_qa_pairs)
  }

  cat("Generated", length(all_qa_pairs), "QA pairs from", length(docs), "functions\n")
  cat("Skipped", length(skipped_functions), "functions due to missing description/title\n")

  return(all_qa_pairs)
}


# Helper Functions --------------------------------------------------------

#' @noRd
get_description <- function(func_data) {
  if (!is.na(func_data$description) && nchar(func_data$description) > 0) {
    return(func_data$description)
  } else if (!is.na(func_data$title) && nchar(func_data$title) > 0) {
    return(func_data$title)
  } else {
    return(NA)
  }
}


#' @noRd
generate_function_qa_pairs <- function(func_data, description) {
  qa_pairs <- list()
  pkg <- func_data$package
  fn <- func_data$function_name

  # 1. Package context - critical for imports
  qa_pairs[["package_context"]] <- list(
    question = paste0("I want to use the ", fn, "() function. Which package do I need to load?"),
    answer = paste0("You need to load the ", pkg, " package: library(", pkg, ")"),
    package = pkg,
    function_name = fn,
    qa_type = "package_context"
  )

  # 2. Purpose - understanding when to use this function
  qa_pairs[["purpose"]] <- list(
    question = paste0("When should I use ", pkg, "::", fn, "() in my R code?"),
    answer = description,
    package = pkg,
    function_name = fn,
    qa_type = "purpose"
  )

  # 3. Syntax - exact code structure
  if (!is.na(func_data$usage) && nchar(func_data$usage) > 0) {
    qa_pairs[["syntax"]] <- list(
      question = paste0("What is the correct syntax to call ", pkg, "::", fn, "()?"),
      answer = func_data$usage,
      package = pkg,
      function_name = fn,
      qa_type = "syntax"
    )
  }

  # 4. Essential parameters - minimum code to work
  if (!is.na(func_data$sample_arguments) && nchar(func_data$sample_arguments) > 0) {
    qa_pairs[["essential_parameters"]] <- list(
      question = paste0("What are the required arguments I must provide to ", pkg, "::", fn, "()?"),
      answer = func_data$sample_arguments,
      package = pkg,
      function_name = fn,
      qa_type = "essential_parameters"
    )
  }

  # 5. All parameters - complete reference for advanced usage
  if (!is.na(func_data$formatted_arguments) && nchar(func_data$formatted_arguments) > 0) {
    qa_pairs[["all_parameters"]] <- list(
      question = paste0("What parameters can I use with ", pkg, "::", fn, "() to customize its behavior?"),
      answer = func_data$formatted_arguments,
      package = pkg,
      function_name = fn,
      qa_type = "all_parameters"
    )
  }

  # 6. Examples - practical code patterns
  if (!is.na(func_data$examples) && nchar(func_data$examples) > 0) {
    qa_pairs[["examples"]] <- list(
      question = paste0("Show me working R code examples using ", pkg, "::", fn, "()."),
      answer = func_data$examples,
      package = pkg,
      function_name = fn,
      qa_type = "examples"
    )
  }

  # 7. Code completion pattern - helps with autocomplete-like behavior
  if (!is.na(func_data$usage) && nchar(func_data$usage) > 0) {
    # Extract just the function call part for completion
    basic_call <- extract_basic_call(func_data$usage, fn)
    if (!is.na(basic_call)) {
      qa_pairs[["code_completion"]] <- list(
        question = paste0("Complete this R code: ", pkg, "::", fn, "("),
        answer = basic_call,
        package = pkg,
        function_name = fn,
        qa_type = "code_completion"
      )
    }
  }

  # 8. Return value - what to expect from the code
  if (!is.na(func_data$return_value) && nchar(func_data$return_value) > 0) {
    qa_pairs[["return_value"]] <- list(
      question = paste0("What will my R code return when I use ", pkg, "::", fn, "()?"),
      answer = func_data$return_value,
      package = pkg,
      function_name = fn,
      qa_type = "return_value"
    )
  }

  return(qa_pairs)
}


#' @noRd
extract_basic_call <- function(usage_text, func_name) {
  tryCatch(
    {
      # Find the specific function in usage
      pattern <- paste0(func_name, "\\s*\\([^)]*\\)")
      match <- regexpr(pattern, usage_text)
      if (match > 0) {
        call_text <- substr(usage_text, match, match + attr(match, "match.length") - 1)
        # Remove the function name, keep just the parameters
        param_part <- sub(paste0("^", func_name, "\\s*\\("), "", call_text)
        param_part <- sub("\\)$", "", param_part)
        return(param_part)
      }
      return(NA)
    },
    error = function(e) {
      return(NA)
    }
  )
}
