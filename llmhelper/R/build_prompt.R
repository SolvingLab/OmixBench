#' Build a templated prompt for LLM interaction using glue
#'
#' This function constructs a structured prompt string by injecting user-supplied
#' parameters into a predefined template. It leverages the glue package to replace
#' named placeholders in the template with actual values, enabling dynamic prompt
#' creation for LLM workflows.
#'
#' @param template A character string containing the prompt template. Placeholders
#'   should be wrapped in `{}` and correspond to names provided in `...`.
#' @param ... Named arguments matching placeholders in `template`. Each name–value
#'   pair will be substituted into the template at runtime.
#'
#' @return A single character string with all `{placeholder}` fields in `template`
#'   replaced by the corresponding values from `...`.
#'
#' @details
#' The `build_prompt()` function uses `glue::glue_data()` internally. Placeholders
#' in `template` (e.g., `{filename}`, `{threshold}`) are resolved by passing a named
#' list of parameters via `...`. You can include any number of placeholders in
#' the template, as long as the corresponding argument is supplied when calling
#' this function.
#'
#' @examples
#' \dontrun{
#' # Define a template with placeholders
#' prompt_template <- "
#' Perform the following analysis on dataset at '{filepath}':
#' 1. Load data from '{filepath}'
#' 2. Normalize using method '{norm_method}'
#' 3. Save results to '{output_dir}'
#'
#' IMPORTANT: Use package::function notation for all function calls."
#'
#' # Build the prompt by supplying named arguments
#' filled_prompt <- build_prompt(
#'   template     = prompt_template,
#'   filepath     = "/path/to/data.csv",
#'   norm_method  = "quantile",
#'   output_dir   = "/path/to/output/"
#' )
#' cat(filled_prompt)
#' }
#'
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
build_prompt <- function(template, ...) {
  args <- list(...)
  glue::glue_data(.x = args, template)
}
