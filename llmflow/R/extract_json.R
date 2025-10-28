#' Extract and parse JSONs from a string (LLM response)
#'
#' This function extracts JSON blocks from a string and parses them using `jsonlite::fromJSON()`.
#' This can be used to extract all JSONs from LLM responses, immediately converting
#' them to R objects.
#'
#' @param llm_response A character string
#'
#' @return A list of parsed JSON objects
#'
#' @export
extract_json <- function(llm_response) {
  positions <- gregexpr("[{}]", llm_response)[[1]]
  matches <- regmatches(llm_response, gregexpr("[{}]", llm_response))[[1]]

  stack <- c()
  blocks <- list()
  block_start <- c()

  for (i in seq_along(positions)) {
    char <- matches[i]
    pos <- positions[i]

    if (is.na(char)) next

    if (char == "{") {
      if (length(stack) == 0) {
        # Start of a new block
        block_start <- c(block_start, pos)
      }
      stack <- c(stack, char)
    } else if (char == "}") {
      if (length(stack) == 0) {
        # Unmatched closing brace
        next
      }
      stack <- stack[-length(stack)]
      if (length(stack) == 0) {
        # End of a block
        start_pos <- block_start[1]
        end_pos <- pos
        block_text <- substr(llm_response, start_pos, end_pos)
        blocks <- c(blocks, block_text)
        block_start <- block_start[-1]
      }
    }
  }

  # Now parse the blocks
  parsed_jsons <- lapply(blocks, function(json_candidate) {
    tryCatch(
      {
        jsonlite::fromJSON(json_candidate, simplifyDataFrame = FALSE)
      },
      error = function(e) NULL
    )
  })

  # Remove NULL entries
  parsed_jsons <- parsed_jsons[!sapply(parsed_jsons, is.null)]
  return(parsed_jsons)
}
