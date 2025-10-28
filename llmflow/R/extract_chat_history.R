#' Extract chat history from ellmer chat object
#'
#' @param chat_obj An ellmer chat object
#' @param include_tokens Whether to include token information
#' @param include_time Whether to include timestamp information
#' @param tz Time zone for timestamps (default "Asia/Shanghai" for CST)
#' @return A data frame with chat history
#' @export
extract_chat_history <- function(chat_obj,
                                 include_tokens = TRUE,
                                 include_time = TRUE,
                                 tz = "Asia/Shanghai") {
  # Get turns
  turns <- chat_obj$get_turns()

  if (length(turns) == 0) {
    return(data.frame())
  }

  # Extract information from each turn
  history_list <- list()

  for (i in seq_along(turns)) {
    turn <- turns[[i]]

    # Start with basic info
    row_data <- data.frame(
      turn = i,
      role = turn@role,
      text = turn@text,
      stringsAsFactors = FALSE
    )

    # Add tokens if requested
    if (include_tokens) {
      if (length(turn@tokens) >= 3) {
        row_data$tokens_input <- turn@tokens[1]
        row_data$tokens_output <- turn@tokens[2]
        row_data$tokens_cached <- turn@tokens[3]

        # Calculate total tokens for assistant responses
        if (turn@role == "assistant") {
          row_data$tokens_total <- turn@tokens[1] + turn@tokens[2] + turn@tokens[3]
        } else {
          row_data$tokens_total <- NA
        }
      } else {
        row_data$tokens_input <- NA
        row_data$tokens_output <- NA
        row_data$tokens_cached <- NA
        row_data$tokens_total <- NA
      }
    }

    # Add timestamp if requested
    if (include_time) {
      if (length(turn@json) > 0 && !is.null(turn@json$created)) {
        row_data$timestamp <- as.POSIXct(turn@json$created,
          origin = "1970-01-01",
          tz = tz
        )
      } else {
        row_data$timestamp <- as.POSIXct(NA)
      }
    }

    history_list[[i]] <- row_data
  }

  # Combine all rows
  df <- do.call(rbind, history_list)

  # Calculate time differences between consecutive assistant responses
  if (include_time && "timestamp" %in% names(df)) {
    df$time_since_last <- NA

    assistant_idx <- which(df$role == "assistant" & !is.na(df$timestamp))

    if (length(assistant_idx) > 1) {
      for (i in 2:length(assistant_idx)) {
        idx <- assistant_idx[i]
        prev_idx <- assistant_idx[i - 1]
        time_diff <- as.numeric(difftime(df$timestamp[idx],
          df$timestamp[prev_idx],
          units = "secs"
        ))
        df$time_since_last[idx] <- round(time_diff, 1)
      }
    }
  }

  # Fix user tokens showing as 0 instead of NA
  if (include_tokens) {
    user_rows <- df$role == "user"
    df$tokens_input[user_rows & df$tokens_input == 0] <- NA
    df$tokens_output[user_rows & df$tokens_output == 0] <- NA
    df$tokens_cached[user_rows & df$tokens_cached == 0] <- NA
  }

  return(df)
}
