#' List available models from Ollama API
#'
#' This function retrieves information about available models from the Ollama API
#' and returns it as a tibble with simplified data extraction.
#'
#' @param .ollama_server The URL of the Ollama server (default: "http://localhost:11434")
#' @return A tibble containing model information, or NULL if no models are found
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
ollama_list_models <- function(.ollama_server = "http://localhost:11434") {
  response <- httr2::request(.ollama_server) |>
    httr2::req_url_path("/api/tags") |>
    httr2::req_perform() |>
    httr2::resp_body_json()

  if (is.null(response$models)) {
    return(NULL)
  }

  # Use map_dfr for simplified data extraction
  purrr::map_dfr(response$models, ~ {
    tibble::tibble(
      name = .x$name,
      modified_at = .x$modified_at,
      size = .x$size,
      format = .x$details$format,
      family = .x$details$family,
      parameter_size = .x$details$parameter_size,
      quantization_level = .x$details$quantization_level
    )
  })
}

#' Download a model from Ollama API
#'
#' This function sends a request to download a specified model from Ollama's
#' model library with progress tracking.
#'
#' @param .model The name of the model to download
#' @param .ollama_server The URL of the Ollama server (default: "http://localhost:11434")
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
ollama_download_model <- function(.model, .ollama_server = "http://localhost:11434") {
  progress_bar <- cli::cli_progress_bar(auto_terminate = FALSE, type = "download")

  # Simplified stream callback function
  stream_callback <- function(.stream) {
    line <- trimws(rawToChar(.stream))
    if (nchar(line) == 0) {
      return(TRUE)
    }

    data <- tryCatch(jsonlite::fromJSON(line), error = function(e) NULL)
    if (is.null(data)) {
      return(TRUE)
    }

    # Update progress bar
    if (grepl("pulling", data$status) && !is.null(data$total) && !is.null(data$completed)) {
      if (data$total > 0 && data$completed >= 0) {
        cli::cli_progress_update(set = data$completed, total = data$total, id = progress_bar)
      }
    }
    TRUE
  }

  httr2::request(.ollama_server) |>
    httr2::req_url_path("/api/pull") |>
    httr2::req_body_json(list(name = .model)) |>
    httr2::req_perform_stream(stream_callback, buffer_kb = 0.05, round = "line")

  cli::cli_progress_done(id = progress_bar)
  invisible(NULL)
}

#' Delete a model from Ollama API
#'
#' This function sends a DELETE request to remove a specified model from
#' the Ollama API and returns the updated model list.
#'
#' @param .model The name of the model to delete
#' @param .ollama_server The URL of the Ollama server (default: "http://localhost:11434")
#' @return Updated tibble of available models after deletion
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
ollama_delete_model <- function(.model, .ollama_server = "http://localhost:11434") {
  httr2::request(.ollama_server) |>
    httr2::req_url_path("/api/delete") |>
    httr2::req_method("DELETE") |>
    httr2::req_body_json(list(name = .model)) |>
    httr2::req_error(body = function(resp) httr2::resp_body_json(resp)$error) |>
    httr2::req_perform()

  cli::cli_alert_success(paste("Model", .model, "has been deleted"))
  ollama_list_models(.ollama_server)
}
