#' Get LLM Response with Text or JSON Output
#'
#' This function sends a prompt to a Language Learning Model (LLM) and returns either
#' a text response or a JSON-structured response based on the provided parameters.
#' It handles retries, validation, and response formatting automatically.
#'
#' @param prompt A character string or tidyprompt object containing the prompt to send to the LLM.
#'   This is the main input that the LLM will respond to.
#'
#' @param llm_client An LLM provider object created by functions like \code{llm_openai()} or
#'   \code{llm_ollama()}. This object contains the configuration for connecting to and
#'   communicating with the specific LLM service.
#'
#' @param max_retries Integer. Maximum number of retry attempts if the LLM fails to provide
#'   a valid response (default: 5). The function will retry if:
#'   - The response doesn't meet validation criteria
#'   - JSON parsing fails (when using json_schema)
#'   - Network or API errors occur
#'   If max_retries is exceeded, NULL is returned.
#'
#' @param max_words Integer or NULL. Maximum number of words allowed in the response
#'   (default: NULL, no limit). Only applies when json_schema is NULL (text responses).
#'   If specified, responses exceeding this limit will trigger a retry.
#'   Example: max_words = 50 limits response to 50 words or fewer.
#'
#' @param max_characters Integer or NULL. Maximum number of characters allowed in the response
#'   (default: NULL, no limit). Only applies when json_schema is NULL (text responses).
#'   If specified, responses exceeding this limit will trigger a retry.
#'   Example: max_characters = 280 limits response to Twitter-like length.
#'
#' @param json_schema List or NULL. JSON schema specification for structured responses
#'   (default: NULL for text responses). When provided, the LLM will be forced to
#'   return a valid JSON object matching the schema. The schema should be a list
#'   representing a JSON schema structure with:
#'   - name: Schema identifier
#'   - description: Schema description
#'   - schema: The actual JSON schema with type, properties, required fields, etc.
#'   Example: list(name = "person", schema = list(type = "object", properties = ...))
#'
#' @param schema_strict Logical. Whether to enforce strict schema validation
#'   (default: FALSE). When TRUE:
#'   - JSON responses must exactly match the schema
#'   - No additional properties are allowed beyond those specified
#'   - All required fields must be present
#'   Only applicable when json_schema is provided.
#'
#' @param schema_type Character. Method for enforcing JSON response format
#'   (default: 'auto'). Options:
#'   - 'auto': Automatically detect best method based on LLM provider
#'   - 'text-based': Add JSON instructions to prompt (works with any provider)
#'   - 'openai': Use OpenAI's native JSON mode (requires compatible OpenAI API)
#'   - 'ollama': Use Ollama's native JSON mode (requires compatible Ollama model)
#'   - 'openai_oo': OpenAI mode without schema enforcement in API
#'   - 'ollama_oo': Ollama mode without schema enforcement in API
#'
#' @param verbose Logical or NULL. Whether to print detailed interaction logs to console
#'   (default: NULL, uses LLM client's setting). When TRUE:
#'   - Shows the prompt being sent
#'   - Displays the LLM's response
#'   - Reports retry attempts and validation failures
#'   Useful for debugging and monitoring LLM interactions.
#'
#' @param stream Logical or NULL. Whether to stream the response in real-time
#'   (default: NULL, uses LLM client's setting). When TRUE:
#'   - Response appears progressively as the LLM generates it
#'   - Provides faster perceived response time
#'   - Only works if the LLM provider supports streaming
#'   Note: Streaming is automatically disabled when verbose = FALSE.
#'
#' @param clean_chat_history Logical. Whether to clean chat history between retries
#'   (default: TRUE). When TRUE:
#'   - Keeps only essential messages in context (first/last user message, last assistant message, system messages)
#'   - Reduces context window usage on retries
#'   - May improve performance with repeatedly failing responses
#'   When FALSE, full conversation history is maintained.
#'
#' @param return_mode Character. What information to return (default: "only_response"). Options:
#'   - "only_response": Returns only the processed LLM response (character string or parsed JSON)
#'   - "full": Returns a comprehensive list containing:
#'     * response: The processed LLM response
#'     * interactions: Number of interactions with the LLM
#'     * chat_history: Complete conversation history
#'     * chat_history_clean: Cleaned conversation history
#'     * start_time: When the function started
#'     * end_time: When the function completed
#'     * duration_seconds: Total execution time
#'     * http_list: Raw HTTP responses from the API
#'
#' @return Depends on return_mode parameter:
#'   - If return_mode = "only_response": Character string (text mode) or parsed list (JSON mode)
#'   - If return_mode = "full": Named list with response and metadata
#'   - NULL if all retry attempts fail
#'
#' @details
#' This function serves as a unified interface for getting responses from LLMs with
#' automatic handling of different response formats and validation. It internally uses
#' the tidyprompt package's \code{answer_as_text()} or \code{answer_as_json()} functions
#' depending on whether a JSON schema is provided.
#'
#' **Text Mode (json_schema = NULL):**
#' - Uses \code{answer_as_text()} with optional word/character limits
#' - Returns plain text responses
#' - Validates response length constraints
#'
#' **JSON Mode (json_schema provided):**
#' - Uses \code{answer_as_json()} with schema validation
#' - Forces structured JSON responses
#' - Validates against provided schema
#' - Returns parsed R objects (lists)
#'
#' **Error Handling:**
#' The function automatically retries on various failure conditions including
#' validation errors, JSON parsing errors, and network issues.
#'
#' @examples
#' \dontrun{
#' # Basic text response
#' client <- llm_ollama()
#' response <- get_llm_response("What is R?", client)
#'
#' # Text response with word limit
#' short_response <- get_llm_response(
#'   "Explain machine learning",
#'   client,
#'   max_words = 50
#' )
#'
#' # JSON response with schema
#' schema <- list(
#'   name = "person_info",
#'   schema = list(
#'     type = "object",
#'     properties = list(
#'       name = list(type = "string"),
#'       age = list(type = "integer")
#'     ),
#'     required = c("name", "age")
#'   )
#' )
#'
#' json_response <- get_llm_response(
#'   "Create a person with name and age",
#'   client,
#'   json_schema = schema
#' )
#'
#' # Full response with metadata
#' full_result <- get_llm_response(
#'   "Hello",
#'   client,
#'   return_mode = "full",
#'   verbose = TRUE
#' )
#' }
#'
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
get_llm_response <- function(prompt,
                             llm_client,
                             max_retries = 5,
                             max_words = NULL,
                             max_characters = NULL,
                             json_schema = NULL,
                             schema_strict = FALSE,
                             schema_type = "auto",
                             verbose = NULL,
                             stream = NULL,
                             clean_chat_history = TRUE,
                             return_mode = c("only_response", "full")) {
  # Validate inputs
  return_mode <- match.arg(return_mode)

  if (!is.null(max_retries) && (!is.numeric(max_retries) || max_retries < 1)) {
    stop("max_retries must be a positive integer")
  }

  if (!is.null(max_words) && (!is.numeric(max_words) || max_words < 1)) {
    stop("max_words must be a positive integer")
  }

  if (!is.null(max_characters) && (!is.numeric(max_characters) || max_characters < 1)) {
    stop("max_characters must be a positive integer")
  }

  # Route to appropriate response handler
  if (is.null(json_schema)) {
    # Text response mode
    res <- prompt |>
      tidyprompt::answer_as_text(
        max_words = max_words,
        max_characters = max_characters
      ) |>
      tidyprompt::send_prompt(
        llm_provider = llm_client,
        max_interactions = max_retries,
        verbose = verbose,
        stream = stream,
        clean_chat_history = clean_chat_history,
        return_mode = return_mode
      )
  } else {
    # JSON response mode
    res <- prompt |>
      tidyprompt::answer_as_json(
        schema = json_schema,
        schema_strict = schema_strict,
        type = schema_type
      ) |>
      tidyprompt::send_prompt(
        llm_provider = llm_client,
        max_interactions = max_retries,
        verbose = verbose,
        stream = stream,
        clean_chat_history = clean_chat_history,
        return_mode = return_mode
      )
  }

  return(res)
}
