#' Set system and user prompts for LLM interaction
#'
#' This function creates a prompt object with system and user prompts
#' using the tidyprompt package for structured LLM communication.
#'
#' @param system The system prompt to set context and behavior (default: bioinformatics assistant)
#' @param user The user prompt or question
#' @return A prompt object configured with system and user prompts
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @export
set_prompt <- function(system = "You are an AI assistant specialized in bioinformatics.",
                       user = "Hi") {
  prompt <- tidyprompt::set_system_prompt(
    prompt = user,
    system_prompt = system
  )

  return(prompt)
}
