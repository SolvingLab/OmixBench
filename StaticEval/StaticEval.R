library(tidyverse)
library(readxl)
library(llmhelper)
library(OmixBenchR)
library(glue)

# ============================================================================
# Load Configuration
# ============================================================================

load('models_config.rda')
tasks_df <- read_excel("OmixTask1002.xlsx")

dir.create("results", showWarnings = FALSE)

# ============================================================================
# Helper Functions
# ============================================================================

create_client <- function(model_name, type, provider) {
  if (type == "API") {
    env_prefix <- toupper(gsub("-", "_", provider))
    base_url <- Sys.getenv(paste0(env_prefix, "_BASE_URL"))
    api_key <- Sys.getenv(paste0(env_prefix, "_API_KEY"))
    
    if (base_url == "" || api_key == "") return(NULL)
    
    llm_provider(
      base_url = base_url,
      api_key = api_key,
      model = model_name,
      temperature = 0.2,
      max_tokens = 15000,
      timeout = 500,
      stream = FALSE,
      verbose = FALSE
    )
  } else {
    llm_ollama(
      model = model_name,
      temperature = 0.2,
      max_tokens = 15000,
      timeout = 500,
      stream = FALSE,
      verbose = FALSE,
      skip_test = FALSE,
      auto_download = TRUE
    )
  }
}

generate_response <- function(task, client) {
  tryCatch({
    get_llm_response(
      set_prompt(
        system = "You are an AI assistant specialized in bioinformatics. Provide complete R code solutions.",
        user = glue("Please provide R code to solve the following bioinformatics task:\n\n{task}")
      ),
      client,
      max_retries = 1,
      verbose = FALSE
    )
  }, error = function(e) NA_character_)
}

# ============================================================================
# Create Evaluators
# ============================================================================

eval_ds <- tryCatch({
  llm_provider(
    base_url = Sys.getenv('DEEPSEEK_BASE_URL'),
    api_key = Sys.getenv('DEEPSEEK_API_KEY'),
    model = 'deepseek-r1-250528',
    temperature = 0.2,
    timeout = 500,
    stream = FALSE,
    verbose = FALSE
  )
}, error = function(e) NULL)

eval_o3 <- tryCatch({
  llm_provider(
    base_url = Sys.getenv('OPENAI_BASE_URL'),
    api_key = Sys.getenv('OPENAI_API_KEY'),
    model = 'o3',
    temperature = 0.2,
    timeout = 500,
    stream = FALSE,
    verbose = FALSE
  )
}, error = function(e) NULL)

eval_claude <- tryCatch({
  llm_provider(
    base_url = Sys.getenv('CLAUDE_BASE_URL'),
    api_key = Sys.getenv('CLAUDE_API_KEY'),
    model = 'claude-sonnet-4-20250514',
    temperature = 0.2,
    timeout = 500,
    stream = FALSE,
    verbose = FALSE
  )
}, error = function(e) NULL)

evaluators <- list(
  deepseek = eval_ds,
  o3 = eval_o3,
  claude = eval_claude
) %>% compact()

cat(glue("Initialized {length(evaluators)} evaluators\n"))

# ============================================================================
# Main Loop
# ============================================================================

for (i in seq_len(nrow(models_config))) {
  
  model_name <- models_config$model_name[i]
  model_type <- models_config$type[i]
  provider <- models_config$provider[i]
  
  cat(glue("\n{'='*60}\n"))
  cat(glue("Model {i}/{nrow(models_config)}: {model_name}\n"))
  cat(glue("{'='*60}\n"))
  
  # File paths
  gen_file <- glue("results/{model_name}_results.rds")
  eval_file <- glue("results/{model_name}_evaluation.rds")
  
  # ============================================================================
  # Generation
  # ============================================================================
  
  if (file.exists(gen_file)) {
    cat("Loading existing generation results...\n")
    results <- readRDS(gen_file)
  } else {
    cat("Starting generation...\n")
    
    client <- create_client(model_name, model_type, provider)
    if (is.null(client)) {
      cat("Failed to initialize client, skipping...\n")
      next
    }
    
    results <- tasks_df %>%
      mutate(
        task_id = row_number(),
        response = map_chr(task, ~generate_response(.x, client)),
        model_name = model_name,
        timestamp = as.character(Sys.time())
      )
    
    saveRDS(results, gen_file)
    cat(glue("Saved generation results: {sum(!is.na(results$response))}/{nrow(results)}\n"))
  }
  
  # ============================================================================
  # Evaluation
  # ============================================================================
  
  if (file.exists(eval_file)) {
    cat("Evaluation already exists, skipping...\n")
    next
  }
  
  if (length(evaluators) == 0) {
    cat("No evaluators available, skipping evaluation...\n")
    next
  }
  
  cat("Starting evaluation...\n")
  
  # Store all evaluation results
  all_eval_results <- list()
  
  for (j in seq_len(nrow(results))) {
    task <- results$task[j]
    task_id <- results$task_id[j]
    response <- results$response[j]
    
    if (is.na(response)) next
    
    # Evaluate with each evaluator
    for (eval_name in names(evaluators)) {
      eval_client <- evaluators[[eval_name]]
      
      eval_result <- tryCatch({
        res <- CodeEval(task, response, eval_client, 
                        schema_type = 'text-based', verbose = FALSE)
        
        # Store complete result
        list(
          dataset = "OmixTask1002",
          model = eval_name,
          task_id = task_id,
          total_score = res$total_score %||% NA,
          breakdown = res$breakdown,
          success = !is.null(res$total_score)
        )
      }, error = function(e) {
        list(
          dataset = "OmixTask1002",
          model = eval_name,
          task_id = task_id,
          total_score = NA,
          breakdown = list(),
          success = FALSE
        )
      })
      
      all_eval_results <- append(all_eval_results, list(eval_result))
    }
    
    if (j %% 10 == 0) {
      cat(glue("Evaluated {j}/{nrow(results)} tasks\n"))
    }
  }
  
  # Convert to dataframe with proper structure
  eval_df <- map_dfr(all_eval_results, function(x) {
    tibble(
      dataset      = x$dataset      %||% NA,
      model        = x$model        %||% NA,
      task_id      = x$task_id      %||% NA,
      total_score  = x$total_score  %||% NA,
      
      problem_solving          = pluck(x, "breakdown", "problem_solving",        .default = NA),
      technical_implementation = pluck(x, "breakdown", "technical_implementation", .default = NA),
      data_handling            = pluck(x, "breakdown", "data_handling",          .default = NA),
      statistical_rigor        = pluck(x, "breakdown", "statistical_rigor",      .default = NA),
      domain_knowledge         = pluck(x, "breakdown", "domain_knowledge",       .default = NA),
      robustness               = pluck(x, "breakdown", "robustness",             .default = NA),
      documentation            = pluck(x, "breakdown", "documentation",          .default = NA),
      
      success = x$success %||% FALSE
    )
  })
  
  # Filter out invalid scores
  eval_df <- eval_df %>%
    filter(
      total_score <= 100,
      problem_solving <= 20,
      technical_implementation <= 25,
      data_handling <= 15,
      statistical_rigor <= 15,
      domain_knowledge <= 10,
      robustness <= 10,
      documentation <= 5
    )
  
  # Calculate average scores across evaluators
  eval_summary <- eval_df %>%
    group_by(task_id) %>%
    summarise(
      n_evaluators = n(),
      avg_total_score = mean(total_score, na.rm = TRUE),
      avg_problem_solving = mean(problem_solving, na.rm = TRUE),
      avg_technical_implementation = mean(technical_implementation, na.rm = TRUE),
      avg_data_handling = mean(data_handling, na.rm = TRUE),
      avg_statistical_rigor = mean(statistical_rigor, na.rm = TRUE),
      avg_domain_knowledge = mean(domain_knowledge, na.rm = TRUE),
      avg_robustness = mean(robustness, na.rm = TRUE),
      avg_documentation = mean(documentation, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Combine with original results
  final_results <- results %>%
    left_join(eval_summary, by = "task_id")
  
  # Save both detailed and summary results
  saveRDS(list(
    detailed = eval_df,
    summary = eval_summary,
    full = final_results
  ), eval_file)
  
  cat(glue("Saved evaluation results\n"))
  cat(glue("Valid evaluations: {nrow(eval_df)}/{length(all_eval_results)}\n"))
  cat(glue("Average score: {round(mean(eval_summary$avg_total_score, na.rm=TRUE), 2)}\n"))
}
