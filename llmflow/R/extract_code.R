#' Extract R code from a string
#'
#' This function extracts R code from a string by matching all content between
#' '```r' or '```R' and '```'.
#'
#' @param input_string A string containing R code blocks, typically a response
#' from an LLM
#'
#' @return A character vector containing the extracted R code
#'
#' @export
#'
#' @examples
#' # Simple example
#' text <- "Here is some R code:\n```r\nprint('Hello')\n```"
#' extract_r_code(text)
#'
#' # Multiple code blocks
#' response <- "
#' First block:
#' ```r
#' x <- 1:10
#' mean(x)
#' ```
#'
#' Second block:
#' ```R
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, hp)) + geom_point()
#' ```
#' "
#' codes <- extract_r_code(response)
#' length(codes) # Returns 2
#'
#' # With surrounding text
#' llm_response <- "
#' To calculate the mean, use this code:
#' ```r
#' data <- c(1, 2, 3, 4, 5)
#' result <- mean(data)
#' print(result)
#' ```
#' The result will be 3.
#' "
#' extract_r_code(llm_response)
#'
extract_r_code <- function(input_string) {
  matches <- gregexpr("(?s)```[rR]\\s*(.*?)\\s*```", input_string, perl = TRUE)
  extracted_code <- regmatches(input_string, matches)

  extracted_code <- lapply(extracted_code, function(x) {
    sub("(?s)```[rR]\\s*(.*?)\\s*```", "\\1", x, perl = TRUE)
  })

  return(unlist(extracted_code))
}

#' Extract Python code from a string
#'
#' This function extracts Python code from a string by matching all content between
#' '```python', '```py' and '```'.
#'
#' @param input_string A string containing Python code blocks, typically a response
#' from an LLM
#'
#' @return A character vector containing the extracted Python code
#'
#' @export
#'
#' @examples
#' # Simple example
#' text <- "Python code:\n```python\nprint('Hello World')\n```"
#' extract_python_code(text)
#'
#' # Using 'py' tag
#' text <- "```py\nimport numpy as np\n```"
#' extract_python_code(text)
#'
#' # Multiple blocks with different tags
#' response <- "
#' Data processing:
#' ```python
#' import pandas as pd
#' df = pd.read_csv('data.csv')
#' df.head()
#' ```
#'
#' Visualization:
#' ```py
#' import matplotlib.pyplot as plt
#' plt.plot([1, 2, 3], [4, 5, 6])
#' plt.show()
#' ```
#' "
#' codes <- extract_python_code(response)
#' length(codes) # Returns 2
#'
#' # Complex example with classes and functions
#' llm_response <- "
#' Here's a complete Python solution:
#' ```python
#' class DataProcessor:
#'     def __init__(self, data):
#'         self.data = data
#'
#'     def process(self):
#'         return [x * 2 for x in self.data]
#'
#' processor = DataProcessor([1, 2, 3])
#' result = processor.process()
#' print(result)
#' ```
#' "
#' extract_python_code(llm_response)
#'
extract_python_code <- function(input_string) {
  matches <- gregexpr("(?s)```(?:python|py)\\s*(.*?)\\s*```", input_string, perl = TRUE)
  extracted_code <- regmatches(input_string, matches)

  extracted_code <- lapply(extracted_code, function(x) {
    sub("(?s)```(?:python|py)\\s*(.*?)\\s*```", "\\1", x, perl = TRUE)
  })

  return(unlist(extracted_code))
}

#' Extract Bash/Shell code from a string
#'
#' This function extracts Bash/Shell code from a string by matching all content between
#' '```bash', '```sh', '```shell' and '```'.
#'
#' @param input_string A string containing Bash/Shell code blocks, typically a response
#' from an LLM
#'
#' @return A character vector containing the extracted Bash/Shell code
#'
#' @export
#'
#' @examples
#' # Simple bash example
#' text <- "Run this:\n```bash\necho 'Hello'\n```"
#' extract_bash_code(text)
#'
#' # Using 'sh' tag
#' text <- "```sh\nls -la\npwd\n```"
#' extract_bash_code(text)
#'
#' # Using 'shell' tag
#' text <- "```shell\nfor i in {1..5}; do echo $i; done\n```"
#' extract_bash_code(text)
#'
#' # Multiple blocks with different tags
#' response <- "
#' Setup script:
#' ```bash
#' #!/bin/bash
#' mkdir -p /tmp/test
#' cd /tmp/test
#' ```
#'
#' Installation:
#' ```sh
#' apt-get update
#' apt-get install -y git
#' ```
#'
#' Configuration:
#' ```shell
#' export PATH=$PATH:/usr/local/bin
#' source ~/.bashrc
#' ```
#' "
#' codes <- extract_bash_code(response)
#' length(codes) # Returns 3
#'
#' # Complex script example
#' script_response <- "
#' Here's a backup script:
#' ```bash
#' #!/bin/bash
#'
#' # Set variables
#' BACKUP_DIR='/backup'
#' DATE=$(date +%Y%m%d)
#'
#' # Create backup
#' tar -czf ${BACKUP_DIR}/backup_${DATE}.tar.gz /home/user/
#'
#' # Check if successful
#' if [ $? -eq 0 ]; then
#'     echo 'Backup completed successfully'
#' else
#'     echo 'Backup failed'
#'     exit 1
#' fi
#' ```
#' "
#' extract_bash_code(script_response)
#'
extract_bash_code <- function(input_string) {
  matches <- gregexpr("(?s)```(?:bash|sh|shell)\\s*(.*?)\\s*```", input_string, perl = TRUE)
  extracted_code <- regmatches(input_string, matches)

  extracted_code <- lapply(extracted_code, function(x) {
    sub("(?s)```(?:bash|sh|shell)\\s*(.*?)\\s*```", "\\1", x, perl = TRUE)
  })

  return(unlist(extracted_code))
}

#' Extract SQL code from a string
#'
#' This function extracts SQL code from a string by matching all content between
#' '```sql' and '```' (case-insensitive).
#'
#' @param input_string A string containing SQL code blocks, typically a response
#' from an LLM
#'
#' @return A character vector containing the extracted SQL code
#'
#' @export
#'
#' @examples
#' # Simple SQL query
#' text <- "Query:\n```sql\nSELECT * FROM users;\n```"
#' extract_sql_code(text)
#'
#' # Case-insensitive matching
#' text <- "```SQL\nSELECT COUNT(*) FROM orders;\n```"
#' extract_sql_code(text)
#'
#' # Multiple SQL blocks
#' response <- "
#' Create table:
#' ```sql
#' CREATE TABLE employees (
#'     id INT PRIMARY KEY,
#'     name VARCHAR(100),
#'     department VARCHAR(50),
#'     salary DECIMAL(10, 2)
#' );
#' ```
#'
#' Insert data:
#' ```sql
#' INSERT INTO employees (id, name, department, salary)
#' VALUES
#'     (1, 'John Doe', 'IT', 75000),
#'     (2, 'Jane Smith', 'HR', 65000);
#' ```
#'
#' Query data:
#' ```sql
#' SELECT name, salary
#' FROM employees
#' WHERE department = 'IT'
#' ORDER BY salary DESC;
#' ```
#' "
#' codes <- extract_sql_code(response)
#' length(codes) # Returns 3
#'
#' # Complex query with joins
#' complex_query <- "
#' Here's the analysis query:
#' ```sql
#' WITH monthly_sales AS (
#'     SELECT
#'         DATE_TRUNC('month', order_date) as month,
#'         SUM(total_amount) as total_sales,
#'         COUNT(DISTINCT customer_id) as unique_customers
#'     FROM orders
#'     WHERE order_date >= '2024-01-01'
#'     GROUP BY DATE_TRUNC('month', order_date)
#' )
#' SELECT
#'     month,
#'     total_sales,
#'     unique_customers,
#'     total_sales / unique_customers as avg_per_customer
#' FROM monthly_sales
#' ORDER BY month;
#' ```
#' "
#' extract_sql_code(complex_query)
#'
extract_sql_code <- function(input_string) {
  matches <- gregexpr("(?si)```sql\\s*(.*?)\\s*```", input_string, perl = TRUE)
  extracted_code <- regmatches(input_string, matches)

  extracted_code <- lapply(extracted_code, function(x) {
    sub("(?si)```sql\\s*(.*?)\\s*```", "\\1", x, perl = TRUE)
  })

  return(unlist(extracted_code))
}

#' Extract JavaScript code from a string
#'
#' This function extracts JavaScript code from a string by matching all content between
#' '```javascript', '```js', '```jsx' and '```'.
#'
#' @param input_string A string containing JavaScript code blocks, typically a response
#' from an LLM
#'
#' @return A character vector containing the extracted JavaScript code
#'
#' @export
#'
#' @examples
#' # Simple JavaScript example
#' text <- "Code:\n```javascript\nconsole.log('Hello');\n```"
#' extract_javascript_code(text)
#'
#' # Using 'js' tag
#' text <- "```js\nconst x = 42;\n```"
#' extract_javascript_code(text)
#'
#' # Using 'jsx' tag for React
#' text <- "```jsx\n<div>Hello World</div>\n```"
#' extract_javascript_code(text)
#'
#' # Multiple blocks with different tags
#' response <- "
#' Frontend code:
#' ```javascript
#' function fetchData() {
#'     return fetch('/api/data')
#'         .then(response => response.json());
#' }
#' ```
#'
#' React component:
#' ```jsx
#' const MyComponent = () => {
#'     const [data, setData] = useState([]);
#'
#'     useEffect(() => {
#'         fetchData().then(setData);
#'     }, []);
#'
#'     return (
#'         <div>
#'             {data.map(item => <p key={item.id}>{item.name}</p>)}
#'         </div>
#'     );
#' };
#' ```
#'
#' Node.js backend:
#' ```js
#' const express = require('express');
#' const app = express();
#'
#' app.get('/api/data', (req, res) => {
#'     res.json([{id: 1, name: 'Item 1'}]);
#' });
#'
#' app.listen(3000);
#' ```
#' "
#' codes <- extract_javascript_code(response)
#' length(codes) # Returns 3
#'
extract_javascript_code <- function(input_string) {
  matches <- gregexpr("(?s)```(?:javascript|js|jsx)\\s*(.*?)\\s*```", input_string, perl = TRUE)
  extracted_code <- regmatches(input_string, matches)

  extracted_code <- lapply(extracted_code, function(x) {
    sub("(?s)```(?:javascript|js|jsx)\\s*(.*?)\\s*```", "\\1", x, perl = TRUE)
  })

  return(unlist(extracted_code))
}

#' Generic function to extract code of any specified language
#'
#' This function provides a flexible way to extract code blocks of any language
#' from a string by specifying the language identifier(s).
#'
#' @param input_string A string containing code blocks
#' @param language Language identifier(s) to extract (e.g., "r", "python", c("bash", "sh"))
#' @param case_sensitive Whether the language matching should be case-sensitive (default: FALSE)
#'
#' @return A character vector containing the extracted code
#'
#' @export
#'
#' @examples
#' # Extract R code
#' text <- "```r\nx <- 1:10\n```"
#' extract_code(text, "r")
#'
#' # Extract multiple language variants
#' text <- "```bash\necho 'test'\n```\n```sh\nls -la\n```"
#' extract_code(text, c("bash", "sh"))
#'
#' # Case-sensitive extraction
#' text <- "```R\nplot(1:10)\n```\n```r\nprint('hello')\n```"
#' extract_code(text, "r", case_sensitive = TRUE) # Only matches lowercase 'r'
#' extract_code(text, "r", case_sensitive = FALSE) # Matches both 'R' and 'r'
#'
#' # Extract custom language
#' text <- "```julia\nprintln(\"Julia code\")\n```"
#' extract_code(text, "julia")
#'
#' # Extract YAML configuration
#' config_text <- "
#' Here's the configuration:
#' ```yaml
#' database:
#'   host: localhost
#'   port: 5432
#'   name: mydb
#' ```
#' "
#' extract_code(config_text, "yaml")
#'
#' # Extract multiple TypeScript and JavaScript blocks
#' mixed_text <- "
#' TypeScript:
#' ```typescript
#' interface User {
#'     name: string;
#'     age: number;
#' }
#' ```
#'
#' JavaScript:
#' ```js
#' const user = {name: 'John', age: 30};
#' ```
#' "
#' # Extract TypeScript
#' extract_code(mixed_text, "typescript")
#' # Extract both TypeScript and JavaScript
#' extract_code(mixed_text, c("typescript", "js"))
#'
extract_code <- function(input_string, language, case_sensitive = FALSE) {
  # Build language pattern
  if (length(language) > 1) {
    lang_pattern <- paste0("(?:", paste(language, collapse = "|"), ")")
  } else {
    lang_pattern <- language
  }

  # Build complete regex pattern
  if (case_sensitive) {
    pattern <- sprintf("(?s)```%s\\s*(.*?)\\s*```", lang_pattern)
  } else {
    pattern <- sprintf("(?si)```%s\\s*(.*?)\\s*```", lang_pattern)
  }

  matches <- gregexpr(pattern, input_string, perl = TRUE)
  extracted_code <- regmatches(input_string, matches)

  # Extract code content
  extracted_code <- lapply(extracted_code, function(x) {
    if (case_sensitive) {
      sub(sprintf("(?s)```%s\\s*(.*?)\\s*```", lang_pattern), "\\1", x, perl = TRUE)
    } else {
      sub(sprintf("(?si)```%s\\s*(.*?)\\s*```", lang_pattern), "\\1", x, perl = TRUE)
    }
  })

  return(unlist(extracted_code))
}

#' Save extracted code to file
#'
#' This function saves extracted code to a file with appropriate extension
#' based on the programming language.
#'
#' @param code_string String or character vector containing the code to save
#' @param filename Output filename. If NULL, generates a timestamped filename
#' @param language Programming language for determining file extension (default: "r")
#'
#' @return The path to the saved file
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract and save R code
#' llm_response <- "```r\nplot(1:10)\n```"
#' code <- extract_r_code(llm_response)
#' save_code_to_file(code) # Saves as "code_20240101_120000.R"
#'
#' # Save with custom filename
#' save_code_to_file(code, "my_plot.R")
#'
#' # Save Python code with auto extension
#' py_code <- "import pandas as pd\ndf = pd.DataFrame()"
#' save_code_to_file(py_code, language = "python") # Creates .py file
#'
#' # Save multiple code blocks
#' response <- "```r\nx <- 1\n```\n```r\ny <- 2\n```"
#' codes <- extract_r_code(response)
#' save_code_to_file(codes, "combined_code.R")
#' }
#'
save_code_to_file <- function(code_string,
                              filename = NULL,
                              language = "r") {
  # Define file extensions
  extensions <- list(
    r = ".R",
    python = ".py",
    bash = ".sh",
    sh = ".sh",
    shell = ".sh",
    sql = ".sql",
    javascript = ".js",
    js = ".js",
    jsx = ".jsx",
    typescript = ".ts",
    tsx = ".tsx",
    yaml = ".yml",
    json = ".json",
    html = ".html",
    css = ".css",
    cpp = ".cpp",
    c = ".c",
    java = ".java",
    go = ".go",
    rust = ".rs",
    julia = ".jl"
  )

  # Get extension
  ext <- extensions[[tolower(language)]]
  if (is.null(ext)) ext <- ".txt"

  # Generate filename if not provided
  if (is.null(filename)) {
    filename <- paste0(
      "code_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ext
    )
  } else {
    # Add extension if not present
    if (!grepl("\\.[a-zA-Z]+$", filename)) {
      filename <- paste0(filename, ext)
    }
  }

  # Combine multiple code blocks if necessary
  if (length(code_string) > 1) {
    code_string <- paste(code_string, collapse = "\n\n")
  }

  # Write to file
  writeLines(code_string, filename)

  # Return the filename
  cat("Code saved to:", filename, "\n")
  return(invisible(filename))
}
