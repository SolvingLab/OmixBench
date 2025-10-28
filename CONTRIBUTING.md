# Contributing to OmixBench

Thank you for your interest in contributing to OmixBench! This document provides guidelines for contributing to the project.

## Code of Conduct

This project and everyone participating in it is governed by our [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to liuzaoqu@163.com.

## Security

If you discover a security vulnerability, please follow the guidelines in our [Security Policy](SECURITY.md). Do not report security vulnerabilities through public GitHub issues.

## Table of Contents

- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Package Structure](#package-structure)
- [Testing](#testing)
- [Types of Contributions](#types-of-contributions)
- [Benchmark Tasks](#benchmark-tasks)
- [Code Review Process](#code-review-process)
- [Communication](#communication)
- [Recognition](#recognition)

## How to Contribute

### Reporting Issues

If you encounter bugs or have feature requests:

1. **Search existing issues** to avoid duplicates
2. **Create a new issue** with a descriptive title
3. **Provide details**:
   - Clear description of the issue/feature
   - Steps to reproduce (for bugs)
   - Expected vs actual behavior
   - R version and package versions
   - Operating system
   - Error messages and logs

### Submitting Changes

1. **Fork the repository**
   ```bash
   git clone https://github.com/SolvingLab/OmixBench.git
   cd OmixBench
   ```

2. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   # Or for bug fixes:
   # git checkout -b fix/bug-description
   # Or for documentation:
   # git checkout -b docs/documentation-improvement
   ```

   **Branch Naming Conventions:**
   - `feature/`: New features or enhancements
   - `fix/`: Bug fixes
   - `docs/`: Documentation improvements
   - `refactor/`: Code refactoring
   - `test/`: Adding or updating tests
   - `chore/`: Maintenance tasks

3. **Make your changes**
   - Follow R coding style conventions
   - Add documentation for new functions
   - Include examples in documentation
   - Add tests if applicable

4. **Test your changes**
   ```r
   # Test package installation
   devtools::install_local("llmhelper")
   devtools::install_local("OmixBenchR")
   devtools::install_local("llmflow")
   
   # Run checks (this checks for errors, warnings, and notes)
   devtools::check("llmhelper")
   devtools::check("OmixBenchR")
   devtools::check("llmflow")
   
   # Run tests if available
   devtools::test("llmhelper")
   devtools::test("OmixBenchR")
   devtools::test("llmflow")
   
   # Build documentation
   devtools::document("llmhelper")
   devtools::document("OmixBenchR")
   devtools::document("llmflow")
   ```

5. **Commit your changes**
   ```bash
   git add .
   git commit -m "Add: brief description of changes"
   ```

6. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

7. **Create a Pull Request**
   - Provide a clear title and description
   - Reference related issues (e.g., "Fixes #123" or "Closes #456")
   - Explain the rationale for changes
   - Include examples if relevant
   - Ensure all checks pass
   - Be responsive to feedback during review

## Development Setup

### Prerequisites

- **R**: Version 4.0.0 or higher
- **RStudio**: Recommended but not required
- **Git**: For version control
- **devtools**: R package for development

```r
# Install required development packages
install.packages(c("devtools", "roxygen2", "testthat", "usethis"))
```

### Setting Up Your Development Environment

1. **Clone the repository**
   ```bash
   git clone https://github.com/SolvingLab/OmixBench.git
   cd OmixBench
   ```

2. **Install dependencies**
   ```r
   # Install dependencies for each package
   devtools::install_deps("llmhelper")
   devtools::install_deps("OmixBenchR")
   devtools::install_deps("llmflow")
   ```

3. **Set up API keys** (for testing)
   ```r
   # Create a .Renviron file in your home directory
   usethis::edit_r_environ()
   
   # Add your API keys (never commit these!)
   # OPENAI_API_KEY="your-key-here"
   # ANTHROPIC_API_KEY="your-key-here"
   # DEEPSEEK_API_KEY="your-key-here"
   ```

4. **Load packages for development**
   ```r
   # Load a package for interactive development
   devtools::load_all("llmhelper")
   ```

## Coding Standards

### R Code Style

Follow standard R conventions:

- Use `snake_case` for function names and variables
- Use `<-` for assignment (not `=`)
- Indent with 2 spaces (no tabs)
- Line length should not exceed 80 characters
- Add roxygen2 documentation for all exported functions

Example:
```r
#' Execute a bioinformatics task using LLM
#'
#' @param task A list containing task description and metadata
#' @param llm_client An LLM client object
#' @param max_iterations Maximum number of retry attempts
#' @param verbose Logical; print progress messages
#' @return A list with execution results
#' @export
#' @examples
#' \dontrun{
#' result <- Execute_Task(task, llm_client, max_iterations = 10)
#' }
Execute_Task <- function(task, llm_client, max_iterations = 10, verbose = FALSE) {
  # Function implementation
}
```

### Documentation

- All exported functions must have roxygen2 documentation
- Include `@param`, `@return`, `@export`, and `@examples`
- Provide clear, concise descriptions
- Add examples that users can run
- Update vignettes when adding major features
- Keep README.md files up to date

### Python Code (for Jupyter Notebooks)

While the core framework is in R, we provide Python examples for broader accessibility:

- Follow PEP 8 style guidelines
- Use type hints where appropriate
- Add docstrings to functions
- Keep notebooks well-documented with markdown cells
- Test notebooks before committing
- Clear output before committing (to reduce file size)

```python
# Example: Clear notebook output before committing
jupyter nbconvert --clear-output --inplace *.ipynb
```

### Documentation Standards

#### Package Documentation

- Use roxygen2 for all function documentation
- Include clear descriptions, parameters, return values, and examples
- Use markdown formatting in roxygen comments
- Link related functions with `@seealso`

#### README Files

- Each package should have a clear README
- Include installation instructions
- Provide basic usage examples
- Link to full documentation

#### Tutorials and Vignettes

- Use clear, descriptive titles
- Include code examples that run without errors
- Explain both what and why
- Consider adding diagrams for complex workflows

### Testing

When adding new features:

- Consider edge cases
- Test with different LLM providers
- Verify error handling
- Check compatibility with existing code
- Add unit tests using testthat when applicable

```r
# Example test structure
testthat::test_that("function handles edge cases", {
  expect_error(my_function(NULL), "Input cannot be NULL")
  expect_equal(my_function("test"), expected_result)
})
```

**Test Coverage Areas:**
- Input validation
- Error handling
- Edge cases (empty inputs, NULL values, etc.)
- Integration with LLM providers
- Output format validation

## Package Structure

### llmhelper/
- `R/`: R source code
- `man/`: Auto-generated documentation
- `DESCRIPTION`: Package metadata
- `NAMESPACE`: Exported functions

### OmixBenchR/
- `R/`: R source code
- `inst/`: Additional files (prompts, templates)
- `man/`: Documentation
- `DESCRIPTION`: Package metadata

### llmflow/
- `R/`: R source code
- `man/`: Documentation
- `DESCRIPTION`: Package metadata
- `AutoFlow_Tasks.R`: Example workflows

## Types of Contributions

### Bug Fixes
- Fix errors in existing code
- Improve error messages
- Handle edge cases

### New Features
- Add new LLM provider integrations
- Implement additional evaluation metrics
- Extend RAG/ReAct capabilities
- Add visualization tools

### Documentation
- Improve README files
- Add tutorials and examples
- Clarify function documentation
- Create vignettes

### Performance
- Optimize slow functions
- Reduce memory usage
- Improve API efficiency

### Testing
- Add unit tests
- Create integration tests
- Improve test coverage

## Benchmark Tasks

### Adding New Tasks

If contributing new benchmark tasks:

1. **Task Specification**
   - Clear description
   - Omics domain classification
   - Complexity level (1-3)
   - Required data files
   - Expected output format

2. **Task Metadata**
   - Domain (e.g., transcriptomics, genomics)
   - Complexity dimensions (see Task Complexity Classification.md)
   - Keywords/tags
   - Difficulty level

3. **Validation**
   - Ensure task is solvable
   - Verify data availability
   - Test with multiple LLMs
   - Document expected results

### Task Quality Criteria

- **Clarity**: Unambiguous description
- **Relevance**: Real-world bioinformatics problem
- **Feasibility**: Solvable within reasonable time/resources
- **Diversity**: Covers different aspects of omics analysis
- **Validation**: Clear success criteria

## Code Review Process

All contributions undergo review:

1. Automated checks (if available)
2. Code review by maintainers
3. Testing on different platforms
4. Documentation review
5. Merge approval

### Review Criteria

- Code quality and style
- Documentation completeness
- Test coverage
- Backward compatibility
- Performance impact

## Communication

- **Issues**: For bug reports and feature requests
- **Pull Requests**: For code contributions
- **Discussions**: For questions and ideas (use GitHub Discussions if enabled)
- **Email**: liuzaoqu@163.com for private inquiries or security concerns

### Tips for Effective Communication

- Search existing issues before creating new ones
- Use clear, descriptive titles
- Provide context and examples
- Be respectful and constructive
- Follow up on your issues and PRs

## Recognition

Contributors will be:
- Listed in package documentation
- Acknowledged in release notes
- Invited to co-author future publications (for significant contributions)

## License

By contributing, you agree that your contributions will be licensed under the GNU General Public License v3.0, the same license as the project.

## Questions?

If you have questions about contributing:
- Open a GitHub issue
- Check existing documentation
- Contact the maintainers

## Acknowledgments

Thank you for helping improve OmixBench! Your contributions make bioinformatics AI tools more accessible and reliable for the research community.

## Getting Help

If you need help with your contribution:
- Check the [documentation](README.md)
- Review existing [issues](https://github.com/SolvingLab/OmixBench/issues)
- Ask questions in your pull request or issue
- Contact the maintainers at liuzaoqu@163.com

## Additional Resources

- [R Packages Book](https://r-pkgs.org/) by Hadley Wickham
- [Git and GitHub Guide](https://docs.github.com/en/get-started)
- [Bioconductor Contribution Guidelines](https://contributions.bioconductor.org/)
- [R Style Guide](http://adv-r.had.co.nz/Style.html)

---

**Maintainer**: Zaoqu Liu (liuzaoqu@163.com)  
**Last Updated**: October 2024
