# Changelog

All notable changes to the OmixBench project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- CODE_OF_CONDUCT.md for community guidelines
- SECURITY.md for security policy and vulnerability reporting
- .gitignore for better repository hygiene
- CHANGELOG.md for version tracking
- Root-level LICENSE file

### Changed
- Enhanced CITATION.cff with complete author list and preferred citation
- Improved CONTRIBUTING.md with additional sections
- Enhanced README.md with better structure and badges

### Fixed
- Corrected date in CITATION.cff (was 2025-10-28, now 2024-10-28)

## [1.0.0] - 2024-10-28

### Added
- Initial release of OmixBench framework
- Three R packages: llmhelper, OmixBenchR, and llmflow
- Static evaluation dataset (OmixTask1002) with 1,002 tasks
- Dynamic evaluation dataset (OmixQA) with 405 tasks
- Comprehensive documentation and tutorials
- Support for multiple LLM providers (OpenAI, Anthropic, DeepSeek, Qwen, Google, Ollama)
- RAG-enhanced ReAct framework for complex workflows
- Automated task execution with error recovery
- Evaluation scripts and visualization tools

### Package Versions
- llmhelper: 1.0.0
- OmixBenchR: 1.0.0
- llmflow: 2.0.0

## Notes

### Version Numbering
- **Major version** (X.0.0): Incompatible API changes
- **Minor version** (1.X.0): New functionality in a backwards compatible manner
- **Patch version** (1.0.X): Backwards compatible bug fixes

### Types of Changes
- **Added**: New features
- **Changed**: Changes in existing functionality
- **Deprecated**: Soon-to-be removed features
- **Removed**: Removed features
- **Fixed**: Bug fixes
- **Security**: Vulnerability fixes

---

For more details on each release, see the [Releases page](https://github.com/SolvingLab/OmixBench/releases).
