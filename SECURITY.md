# Security Policy

## Reporting Security Vulnerabilities

The OmixBench team takes the security of our software seriously. If you believe you have found a security vulnerability in OmixBench, we encourage you to let us know right away.

**Please do not report security vulnerabilities through public GitHub issues.**

Instead, please report them via email to:

**liuzaoqu@163.com**

Please include the following information in your report:

- Type of issue (e.g., buffer overflow, SQL injection, cross-site scripting, etc.)
- Full paths of source file(s) related to the manifestation of the issue
- The location of the affected source code (tag/branch/commit or direct URL)
- Any special configuration required to reproduce the issue
- Step-by-step instructions to reproduce the issue
- Proof-of-concept or exploit code (if possible)
- Impact of the issue, including how an attacker might exploit the issue

This information will help us triage your report more quickly.

## Response Timeline

- We will acknowledge receipt of your vulnerability report within 3 business days
- We will provide a detailed response within 7 business days indicating the next steps in handling your report
- We will keep you informed of the progress towards a fix and full announcement
- We may ask for additional information or guidance

## Preferred Languages

We prefer all communications to be in English or Chinese (中文).

## Disclosure Policy

- We will work with you to understand and resolve the issue quickly
- Once the issue is resolved, we will publicly disclose the vulnerability in our release notes
- We will credit you for the discovery unless you prefer to remain anonymous

## Supported Versions

We release patches for security vulnerabilities for the following versions:

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Security Best Practices for Users

### API Keys and Credentials

When using OmixBench with LLM providers:

1. **Never commit API keys** to version control
2. **Use environment variables** to store sensitive credentials
3. **Rotate API keys** regularly
4. **Use separate API keys** for development and production
5. **Monitor API usage** for unusual activity

Example of secure API key usage:

```r
# Good: Using environment variables
llm_client <- llm_provider(
  api_key = Sys.getenv("OPENAI_API_KEY")
)

# Bad: Hardcoding API keys (NEVER DO THIS)
# llm_client <- llm_provider(api_key = "sk-proj-...")
```

### Data Privacy

When working with sensitive bioinformatics data:

1. **Anonymize data** before sending to LLM providers when possible
2. **Review privacy policies** of LLM providers before use
3. **Use local models** (e.g., via Ollama) for highly sensitive data
4. **Be aware** that data sent to cloud-based LLMs may be used for model training unless you opt out
5. **Comply with institutional policies** regarding data sharing and external services

### Code Execution

When using dynamic task execution features:

1. **Review generated code** before execution in production environments
2. **Run in isolated environments** when testing with untrusted inputs
3. **Set appropriate resource limits** to prevent runaway processes
4. **Monitor execution logs** for unexpected behavior
5. **Validate outputs** before using in critical analyses

### Package Dependencies

1. **Keep packages updated** to receive security patches
2. **Review dependencies** periodically for known vulnerabilities
3. **Use `sessionInfo()`** to document package versions for reproducibility
4. **Install from trusted sources** (CRAN, Bioconductor, official GitHub repos)

## Security Updates

Security updates will be released as patch versions (e.g., 1.0.1) and will be announced via:

- GitHub Security Advisories
- Release notes
- Email to the maintainer's list (if applicable)

## Additional Resources

- [R Security Best Practices](https://cran.r-project.org/web/packages/policies.html)
- [Bioconductor Security](https://www.bioconductor.org/)
- [OWASP Top Ten](https://owasp.org/www-project-top-ten/)

## Questions?

If you have questions about this security policy, please contact:

**Zaoqu Liu**  
Email: liuzaoqu@163.com

---

**Last Updated**: October 2024
