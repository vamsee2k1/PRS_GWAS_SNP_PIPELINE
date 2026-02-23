diff --git a//Users/vamseea/Documents/bioinfo-pipeline/claude_test_pipeline/github/SECURITY.md b//Users/vamseea/Documents/bioinfo-pipeline/claude_test_pipeline/github/SECURITY.md
new file mode 100644
--- /dev/null
+++ b//Users/vamseea/Documents/bioinfo-pipeline/claude_test_pipeline/github/SECURITY.md
@@ -0,0 +1,99 @@
+# Security Policy
+
+## Is This Required?
+
+No, a `SECURITY.md` file is not strictly required to publish an open-source repository.
+
+However, it is strongly recommended for public projects because it tells users:
+
+- how to report vulnerabilities safely (without posting them publicly)
+- which versions/branches you actively support with security fixes
+- what response timeline to expect
+
+GitHub also uses this file in the repository Security tab and for private vulnerability reporting workflows.
+
+## Supported Versions
+
+This project is currently in early release stage. Security fixes are provided on a best-effort basis for:
+
+| Version / Branch | Supported | Notes |
+| --- | --- | --- |
+| `main` | Yes | Active development branch |
+| Latest tagged release | Yes | Recommended for reproducible use |
+| Older tags/releases | No | Upgrade to latest release |
+| Unreleased forks / modified copies | No | Maintainers cannot guarantee patch support |
+
+## Reporting a Vulnerability
+
+### Preferred (Private) Reporting
+
+If you discover a potential security issue, **do not open a public GitHub issue with exploit details**.
+
+Please use one of the following:
+
+1. **GitHub Private Vulnerability Reporting** (preferred)
+   - Go to the repository Security tab and use “Report a vulnerability” (if enabled).
+
+2. **If private reporting is not enabled yet**
+   - Open a normal issue with minimal detail (e.g., “Please contact me regarding a security concern”)
+   - Do **not** include exploit steps, secrets, tokens, or sensitive file paths in the public issue
+
+### What to Include in a Report
+
+Please include:
+
+- affected file(s) / script(s)
+- version / commit hash (or branch)
+- reproduction steps (private report only)
+- expected vs actual behavior
+- impact assessment (e.g., command injection, path traversal, secret leakage)
+- suggested fix (optional)
+
+## Response Expectations
+
+Best-effort maintainer targets:
+
+- Acknowledgement: within **5 business days**
+- Initial triage: within **14 business days**
+- Status update cadence: at least every **14 business days** while actively triaging
+
+If accepted:
+
+- a fix will be developed for `main`
+- maintainers may prepare a patch release/tag (if applicable)
+- coordinated disclosure timing will be discussed before public details are posted
+
+If declined:
+
+- maintainers will provide a brief explanation (e.g., out of scope, not reproducible, not a security issue)
+
+## Scope (Examples)
+
+Security-relevant examples for this repository may include:
+
+- command injection via unsanitized shell arguments in workflow rules/scripts
+- path traversal or unsafe file writes from user-supplied input paths
+- accidental exposure of secrets/tokens in logs or config examples
+- unsafe temporary-file handling
+- dependency vulnerabilities that result in code execution risk in the pipeline context
+
+## Out of Scope (Examples)
+
+The following are generally not treated as security vulnerabilities in this repo:
+
+- incorrect biological interpretation due to mismatched reference builds
+- data quality problems in third-party public datasets
+- expected failures caused by missing local tools/resources/indexes
+- performance issues (unless they enable a security impact such as denial-of-service in a hosted deployment)
+
+## Security Best Practices for Users
+
+- Review and validate external input files before running
+- Use trusted reference resources and verify checksums where possible
+- Keep conda environments and dependencies updated
+- Avoid running the pipeline with elevated privileges
+- Do not commit secrets, API keys, or patient-identifiable data to the repository
+
+## Disclosure and Credit
+
+Responsible disclosure is appreciated. If you would like public credit for a valid report, mention that in your private report and maintainers will include credit where appropriate.
