# GitHub Publishing and Governance Guide

## Repository Identity

- Repository: `https://github.com/vamsee2k1/PRS_GWAS_SNP_PIPELINE`
- Primary branch: `main`
- Owner handle: `@vamsee2k1`

## Open-Source and Citation Files

These files are present and should be kept in sync with releases:

- `LICENSE` (Apache-2.0)
- `NOTICE`
- `CITATION.cff`
- `CONTRIBUTING.md`
- `.github/CODEOWNERS`
- `.github/pull_request_template.md`
- `SECURITY.md`

## Ownership and Review Controls

Current code-owner mapping:

```text
* @vamsee2k1
```

Recommended branch-protection settings for `main`:

- Require a pull request before merging
- Require at least 1 approval
- Dismiss stale approvals when new commits are pushed
- Require status checks to pass before merging
- Require conversation resolution before merging
- Restrict direct pushes to maintainers
- Disallow force pushes
- Disallow branch deletion

Optional hardening:

- Require linear history
- Require signed commits

## Security and Disclosure

- Security reporting policy: `SECURITY.md`
- Prefer private vulnerability reporting through GitHub Security tab

## Release Workflow

Suggested first release tag:

```bash
git tag -a v0.1.0 -m "Initial public release"
git push origin v0.1.0
```

After each release:

- Update `CITATION.cff` `version` if needed
- Confirm README references remain accurate

## Publication-Readiness Checklist

- [x] Repository pushed to GitHub
- [x] `CODEOWNERS` set to repository owner (`@vamsee2k1`)
- [x] `CITATION.cff` points to correct repository URL
- [x] README documents canonical real-data test runs
- [x] Test runs documented for `gwas_summary`, `variant_only`, and `full`
  - `FINAL_TEST.md`
  - `docs/REAL_DATASET_EXAMPLES.md`
- [ ] Branch protection enabled on `main` (GitHub UI)
- [ ] First release tag created (for example `v0.1.0`)

