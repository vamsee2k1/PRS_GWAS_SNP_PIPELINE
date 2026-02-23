# GitHub Publishing and Governance Guide

## What You Asked For

You want this project to be:

- open source
- reusable with proper attribution/citation
- hard to modify directly on your main branch without review
- publication-ready

This is achievable with a combination of repository files (already added) and GitHub web settings (you must enable them in GitHub UI).

## Files Added for Open-Source + Citation Readiness

- `LICENSE` (Apache-2.0)
- `NOTICE`
- `CITATION.cff`
- `CONTRIBUTING.md`
- `.github/CODEOWNERS`
- `.github/pull_request_template.md`

These make the repository clearly reusable while documenting expectations for attribution and review.

## What You Still Need to Customize Before Publishing

Update these placeholders:

- `.github/CODEOWNERS`
  - replace `@your-github-username` with your real GitHub handle
- `CITATION.cff`
  - `repository-code`
  - author name formatting (if needed)
  - version and release/tag once you publish

## Git Initialization (Local)

This folder is currently not a git repository. Run:

```bash
cd /Users/vamseea/Documents/bioinfo-pipeline/claude_test_pipeline
git init
git add .
git commit -m "Initial bioinformatics pipeline with multi-mode support and publication metadata"
```

## Create GitHub Repository and Push

### Option A: GitHub CLI (`gh`) (recommended)

```bash
gh auth login
gh repo create <your-repo-name> --private --source=. --remote=origin --push
```

If you want it public:

```bash
gh repo create <your-repo-name> --public --source=. --remote=origin --push
```

### Option B: GitHub Web UI + git remote

1. Create a new repo in GitHub web UI.
2. Then run:

```bash
git remote add origin git@github.com:<your-username>/<your-repo-name>.git
git branch -M main
git push -u origin main
```

## Prevent “Anyone Can Edit Main” (GitHub Settings)

Go to: `Settings -> Branches -> Add branch protection rule`

Recommended rule for `main`:

- Require a pull request before merging
- Require approvals (at least 1)
- Dismiss stale approvals when new commits are pushed
- Require status checks to pass before merging
- Require conversation resolution before merging
- Restrict who can push to matching branches (maintainers only)
- Do not allow force pushes
- Do not allow deletions

Optional:
- Require linear history
- Require signed commits (if your workflow uses them)

## Enforce Reviewer Ownership (CODEOWNERS)

`.github/CODEOWNERS` ensures you (or listed maintainers) are automatically requested for review. Example:

```text
* @vamseea
Snakefile @vamseea
workflow/scripts/ @vamseea
envs/ @vamseea
config/ @vamseea
```

## Open Source + Attribution Language (What to Say in README)

Recommended wording:

- This project is open source under Apache-2.0.
- Reuse is permitted under the license terms.
- If used in research/benchmarking/publications, cite the repository release and report version/configuration/reference build.

The repo already includes `NOTICE` and `CITATION.cff` to support this.

## Release and Publication Readiness Checklist

- [ ] `CITATION.cff` placeholders updated
- [ ] `CODEOWNERS` updated with real GitHub usernames
- [ ] repo pushed to GitHub
- [ ] branch protection enabled on `main`
- [ ] first release tag created (e.g., `v0.1.0`)
- [ ] README updated with exact repository URL
- [ ] test runs documented for `gwas_summary`, `variant_only`, and `full`

## Suggested First Release Tag

```bash
git tag -a v0.1.0 -m "Initial public release"
git push origin v0.1.0
```

After tagging, update `CITATION.cff` `version` and release metadata if desired.
