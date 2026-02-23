# Contributing

This repository is open source, but direct modification access should be limited to maintainers.

## How Contributions Should Happen

- Open an issue before large feature changes.
- Use pull requests; do not commit directly to the default branch.
- Include a reproducible test command and expected outputs.
- For pipeline behavior changes, document:
  - affected mode(s): `full`, `variant_only`, `gwas_summary`, `vcf_interpretation`
  - config keys changed
  - output schema changes (if any)
- For new resources/datasets, document provenance and license/terms.

## Coding Expectations

- Prefer small, reviewable PRs.
- Preserve backward compatibility for config keys when possible.
- Add preflight validation for new required inputs/resources.
- Add/update docs and example configs when adding new modes/branches.

## Reproducibility Requirements

Every PR that changes results should include:

- command(s) used
- config file used
- reference build (`GRCh38` / `GRCh37`)
- input dataset description (or synthetic fixture)

## Citation and Reuse

If you publish results generated with this pipeline, please cite the software repository release and describe the exact version/configuration used. See `CITATION.cff`.

## Maintainer Controls (GitHub)

GitHub settings (configured in the web UI) should enforce:

- protected default branch
- required pull request reviews
- status checks before merge
- no force pushes to default branch
- write access only for maintainers

