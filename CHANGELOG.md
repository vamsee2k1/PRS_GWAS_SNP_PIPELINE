# Changelog

All notable changes to this project are documented in this file.

The format is based on Keep a Changelog, and this project follows Semantic Versioning.

## [0.1.0] - 2026-03-02

### Added
- Unified Snakemake pipeline supporting:
  - `full` mode (FASTQ DNA/RNA workflow branches)
  - `variant_only` mode with `vcf_interpretation` and `gwas_summary` routing
- AI add-on modules:
  - QC anomaly detection
  - variant prioritization
  - explainer report generation
- Bundled public-data-derived `test_data` pack for the 3 main modes with matching configs under `config/test_data/`.
- Documentation set for:
  - biological significance deep dive
  - reproducible execution
  - input formats
  - AI integration
  - real dataset examples
- Publication/governance files:
  - `CITATION.cff`
  - `CONTRIBUTING.md`
  - `SECURITY.md`
  - Apache-2.0 `LICENSE` + `NOTICE`

### Changed
- README expanded with:
  - clone-to-results instructions
  - generic YAML templates
  - real-run runtime/quality tables
  - PRS/AI usage guardrails
  - external resource references
- `run_pipeline.sh` default terminal output simplified via `run_with_loading.py`, with full verbose mode available using `PIPELINE_VERBOSE=1`.

### Fixed
- External GWAS TSV conversion path now supports REF derivation from FASTA when configured.
- Refined variant-only normalization/filtering behavior for ADVP-style GWAS input tables.
