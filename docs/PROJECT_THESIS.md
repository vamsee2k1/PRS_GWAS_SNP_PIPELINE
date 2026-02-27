# Project Thesis: Reproducible Multi-Mode Bioinformatics Pipeline for Variant Interpretation, GWAS Signal Prioritization, and PRS Reporting

## Abstract

This thesis presents the design, implementation, and validation of a reproducible bioinformatics pipeline that unifies three practical workflows in a single Snakemake framework: (1) full FASTQ processing for alignment and variant calling, (2) variant-only interpretation from external VCF/VCF.GZ, and (3) GWAS summary harmonization from TSV/CSV tables. The system was engineered for publication-grade traceability with preflight validation, run manifests, mode-specific configuration, and isolated software environments. Final validation was completed on real datasets on February 27, 2026, with successful execution in all three target modes. Results demonstrate stable generation of quality-control metrics, biologically interpretable association and enrichment outputs, and PRS reports with explicit locus-coverage diagnostics. Historical runtime benchmarks from February 24, 2026 are retained for context, including the prior ~37-minute full-mode benchmark run.

## Keywords

- bioinformatics workflow engineering
- Snakemake reproducibility
- GWAS summary harmonization
- variant interpretation
- polygenic risk score (PRS)
- Alzheimer disease (AD) example analysis

## 1. Introduction and Motivation

Modern genomics analysis often fragments into notebooks, one-off scripts, and dataset-specific patches. This creates four recurring problems:

- low reproducibility across machines and users
- fragile input handling for heterogeneous public datasets
- weak operational traceability for publication/reporting
- difficult end-to-end validation when switching between sequencing and external variant inputs

This project addresses those constraints with a single, mode-aware workflow that preserves methodological rigor while remaining practical for real-world data ingestion.

## 2. Research Objectives and Questions

### 2.1 Objectives

- Build a single workflow engine that supports full FASTQ and external-variant analysis paths.
- Enforce preflight validation before expensive compute.
- Standardize output artifacts (QC, association, enrichment, PRS, manifest).
- Demonstrate reproducible reruns on real datasets across all modes.

### 2.2 Research Questions

1. Can one pipeline reliably support full sequencing, external VCF interpretation, and GWAS summary harmonization without branch-specific rewrites?
2. Does the workflow preserve biologically meaningful outputs under realistic public-data imperfections?
3. Can runtime and quality behavior be clearly documented for publication readiness?

## 3. System Architecture and Design

The workflow is implemented in Snakemake (`Snakefile`) with explicit mode routing:

- `full`: FASTQ -> QC -> trimming -> alignment -> variant calling/filtering -> interpretation (+ optional PRS)
- `variant_only` with `vcf_interpretation`: external VCF normalization/filtering -> interpretation (+ optional PRS)
- `variant_only` with `gwas_summary`: TSV/CSV alias mapping + REF derivation + row validation -> interpretation

Core design controls:

- preflight gate: `preflight_resources`
- mode-specific config files in `config/`
- environment isolation via `envs/*.yaml`
- deterministic output trees by run directory
- run-level provenance in `docs/run_manifest.txt`

## 4. Data, Materials, and Canonical Validation Runs

## 4.1 Canonical Real-Data Runs (Feb 27, 2026)

| Mode | Config | Input | Output Directory |
| --- | --- | --- | --- |
| `variant_only` + `gwas_summary` | `config/final_tests/final_test_gwas_advp.yaml` | `INPUT_TEST_FILES/advp.variant.records.hg38.tsv` | `results_final_test_gwas_advp_20260227` |
| `variant_only` + `vcf_interpretation` | `config/final_tests/final_test_ad_high_prs.yaml` | `data/variants/alzheimers/final_tests/1kg_chr19_high_prs_top40.vcf.gz` | `results_final_test_ad_high_prs_20260227` |
| `full` (DNA short-read) | `config/final_tests/final_test_full_giab.yaml` | `data/giab_hg002/HG002_sub_R1.fastq.gz`, `data/giab_hg002/HG002_sub_R2.fastq.gz` | `results_final_test_full_giab_20260227` |

## 4.2 Historical Benchmark Runs (Feb 24, 2026)

These are retained for runtime context:

- `results_advp_gwas_20260224_retest` (`22 sec`)
- `results_alzheimers_prs_1kg_chr19_20260224_retest` (`8 min 58 sec`)
- `results_giab_hg002_full_20260224_practical_depth_hg002` (`37 min 10 sec` wall span)

## 5. Methods

## 5.1 Preflight Validation

Every mode enters preflight checks for:

- required file availability
- build compatibility logic
- mode-specific parameter consistency
- format readiness for external variant inputs

Primary artifact:

- `results.../docs/preflight_checks.tsv`

## 5.2 Full FASTQ Methodology

Processing stages:

1. Raw read QC (`fastqc`)
2. Read trimming/filtering (`fastp`)
3. Alignment (`bwa-mem2`, `bwa`, or configured fallback)
4. BAM indexing and duplicate metrics
5. Alignment QC (`samtools flagstat`, `samtools stats`)
6. Depth profiling (`mosdepth`, optional `samtools depth`)
7. Variant calling/filtering (`bcftools`)
8. Variant interpretation, enrichment, optional PRS

## 5.3 Variant-Only VCF Methodology

Processing stages:

1. External VCF normalization
2. Filter canonicalization
3. Final VCF derivation
4. Variant QC metrics and plots
5. Variant-gene mapping and pathway enrichment
6. Optional PRS scoring and QC

## 5.4 GWAS Summary Methodology

Processing stages:

1. Header alias mapping to required fields
2. Reference-based `REF` derivation
3. Row validation and invalid-row reporting
4. Table-to-VCF normalization for downstream consistency
5. Association visuals (Manhattan/QQ/volcano), enrichment

Diagnostic artifact:

- `variants/external.input.field_report.tsv`

## 5.5 PRS Methodology

PRS branch behavior:

- model file ingestion from configured weights
- match by `rsID` and position/allele harmonization
- optional strand handling policy
- per-sample score table + percentile/z-score outputs
- explicit model-coverage diagnostics (`prs_qc.tsv`)

Interpretation rule:

- low matched-model-variant rate indicates limited interpretability of absolute PRS values.

## 6. Results

## 6.1 Runtime Results with Timestamps and Historical Context

Final run timestamps were read from Snakemake logs (`.snakemake/log/*.snakemake.log`).

| Mode | Final Run Timestamp Window (Feb 27, 2026) | Final Runtime | Historical Runtime (Feb 24, 2026) | Historical Reference |
| --- | --- | --- | --- | --- |
| `variant_only` + `gwas_summary` | `20:50:16 -> 20:50:37` | `21 sec` | `22 sec` | `results_advp_gwas_20260224_retest` |
| `variant_only` + `vcf_interpretation` | `20:48:53 -> 20:50:02` | `1 min 9 sec` | `8 min 58 sec` | `results_alzheimers_prs_1kg_chr19_20260224_retest` |
| `full` (DNA short-read) | `20:50:49 -> 20:57:13` | `6 min 24 sec` | `37 min 10 sec` | `results_giab_hg002_full_20260224_practical_depth_hg002` |

Interpretation note:

- Historical and final runtimes are informative but not strict like-for-like performance claims because run conditions, data slices, and settings differ.

## 6.2 Quality Metrics by Mode

| Mode | Preflight | Core Variant Metrics | Additional Mode-Specific Quality |
| --- | --- | --- | --- |
| `gwas_summary` | `0 ERROR`, `2 WARNING` | `total_variants=1,981`, `Ts/Tv=2.22114` | input rows `6,346` -> retained `1,981`; AD-associated hits in `APOE/TOMM40` region |
| `vcf_interpretation` (AD high-PRS-focused) | `0 ERROR`, `0 WARNING` | `total_variants=1,625,698`, `sample_count=40`, `Ts/Tv=2.35086` | PRS model variants `83`, matched `3` (`3.614%`), build match `True` |
| `full` (GIAB HG002) | `0 ERROR`, `2 WARNING` | `total_variants=13,520`, `SNP=12,556`, `INDEL=905`, `Ts/Tv=1.53657` | read retention `94.167%`; mean depth `0.089x`; PRS matched loci `0/83` |

Preflight warning context:

- GWAS run: missing raw `REF` column and intentionally disabled PRS weights in summary mode.
- Full run: missing `bwa-mem2` index at path; automatic fallback to classic `bwa`.

## 6.3 Figure Portfolio and Interpretation

### 6.3.1 GWAS Summary Mode Figures

Figure 1. Manhattan plot from ADVP-style GWAS summary harmonization:

![Figure 1: GWAS Manhattan](assets/readme/gwas/variant_manhattan.png)

Figure 2. Pathway enrichment dotplot from mapped variant genes:

![Figure 2: GWAS Enrichment Dotplot](assets/readme/gwas/variant_enrichment_dotplot.png)

Figure 3. Functional class distribution:

![Figure 3: GWAS Functional Classes](assets/readme/gwas/variant_functional_classes.png)

Interpretation:

- strong chr19 peak structure is consistent with known AD loci (`APOE/TOMM40` neighborhood)
- enrichment terms show biologically plausible AD-related lipid/cholesterol pathways

### 6.3.2 Variant-Only + PRS Figures

Figure 4. Variant class counts in high-PRS-focused chr19 subset:

![Figure 4: Variant Type Counts](assets/readme/variant_only/variant_type_counts.png)

Figure 5. Variant-gene enrichment:

![Figure 5: Variant Enrichment Dotplot](assets/readme/variant_only/variant_enrichment_dotplot.png)

Figure 6. PRS distribution:

![Figure 6: PRS Distribution](assets/readme/variant_only/prs_distribution.png)

Interpretation:

- cohort-level PRS spread is visible
- model-locus overlap is limited (`3/83`), so this run is strong for workflow validation but limited for absolute risk interpretation

### 6.3.3 Full FASTQ Mode Figures

Figure 7. QC before/after trimming:

![Figure 7: QC Before/After](assets/readme/full/qc_before_after.png)

Figure 8. Variant class counts:

![Figure 8: Full-Mode Variant Type Counts](assets/readme/full/variant_type_counts.png)

Figure 9. Variant-gene enrichment:

![Figure 9: Full-Mode Enrichment Dotplot](assets/readme/full/variant_enrichment_dotplot.png)

Interpretation:

- post-trim quality improvements are explicit (Q30 increase and high read retention)
- output quality profile is appropriate for technical benchmarking on HG002 subset data

## 7. Biological Significance

This project is biologically relevant in two distinct ways:

- hypothesis-prioritization workflows:
  - identifies high-confidence associated loci from GWAS-style inputs
  - maps loci to genes/pathways for biological interpretation
- technical validation workflows:
  - benchmarks complete sequencing-to-variant pipelines on known reference samples
  - verifies that QC, calling, and annotation branches produce coherent outputs

AD-specific notes:

- GWAS mode repeatedly surfaces known AD loci around `APOE/TOMM40`.
- PRS mode behaves correctly with transparent coverage diagnostics; low model overlap is surfaced rather than hidden.

## 8. Reproducibility and Publication Artifacts

Mandatory artifacts per run:

- `docs/preflight_checks.tsv`
- `docs/run_manifest.txt`
- mode-specific core outputs in `variants/`, `prs/`, `qc/`, and `depth/`

Publication-ready summary files:

- `FINAL_TEST.md`
- `docs/REAL_DATASET_EXAMPLES.md`
- `README.md` final validation table with historical runtime context

## 9. Limitations and Threats to Validity

- Runtime comparisons across runs can be confounded by dataset scale, caching, and configuration differences.
- PRS interpretation depends on locus overlap between model weights and callable variants.
- GWAS summary sources may require row-level filtering due to heterogeneity in raw public tables.
- Benchmark-subset depth values should not be over-interpreted as production whole-genome coverage.

## 10. Conclusions and Future Work

The project successfully establishes a reproducible, multi-mode workflow that bridges sequencing pipelines and external variant interpretation. Final real-data reruns confirm operational readiness across all target modes with consistent quality outputs and traceable provenance. The framework is suitable for publication-oriented reporting, benchmarking, and method extension.

Planned next steps:

- add optional HPC profile templates (e.g., SLURM)
- add stricter production policy for GWAS row-drop thresholds
- add containerized execution profiles (Docker/Singularity)
- add compact CI fixtures for automated regression on every commit

## 11. References

## 11.1 Repository-Internal References

1. `Snakefile` (workflow orchestration)
2. `workflow/scripts/` (analysis and plotting scripts)
3. `FINAL_TEST.md` (final rerun metrics and pass/fail summary)
4. `docs/REAL_DATASET_EXAMPLES.md` (canonical run interpretation)
5. `config/final_tests/*.yaml` (publish-ready config snapshots)

## 11.2 Data and Tool References

1. Snakemake documentation: https://snakemake.readthedocs.io/
2. fastp: https://github.com/OpenGene/fastp
3. FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
4. MultiQC: https://multiqc.info/
5. samtools: http://www.htslib.org/doc/samtools.html
6. bcftools: http://www.htslib.org/doc/bcftools.html
7. bedtools: https://bedtools.readthedocs.io/
8. mosdepth: https://github.com/brentp/mosdepth
9. 1000 Genomes / IGSR FTP: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
10. GIAB (Genome in a Bottle): https://www.nist.gov/programs-projects/genome-bottle
11. GWAS Catalog summary statistics: https://www.ebi.ac.uk/gwas/downloads/summary-statistics
12. PGS Catalog downloads: https://www.pgscatalog.org/downloads/
13. MSigDB collections: https://www.gsea-msigdb.org/gsea/msigdb
14. Bellenguez et al. AD PRS model reference (as captured in `prs_qc.tsv`): doi:10.1038/s41588-022-01024-z

