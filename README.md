# Reproducible Bioinformatics Pipeline (DNA/RNA + PRS)

## Biological Significance

This pipeline is designed to connect **raw sequencing and external variant data to biologically interpretable outputs** that are commonly used in genomics and translational research:

- variant calling and functional annotation for DNA sequencing data
- GWAS summary interpretation (Manhattan/QQ/volcano plots)
- variant-to-gene mapping and pathway enrichment
- optional PRS scoring from genotype-bearing VCFs

In practice, this supports biologically meaningful workflows such as:

- prioritizing disease-associated loci (for example, AD GWAS signals in `APOE/TOMM40`)
- identifying pathways enriched in variant-mapped genes (for example, lipid/cholesterol transport pathways in Alzheimer-related GWAS summaries)
- comparing preprocessing quality changes before/after trimming in FASTQ pipelines
- harmonizing external cohort VCFs for downstream interpretation and PRS compatibility checks

This project provides a **fundamentals-first, reproducible Snakemake pipeline** that supports:

1. Pre-processing with FastQC
2. Quality control before/after filtering with summary plots
3. Filtering (read-level and variant-level)
4. Alignment to reference genome (GRCh38 recommended, hg19 supported)
5. Transcript and SNP identification + functional variant annotation
6. Differential gene expression (RNA mode)
7. Enrichment analysis (RNA mode + variant-gene mode)
8. PRS scoring

All outputs are written to `results/`, and detailed documentation is in `docs/`.

## Observed Runtime (Example Test Runs)

The following runtimes were observed on a local laptop test environment and are intended as **practical expectations**, not guarantees. They vary with CPU/RAM, storage, dataset size, reference indexes, and whether conda environments are already built.

These times reflect the **main workflow runs** (not one-time dataset downloads or one-time reference index creation).

| Mode | Example Input | Observed Time | Notes |
| --- | --- | --- | --- |
| `gwas_summary` (`variant_only`) | ADVP Alzheimer GWAS TSV (`advp.variant.records.hg38.tsv`) | ~18 sec | Completed successfully; summary-stat interpretation/enrichment path |
| `variant_only` (`vcf_interpretation`) | 1000G Phase 3 chr19 VCF | ~10 min 42 sec | Completed successfully; includes PRS branch (0 matched loci in this test) |
| `full` (DNA smoke test) | Local FASTQ subset (`NIST7035` subset) | ~1 min 44 sec to complete variant interpretation outputs | Core QC/alignment/calling/annotation outputs completed; depth branch later became bottleneck |
| `full` (same run, depth enabled) | Same FASTQ subset | `depth_per_sample` finished by ~2 min 49 sec; `depth_summary` failed at ~6 min 55 sec | Generated a ~39 GB depth file (`samtools depth -a`), then `depth_summary` exited (code 137) |

Recommended local testing settings for `full` mode:

- `run.enable_depth: false` for fast functional validation
- or `depth.emit_all_positions: false` to avoid huge `samtools depth -a` outputs

Additional runtime guidance for local/HPC/container setups is in `docs/REPRODUCIBLE_EXECUTION.md`.
Alzheimer PRS model details and calculation notes are in `docs/ALZHEIMERS_PRS.md`.
Alzheimer dataset files and caveats are in `docs/ALZHEIMERS_DATASETS.md`.
Statistical threshold rationale is documented in `docs/STATISTICAL_THRESHOLDS.md`.
Project architecture/thesis-style design rationale is in `docs/PROJECT_THESIS.md`.
GitHub publishing, branch protection, and citation setup is in `docs/GITHUB_PUBLISHING.md`.
Recent regression-test and example-file validation evidence is in `docs/validation/2026-02-23_regression_checks.md`.

## Pipeline Architecture (High-Level)

```mermaid
flowchart TD
    A([Start]) --> B["Load config and select mode"]
    B --> C["Preflight validation"]
    C --> D{"Mode"}

    D -->|Full mode| E["FASTQ workflow"]
    E --> E1["QC and trimming"]
    E1 --> E2["Alignment and BAM processing"]
    E2 --> E3["Variant calling and filtering"]
    E3 --> J["Final filtered VCF"]

    D -->|Variant-only mode| F["External variants input"]
    F --> F1{"Input type"}
    F1 -->|VCF/VCF.GZ| F2["Normalize and filter variants"]
    F1 -->|CSV/TSV| F3["Convert table to VCF and normalize"]
    F2 --> J
    F3 --> J

    J --> K["Variant annotation and gene mapping"]
    K --> L["Variant QC plots and interpretation visuals"]
    K --> M["Pathway enrichment"]

    J --> N{"PRS weights configured?"}
    N -->|Yes| O["PRS scoring and PRS QC"]
    N -->|No| P["Skip PRS"]

    E2 --> Q{"RNA mode?"}
    Q -->|Yes| R["RNA quantification, DE analysis, enrichment"]
    Q -->|No| S["Skip RNA branch"]

    L --> T["Run manifest and final outputs"]
    M --> T
    O --> T
    R --> T
    C --> T
    T --> U([End])
```

Detailed split flowcharts (preflight/mode routing, full mode, and variant-only + interpretation + PRS) are documented in `docs/PIPELINE_DETAILS.md`.

## Quick Start

1. Prepare inputs using `docs/USER_INPUT_FORMATS.md`.
2. Update `config/config.yaml`.
3. Run:

```bash
snakemake --use-conda --cores 8
```

If you want a terminal loading spinner between rule logs:

```bash
./run_pipeline.sh --use-conda --cores 8
```

If Snakemake is not installed, create it from `envs/workflow.yaml` first.

For dry-run:

```bash
snakemake -n
```

For an explicit preflight validation-only check (resources, formats, build compatibility):

```bash
snakemake --use-conda --cores 1 --until preflight_resources
```

Mode-specific example configs (ready to copy and edit):

- `config/examples.mode_gwas_summary.yaml`
- `config/examples.mode_variant_only_vcf.yaml`
- `config/examples.mode_full_dna_fastq.yaml`

Variant-only mode (VCF/CSV/TSV to filtered variants + interpretation):

1. Set `run.mode: variant_only`
2. Set `run.variant_data_mode` to `vcf_interpretation` or `gwas_summary`
3. Set `paths.variants_input` to `.vcf`, `.vcf.gz`, `.csv`, or `.tsv`
4. Run `snakemake --use-conda --cores 8`

Direct shorthand modes are also supported:
- `run.mode: vcf_interpretation`
- `run.mode: gwas_summary`

Alzheimer PRS production config:

```bash
snakemake --use-conda --cores 8 --configfile config/config.alzheimers_prs.yaml
```

Real Alzheimer disease-annotated ClinVar variant run (no synthetic sample genotypes):

```bash
snakemake --use-conda --cores 8 --configfile config/config.alzheimers_clinvar.yaml
```

## Main Configuration

- `run.mode`: `full`, `variant_only`, `vcf_interpretation`, or `gwas_summary`
- `run.variant_data_mode`: `auto`, `vcf_interpretation`, `gwas_summary`
- `run.assay`: `dna` or `rna`
- `run.read_type`: `short` or `long`
- `run.enable_depth`: enable/disable depth outputs for quick tests vs production runs
- `reference.fasta`: reference genome FASTA
- `reference.gtf`: gene annotation GTF
- `reference.star_index`: required for RNA short-read mode
- `alignment.dna_short_aligner`: `auto`, `bwa_mem2`, `bwa`, `minimap2_sr`
- `depth.emit_all_positions`: use `samtools depth -a` (large output) vs covered positions only
- `paths.variants_input`: optional external `VCF/VCF.GZ/CSV/TSV`
- `paths.prs_weights`: PRS model file (PGS Catalog or simple TSV)
- `paths.gene_sets`: enrichment GMT file
- `annotation.method`: `overlap` or `snpeff`
- `annotation.snpeff_database`: required when `annotation.method=snpeff`
- `prs.weights_format`: `auto`, `pgs_catalog`, or `simple_tsv`
- `prs.allow_ambiguous_strand`: include/exclude ambiguous A/T and C/G SNPs
- `validation.enforce_build_match`: fail preflight when detected VCF build mismatches expected reference build
- `validation.fail_on_warning`: treat preflight warnings as hard failures
- `thresholds.min_variant_depth`: optional INFO/DP filter for variant records
- `thresholds.require_filter_pass`: require `FILTER=PASS` for external VCF filtering
- `thresholds.de_padj`: DE adjusted p-value cutoff (default `0.05`)
- `thresholds.de_abs_log2fc`: DE effect-size cutoff (default `1.0`)
- `thresholds.enrichment_qvalue`: enrichment q-value/FDR cutoff (default `0.05`)
- `thresholds.variant_assoc_pvalue`: variant association cutoff (default `5e-8`)
- `thresholds.variant_assoc_abs_effect`: variant association effect cutoff (default `0.0`)

## Reproducibility Features

- Workflow engine: Snakemake
- Isolated environments: per-step conda envs in `envs/`
- Deterministic output tree under `results/`
- Early preflight validation report: `results/docs/preflight_checks.tsv`
- Run metadata: `results/docs/run_manifest.txt`

## Output Structure

- `results/qc/`: FastQC + MultiQC + before/after QC summary
- `results/qc/alignment/`: per-sample alignment QC (`samtools flagstat`, `samtools stats`)
- `results/trimmed/`: filtered reads and fastp reports
- `results/alignment/`: sorted BAM and BAM index
- `results/depth/`: depth files, summary table, depth plot, mosdepth summaries
- `results/variants/`: called/filtered/final VCFs
- `results/variants/snpeff_summary.html`: optional `snpEff` consequence annotation report (`annotation.method=snpeff`)
- `results/variants/variant_qc_metrics.tsv`: SNP/indel counts, Ts/Tv, QUAL/DP/GQ summary metrics
- `results/variants/sample_missingness.tsv`: per-sample genotype missingness
- `results/variants/variant_type_counts.png`: variant class counts
- `results/variants/variant_qual_distribution.png`: QUAL histogram
- `results/variants/variant_dp_distribution.png`: DP histogram
- `results/variants/variant_gq_distribution.png`: GQ histogram
- `results/variants/variant_het_allele_balance.png`: heterozygous allele-balance (AD)
- `results/variants/variant_missingness_by_sample.png`: sample missingness plot
- `results/variants/annotated_variants.tsv`: per-variant functional class + mapped genes
- `results/variants/variant_genes.tsv`: genes hit by variants (for biology/enrichment)
- `results/variants/variant_association_hits.tsv`: significant variants by threshold
- `results/variants/variant_functional_classes.png`: functional-class distribution
- `results/variants/variant_chromosome_counts.png`: chromosome-level variant counts
- `results/variants/variant_volcano.png`: association volcano (or proxy prioritization from QUAL/AF/clinical-significance when GWAS p/effect are absent)
- `results/variants/variant_manhattan.png`: Manhattan plot (reported p-values or proxy ranking)
- `results/variants/variant_qq.png`: QQ plot for GWAS-style reported p-values
- `results/variants/variant_heatmap_top.png`: top-variant genotype/dosage heatmap; falls back to numeric variant-feature heatmap when sample genotype columns are unavailable
- `results/variants/variant_enrichment.tsv` + `variant_enrichment_dotplot.png`: pathway enrichment from variant-mapped genes
- `results/transcripts/`: StringTie transcripts + gene counts (RNA)
- `results/dge/`: DESeq2 result tables + volcano + heatmap (RNA)
- `results/enrichment/`: enrichment table + dotplot (RNA)
- `results/prs/`: PRS score output + PRS summary table + PRS QC + PRS distribution plot
- `results/docs/`: run manifest

## Example Results (From Real Test Runs)

The plots below are generated by this pipeline from the test runs executed during validation of the three modes. These are included as examples of output structure and interpretation workflow.

### 1) GWAS Summary Mode (`ADVP` Alzheimer GWAS TSV, GRCh38)

Example run type:
- `variant_only` + `run.variant_data_mode: gwas_summary`
- Input: `advp.variant.records.hg38.tsv` (ADVP/NIAGADS-style GWAS summary table)

Key biological signal observed in this test run:
- Strong association peak around chromosome 19 (`APOE` / `TOMM40` region), consistent with known Alzheimer disease GWAS signals.
- Enrichment highlights cholesterol/lipid transport and related pathways (e.g., `NR1H3/NR1H2`, `ABC transporters`), which is biologically plausible in AD genetics.

Example top hits seen in generated `variant_association_hits.tsv`:
- `rs2075650` (`TOMM40`)
- `rs769449` (`APOE`)

![GWAS Manhattan Plot (ADVP example)](docs/assets/readme/gwas/variant_manhattan.png)

![GWAS Pathway Enrichment Dotplot (ADVP example)](docs/assets/readme/gwas/variant_enrichment_dotplot.png)

![GWAS Functional Class Distribution (ADVP example)](docs/assets/readme/gwas/variant_functional_classes.png)

Caveat:
- ADVP is a summary-statistics dataset, not a sample genotype VCF, so it supports GWAS interpretation and enrichment but not individual-level PRS scoring.

### 2) Variant-Only Mode (External VCF, 1000 Genomes chr19 example)

Example run type:
- `variant_only` + `run.variant_data_mode: vcf_interpretation`
- Input: `ALL.chr19.phase3...vcf.gz` (1000 Genomes chr19)

What this demonstrates:
- External VCF normalization/filtering
- Variant class/QC summaries from genotype-bearing VCF input
- Gene/pathway enrichment from mapped variant genes
- PRS harmonization branch behavior

Interpretation notes from this test:
- The enrichment output on chr19 is dominated by chromosome-content effects (gene density and immune/lipid-related gene clusters on chr19), so this is a mechanics/annotation test rather than disease-specific inference.
- PRS distribution is flat in this run because the configured AD PRS model (`PGS002280`) had `0` matched loci in this chr19 test slice (`prs_qc.tsv` confirms `matched_model_variants = 0`).

![Variant Type Counts (1000G chr19 example)](docs/assets/readme/variant_only/variant_type_counts.png)

![Variant Enrichment Dotplot (1000G chr19 example)](docs/assets/readme/variant_only/variant_enrichment_dotplot.png)

![PRS Distribution (1000G chr19 example)](docs/assets/readme/variant_only/prs_distribution.png)

### 3) Full Mode (FASTQ -> Alignment -> Variant Calling; local smoke test)

Example run type:
- `full` mode with DNA input
- Local subset FASTQ used for fast end-to-end path validation

What this demonstrates:
- Read QC and trimming
- Alignment
- Variant calling/filtering
- Variant interpretation outputs

Observed test-run QC outcome (example):
- Read retention after filtering: `94.574%`
- Q30 fraction improved from `0.923595` to `0.949523`

Observed test-run variant QC outcome (example):
- `total_variants = 1962`
- `snp_count = 1870`
- `ts/tv = 1.46377`

These values are from a small local test subset and are useful for pipeline validation, not clinical or cohort-level interpretation.

![QC Before/After Summary (Full mode example)](docs/assets/readme/full/qc_before_after.png)

![Variant Type Counts (Full mode example)](docs/assets/readme/full/variant_type_counts.png)

![Variant Enrichment Dotplot (Full mode example)](docs/assets/readme/full/variant_enrichment_dotplot.png)

## Notes

- **GRCh38 is recommended** as the current standard reference.
- Use hg19/GRCh37 only when you must match legacy cohorts.
- Differential expression and enrichment are RNA-focused analyses.
- Default PRS model is Alzheimer disease PGS Catalog score `PGS002280` (GRCh38).
- If no Alzheimer model loci are matched in the VCF, PRS output reports score `0` with an explicit non-interpretable coverage note.
- This repository is open source (Apache-2.0). If used in research/benchmarking/publications, cite the repository release and report the exact version/config/reference resources used (see `CITATION.cff`).
