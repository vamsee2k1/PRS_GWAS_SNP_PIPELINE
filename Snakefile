import os
from datetime import datetime
from pathlib import Path

import pandas as pd

configfile: "config/config.yaml"

OUTDIR = config.get("output_dir", "results")
RUN = config.get("run", {})
ASSAY = RUN.get("assay", "rna").lower()
READ_TYPE = RUN.get("read_type", "short").lower()
MODE_RAW = RUN.get("mode", "full").lower()
VARIANT_DATA_MODE_CFG = RUN.get("variant_data_mode", "auto").lower()
THREADS_DEFAULT = int(RUN.get("threads", 8))
PRS_USE_EXTERNAL = bool(RUN.get("use_external_variants_for_prs", True))

SAMPLESHEET = config.get("samplesheet", "config/samples.tsv")
METADATA = config.get("metadata", "config/metadata.tsv")

REFERENCE = config.get("reference", {})
REF_FASTA = REFERENCE.get("fasta", "resources/reference/genome.fa")
REF_GTF = REFERENCE.get("gtf", "resources/annotation/genes.gtf")
STAR_INDEX = REFERENCE.get("star_index", "")
STAR_INDEX_MARKER = f"{STAR_INDEX}/.snakemake_star_index.done" if STAR_INDEX else ""

ALIGNMENT = config.get("alignment", {})
DNA_SHORT_ALIGNER = str(ALIGNMENT.get("dna_short_aligner", "auto")).strip().lower()
LONG_READ_PRESET_DNA = ALIGNMENT.get("long_read_preset_dna", "map-ont")
LONG_READ_PRESET_RNA = ALIGNMENT.get("long_read_preset_rna", "splice")
DEPTH_CFG = config.get("depth", {})
DEPTH_EMIT_ALL_POSITIONS = bool(DEPTH_CFG.get("emit_all_positions", True))
DEPTH_EMIT_PER_BASE_TSV = bool(DEPTH_CFG.get("emit_per_base_tsv", True))
RUN_ENABLE_DEPTH = bool(RUN.get("enable_depth", True))

THRESHOLDS = config.get("thresholds", {})
MIN_BASEQ = int(THRESHOLDS.get("min_base_quality", 20))
MIN_READ_LENGTH = int(THRESHOLDS.get("min_read_length", 50))
MIN_VAR_QUAL = int(THRESHOLDS.get("min_variant_qual", 30))
MIN_VAR_DP = int(THRESHOLDS.get("min_variant_depth", 0))
REQUIRE_FILTER_PASS = bool(THRESHOLDS.get("require_filter_pass", True))
DE_PADJ = float(THRESHOLDS.get("de_padj", 0.05))
DE_ABS_LOG2FC = float(THRESHOLDS.get("de_abs_log2fc", 1.0))
ENRICH_QVALUE = float(THRESHOLDS.get("enrichment_qvalue", 0.05))
VAR_ASSOC_PVALUE = float(THRESHOLDS.get("variant_assoc_pvalue", 5e-8))
VAR_ASSOC_ABS_EFFECT = float(THRESHOLDS.get("variant_assoc_abs_effect", 0.0))
VAR_HEATMAP_TOP_N = int(THRESHOLDS.get("variant_heatmap_top_n", 50))
VAR_QC_MAX_GT_SAMPLES = int(THRESHOLDS.get("variant_qc_max_samples_for_genotype_metrics", 1000))

PATHS = config.get("paths", {})
EXT_VARIANTS = str(PATHS.get("variants_input", "")).strip()
PRS_WEIGHTS = str(PATHS.get("prs_weights", "")).strip()
GENE_SETS = str(PATHS.get("gene_sets", "")).strip()

PRS_CFG = config.get("prs", {})
PRS_WEIGHTS_FORMAT = str(PRS_CFG.get("weights_format", "auto")).strip().lower()
PRS_PREFER_DS = bool(PRS_CFG.get("prefer_ds", True))
PRS_ALLOW_AMBIGUOUS_STRAND = bool(PRS_CFG.get("allow_ambiguous_strand", False))

ANNOTATION_CFG = config.get("annotation", {})
ANNOTATION_METHOD = str(ANNOTATION_CFG.get("method", "overlap")).strip().lower()
SNPEFF_DATABASE = str(ANNOTATION_CFG.get("snpeff_database", "")).strip()

VALIDATION_CFG = config.get("validation", {})
VALIDATE_ENFORCE_BUILD_MATCH = bool(VALIDATION_CFG.get("enforce_build_match", True))
VALIDATE_FAIL_ON_WARNING = bool(VALIDATION_CFG.get("fail_on_warning", False))

if MODE_RAW in {"vcf_interpretation", "gwas_summary"}:
    MODE = "variant_only"
    VARIANT_DATA_MODE = MODE_RAW
else:
    MODE = MODE_RAW
    VARIANT_DATA_MODE = VARIANT_DATA_MODE_CFG


def infer_variant_data_mode(path):
    low = str(path).lower()
    if low.endswith((".csv", ".tsv")):
        return "gwas_summary"
    return "vcf_interpretation"


if VARIANT_DATA_MODE not in {"auto", "vcf_interpretation", "gwas_summary"}:
    raise ValueError("run.variant_data_mode must be one of: auto, vcf_interpretation, gwas_summary")

if MODE not in {"full", "variant_only"}:
    raise ValueError("run.mode must be one of: full, variant_only")
if ASSAY not in {"dna", "rna"}:
    raise ValueError("run.assay must be one of: dna, rna")
if READ_TYPE not in {"short", "long"}:
    raise ValueError("run.read_type must be one of: short, long")

SAMPLES = []
SAMPLE_TO_R1 = {}
SAMPLE_TO_R2 = {}

if MODE != "variant_only":
    samples_df = pd.read_csv(SAMPLESHEET, sep="\t", dtype=str).fillna("")
    required_cols = {"sample", "fastq_1", "fastq_2"}
    if not required_cols.issubset(set(samples_df.columns)):
        missing = ", ".join(sorted(required_cols - set(samples_df.columns)))
        raise ValueError(f"samplesheet is missing required columns: {missing}")

    SAMPLES = samples_df["sample"].tolist()
    if not SAMPLES:
        raise ValueError("samplesheet has no samples; add at least one row")

    SAMPLE_TO_R1 = dict(zip(samples_df["sample"], samples_df["fastq_1"]))
    SAMPLE_TO_R2 = dict(zip(samples_df["sample"], samples_df["fastq_2"]))

if MODE == "variant_only" and not EXT_VARIANTS:
    raise ValueError("run.mode=variant_only requires paths.variants_input (VCF/VCF.GZ/CSV)")

if MODE == "variant_only" and VARIANT_DATA_MODE == "auto":
    ACTIVE_VARIANT_DATA_MODE = infer_variant_data_mode(EXT_VARIANTS)
else:
    ACTIVE_VARIANT_DATA_MODE = "vcf_interpretation" if MODE != "variant_only" else VARIANT_DATA_MODE

if ANNOTATION_METHOD not in {"overlap", "snpeff"}:
    raise ValueError("annotation.method must be one of: overlap, snpeff")
if ANNOTATION_METHOD == "snpeff" and not SNPEFF_DATABASE:
    raise ValueError("annotation.snpeff_database is required when annotation.method=snpeff")
if DNA_SHORT_ALIGNER not in {"auto", "bwa_mem2", "bwa", "minimap2_sr"}:
    raise ValueError("alignment.dna_short_aligner must be one of: auto, bwa_mem2, bwa, minimap2_sr")


def has_paired_end(sample):
    r2 = str(SAMPLE_TO_R2.get(sample, "")).strip()
    return r2 not in {"", "NA", "na", ".", "None", "none"}


def star_index_outputs(wildcards=None):
    if ASSAY == "rna" and READ_TYPE == "short":
        if not STAR_INDEX:
            raise ValueError("reference.star_index is required for RNA short-read mode")
        return [STAR_INDEX_MARKER]
    return []


def variant_bam_for_sample(sample):
    if MODE != "variant_only" and ASSAY == "dna" and READ_TYPE == "short":
        return f"{OUTDIR}/alignment/{sample}.dedup.bam"
    return f"{OUTDIR}/alignment/{sample}.sorted.bam"


def variant_bai_for_sample(sample):
    return f"{variant_bam_for_sample(sample)}.bai"


def variant_bam_outputs(wildcards=None):
    return [variant_bam_for_sample(s) for s in SAMPLES]


def variant_bai_outputs(wildcards=None):
    return [variant_bai_for_sample(s) for s in SAMPLES]


def annotation_vcf_source(wildcards=None):
    if ANNOTATION_METHOD == "snpeff":
        return f"{OUTDIR}/variants/final.snpeff.vcf.gz"
    return f"{OUTDIR}/variants/final.filtered.vcf.gz"


def annotation_vcf_index_source(wildcards=None):
    return f"{annotation_vcf_source()}.tbi"


def has_bwa_mem2_index(ref_path):
    return Path(f"{ref_path}.bwt.2bit.64").exists()


def has_bwa_index(ref_path):
    suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    return all(Path(f"{ref_path}{suf}").exists() for suf in suffixes)


def select_dna_short_aligner(ref_path):
    if DNA_SHORT_ALIGNER == "auto":
        if has_bwa_mem2_index(ref_path):
            return "bwa_mem2"
        if has_bwa_index(ref_path):
            return "bwa"
        return "bwa_mem2"
    return DNA_SHORT_ALIGNER


def core_targets():
    targets = [f"{OUTDIR}/docs/preflight_checks.tsv"]

    if MODE != "variant_only":
        targets.extend(expand(f"{OUTDIR}/qc/raw/{{sample}}.done", sample=SAMPLES))
        targets.extend(expand(f"{OUTDIR}/qc/trimmed/{{sample}}.done", sample=SAMPLES))
        targets.extend(expand(f"{OUTDIR}/alignment/{{sample}}.sorted.bam.bai", sample=SAMPLES))
        targets.extend(variant_bai_outputs())
        targets.extend(
            [
                f"{OUTDIR}/qc/qc_before_after.tsv",
                f"{OUTDIR}/qc/qc_before_after.png",
                f"{OUTDIR}/qc/multiqc_report.html",
            ]
        )
        targets.extend(expand(f"{OUTDIR}/qc/alignment/{{sample}}.flagstat.txt", sample=SAMPLES))
        targets.extend(expand(f"{OUTDIR}/qc/alignment/{{sample}}.stats.txt", sample=SAMPLES))
        if RUN_ENABLE_DEPTH:
            targets.extend(expand(f"{OUTDIR}/depth/{{sample}}.depth.tsv", sample=SAMPLES))
            targets.extend(expand(f"{OUTDIR}/depth/{{sample}}.mosdepth.summary.txt", sample=SAMPLES))
            targets.extend(
                [
                    f"{OUTDIR}/depth/depth_summary.tsv",
                    f"{OUTDIR}/depth/depth_distribution.png",
                ]
            )
        if ASSAY == "dna" and READ_TYPE == "short":
            targets.extend(expand(f"{OUTDIR}/alignment/{{sample}}.dedup.metrics.txt", sample=SAMPLES))

    if MODE != "variant_only":
        targets.append(f"{OUTDIR}/variants/called.filtered.vcf.gz")
        targets.append(f"{OUTDIR}/variants/called.filtered.vcf.gz.tbi")
    elif EXT_VARIANTS:
        targets.append(f"{OUTDIR}/variants/external.filtered.vcf.gz")
        targets.append(f"{OUTDIR}/variants/external.filtered.vcf.gz.tbi")

    targets.append(f"{OUTDIR}/variants/final.filtered.vcf.gz")
    targets.append(f"{OUTDIR}/variants/final.filtered.vcf.gz.tbi")
    if ANNOTATION_METHOD == "snpeff":
        targets.append(f"{OUTDIR}/variants/final.snpeff.vcf.gz")
        targets.append(f"{OUTDIR}/variants/final.snpeff.vcf.gz.tbi")
        targets.append(f"{OUTDIR}/variants/snpeff_summary.html")
    targets.extend(
        [
            f"{OUTDIR}/variants/annotated_variants.tsv",
            f"{OUTDIR}/variants/variant_genes.tsv",
            f"{OUTDIR}/variants/variant_annotation_summary.tsv",
            f"{OUTDIR}/variants/variant_qc_metrics.tsv",
            f"{OUTDIR}/variants/sample_missingness.tsv",
            f"{OUTDIR}/variants/variant_type_counts.png",
            f"{OUTDIR}/variants/variant_qual_distribution.png",
            f"{OUTDIR}/variants/variant_dp_distribution.png",
            f"{OUTDIR}/variants/variant_gq_distribution.png",
            f"{OUTDIR}/variants/variant_het_allele_balance.png",
            f"{OUTDIR}/variants/variant_missingness_by_sample.png",
            f"{OUTDIR}/variants/variant_association_hits.tsv",
            f"{OUTDIR}/variants/variant_functional_classes.png",
            f"{OUTDIR}/variants/variant_chromosome_counts.png",
            f"{OUTDIR}/variants/variant_volcano.png",
            f"{OUTDIR}/variants/variant_manhattan.png",
            f"{OUTDIR}/variants/variant_qq.png",
            f"{OUTDIR}/variants/variant_heatmap_top.png",
            f"{OUTDIR}/variants/variant_enrichment.tsv",
            f"{OUTDIR}/variants/variant_enrichment_dotplot.png",
        ]
    )

    if MODE != "variant_only" and ASSAY == "rna":
        targets.extend(
            [
                f"{OUTDIR}/transcripts/stringtie_merged.gtf",
                f"{OUTDIR}/transcripts/gene_counts.tsv",
                f"{OUTDIR}/dge/deseq2_results.tsv",
                f"{OUTDIR}/dge/normalized_counts.tsv",
                f"{OUTDIR}/dge/volcano.png",
                f"{OUTDIR}/dge/heatmap_top50.png",
                f"{OUTDIR}/enrichment/enrichment_results.tsv",
                f"{OUTDIR}/enrichment/enrichment_dotplot.png",
            ]
        )

    if PRS_WEIGHTS:
        targets.extend(
            [
                f"{OUTDIR}/prs/plink.sscore",
                f"{OUTDIR}/prs/prs_scores.tsv",
                f"{OUTDIR}/prs/prs_qc.tsv",
                f"{OUTDIR}/prs/prs_distribution.png",
            ]
        )

    return targets


rule all:
    input:
        core_targets() + [f"{OUTDIR}/docs/run_manifest.txt"]


rule preflight_resources:
    output:
        report=f"{OUTDIR}/docs/preflight_checks.tsv"
    conda:
        "envs/python_plot.yaml"
    params:
        mode=MODE,
        assay=ASSAY,
        read_type=READ_TYPE,
        samplesheet=SAMPLESHEET,
        metadata=METADATA,
        ref_fasta=REF_FASTA,
        ref_gtf=REF_GTF,
        star_index=STAR_INDEX,
        dna_short_aligner=DNA_SHORT_ALIGNER,
        expected_build=REFERENCE.get("genome_build", ""),
        variants_input=EXT_VARIANTS,
        prs_weights=PRS_WEIGHTS,
        prs_weights_format=PRS_WEIGHTS_FORMAT,
        gene_sets=GENE_SETS,
        annotation_method=ANNOTATION_METHOD,
        snpeff_database=SNPEFF_DATABASE,
        enforce_build_match="--enforce-build-match" if VALIDATE_ENFORCE_BUILD_MATCH else "",
        fail_on_warning="--fail-on-warning" if VALIDATE_FAIL_ON_WARNING else "",
    shell:
        "mkdir -p {OUTDIR}/docs && "
        "python workflow/scripts/preflight_validate_resources.py "
        "--mode {params.mode} --assay {params.assay} --read-type {params.read_type} "
        "--samplesheet '{params.samplesheet}' --metadata '{params.metadata}' "
        "--reference-fasta '{params.ref_fasta}' --reference-gtf '{params.ref_gtf}' "
        "--reference-star-index '{params.star_index}' "
        "--dna-short-aligner '{params.dna_short_aligner}' "
        "--expected-build '{params.expected_build}' "
        "--variants-input '{params.variants_input}' "
        "--prs-weights '{params.prs_weights}' --prs-weights-format '{params.prs_weights_format}' "
        "--gene-sets '{params.gene_sets}' "
        "--annotation-method '{params.annotation_method}' --snpeff-database '{params.snpeff_database}' "
        "{params.enforce_build_match} {params.fail_on_warning} "
        "--out-report {output.report}"


rule fastqc_raw:
    input:
        r1=lambda wc: SAMPLE_TO_R1[wc.sample]
    output:
        done=f"{OUTDIR}/qc/raw/{{sample}}.done"
    threads: 2
    conda:
        "envs/qc.yaml"
    params:
        outdir=f"{OUTDIR}/qc/raw"
    run:
        Path(params.outdir).mkdir(parents=True, exist_ok=True)
        r2 = SAMPLE_TO_R2.get(wildcards.sample, "")
        if has_paired_end(wildcards.sample):
            shell("fastqc -t {threads} -o {params.outdir} {input.r1} {r2}")
        else:
            shell("fastqc -t {threads} -o {params.outdir} {input.r1}")
        shell("touch {output.done}")


rule trim_filter_reads:
    input:
        r1=lambda wc: SAMPLE_TO_R1[wc.sample]
    output:
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        html=f"{OUTDIR}/trimmed/{{sample}}.fastp.html",
        json=f"{OUTDIR}/trimmed/{{sample}}.fastp.json"
    threads: THREADS_DEFAULT
    conda:
        "envs/qc.yaml"
    params:
        min_q=MIN_BASEQ,
        min_len=MIN_READ_LENGTH
    run:
        Path(output.r1).parent.mkdir(parents=True, exist_ok=True)
        r2 = SAMPLE_TO_R2.get(wildcards.sample, "")
        if has_paired_end(wildcards.sample):
            shell(
                "fastp "
                "-i {input.r1} -I {r2} "
                "-o {output.r1} -O {output.r2} "
                "--qualified_quality_phred {params.min_q} "
                "--length_required {params.min_len} "
                "--thread {threads} --html {output.html} --json {output.json}"
            )
        else:
            shell(
                "fastp "
                "-i {input.r1} "
                "-o {output.r1} "
                "--qualified_quality_phred {params.min_q} "
                "--length_required {params.min_len} "
                "--thread {threads} --html {output.html} --json {output.json}"
            )
            shell("gzip -nc /dev/null > {output.r2}")


rule fastqc_trimmed:
    input:
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.trimmed.fastq.gz"
    output:
        done=f"{OUTDIR}/qc/trimmed/{{sample}}.done"
    threads: 2
    conda:
        "envs/qc.yaml"
    params:
        outdir=f"{OUTDIR}/qc/trimmed"
    run:
        Path(params.outdir).mkdir(parents=True, exist_ok=True)
        if has_paired_end(wildcards.sample):
            shell("fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}")
        else:
            shell("fastqc -t {threads} -o {params.outdir} {input.r1}")
        shell("touch {output.done}")


rule qc_before_after_report:
    input:
        jsons=expand(f"{OUTDIR}/trimmed/{{sample}}.fastp.json", sample=SAMPLES)
    output:
        table=f"{OUTDIR}/qc/qc_before_after.tsv",
        plot=f"{OUTDIR}/qc/qc_before_after.png"
    conda:
        "envs/python_plot.yaml"
    script:
        "workflow/scripts/qc_compare.py"


rule multiqc:
    input:
        raw_done=expand(f"{OUTDIR}/qc/raw/{{sample}}.done", sample=SAMPLES),
        trimmed_done=expand(f"{OUTDIR}/qc/trimmed/{{sample}}.done", sample=SAMPLES),
        fastp_json=expand(f"{OUTDIR}/trimmed/{{sample}}.fastp.json", sample=SAMPLES)
    output:
        report=f"{OUTDIR}/qc/multiqc_report.html"
    conda:
        "envs/qc.yaml"
    shell:
        "mkdir -p {OUTDIR}/qc && "
        "multiqc -f -o {OUTDIR}/qc {OUTDIR}/qc/raw {OUTDIR}/qc/trimmed {OUTDIR}/trimmed"


if ASSAY == "rna" and READ_TYPE == "short":
    rule build_star_index:
        input:
            ref=REF_FASTA,
            gtf=REF_GTF
        output:
            marker=STAR_INDEX_MARKER
        threads: THREADS_DEFAULT
        conda:
            "envs/alignment.yaml"
        shell:
            "mkdir -p {STAR_INDEX} && "
            "STAR --runThreadN {threads} --runMode genomeGenerate "
            "--genomeDir {STAR_INDEX} --genomeFastaFiles {input.ref} "
            "--sjdbGTFfile {input.gtf} --sjdbOverhang 100 && "
            "touch {output.marker}"


rule align_reads:
    input:
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        ref=REF_FASTA,
        star_idx=star_index_outputs
    output:
        bam=f"{OUTDIR}/alignment/{{sample}}.sorted.bam"
    threads: THREADS_DEFAULT
    conda:
        "envs/alignment.yaml"
    run:
        Path(output.bam).parent.mkdir(parents=True, exist_ok=True)
        r2_exists = has_paired_end(wildcards.sample)

        if ASSAY == "dna" and READ_TYPE == "short":
            dna_short_aligner = select_dna_short_aligner(input.ref)
            if r2_exists:
                if dna_short_aligner == "bwa_mem2":
                    shell(
                        "bwa-mem2 mem -t {threads} {input.ref} {input.r1} {input.r2} "
                        "| samtools sort -@ {threads} -o {output.bam}"
                    )
                elif dna_short_aligner == "bwa":
                    shell(
                        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} "
                        "| samtools sort -@ {threads} -o {output.bam}"
                    )
                elif dna_short_aligner == "minimap2_sr":
                    shell(
                        "minimap2 -ax sr -t {threads} {input.ref} {input.r1} {input.r2} "
                        "| samtools sort -@ {threads} -o {output.bam}"
                    )
                else:
                    raise ValueError(f"Unsupported DNA short-read aligner: {dna_short_aligner}")
            else:
                if dna_short_aligner == "bwa_mem2":
                    shell(
                        "bwa-mem2 mem -t {threads} {input.ref} {input.r1} "
                        "| samtools sort -@ {threads} -o {output.bam}"
                    )
                elif dna_short_aligner == "bwa":
                    shell(
                        "bwa mem -t {threads} {input.ref} {input.r1} "
                        "| samtools sort -@ {threads} -o {output.bam}"
                    )
                elif dna_short_aligner == "minimap2_sr":
                    shell(
                        "minimap2 -ax sr -t {threads} {input.ref} {input.r1} "
                        "| samtools sort -@ {threads} -o {output.bam}"
                    )
                else:
                    raise ValueError(f"Unsupported DNA short-read aligner: {dna_short_aligner}")

        elif ASSAY == "dna" and READ_TYPE == "long":
            if r2_exists:
                shell(
                    "minimap2 -ax {LONG_READ_PRESET_DNA} -t {threads} {input.ref} {input.r1} {input.r2} "
                    "| samtools sort -@ {threads} -o {output.bam}"
                )
            else:
                shell(
                    "minimap2 -ax {LONG_READ_PRESET_DNA} -t {threads} {input.ref} {input.r1} "
                    "| samtools sort -@ {threads} -o {output.bam}"
                )

        elif ASSAY == "rna" and READ_TYPE == "short":
            if not STAR_INDEX:
                raise ValueError("reference.star_index is required for RNA short-read mode")
            outprefix = f"{OUTDIR}/alignment/{wildcards.sample}.STAR."
            if r2_exists:
                shell(
                    "STAR --runThreadN {threads} --genomeDir {STAR_INDEX} "
                    "--readFilesIn {input.r1} {input.r2} --readFilesCommand zcat "
                    "--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {outprefix}"
                )
            else:
                shell(
                    "STAR --runThreadN {threads} --genomeDir {STAR_INDEX} "
                    "--readFilesIn {input.r1} --readFilesCommand zcat "
                    "--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {outprefix}"
                )
            shell("mv {outprefix}Aligned.sortedByCoord.out.bam {output.bam}")

        elif ASSAY == "rna" and READ_TYPE == "long":
            if r2_exists:
                shell(
                    "minimap2 -ax {LONG_READ_PRESET_RNA} -uf -k14 -t {threads} {input.ref} {input.r1} {input.r2} "
                    "| samtools sort -@ {threads} -o {output.bam}"
                )
            else:
                shell(
                    "minimap2 -ax {LONG_READ_PRESET_RNA} -uf -k14 -t {threads} {input.ref} {input.r1} "
                    "| samtools sort -@ {threads} -o {output.bam}"
                )

        else:
            raise ValueError("Unsupported assay/read_type combination")


rule index_bam:
    input:
        bam=f"{OUTDIR}/alignment/{{sample}}.sorted.bam"
    output:
        bai=f"{OUTDIR}/alignment/{{sample}}.sorted.bam.bai"
    threads: THREADS_DEFAULT
    conda:
        "envs/variant.yaml"
    shell:
        "samtools index -@ {threads} {input.bam}"


if MODE != "variant_only" and ASSAY == "dna" and READ_TYPE == "short":
    rule mark_duplicates:
        input:
            bam=f"{OUTDIR}/alignment/{{sample}}.sorted.bam",
            bai=f"{OUTDIR}/alignment/{{sample}}.sorted.bam.bai"
        output:
            bam=f"{OUTDIR}/alignment/{{sample}}.dedup.bam",
            bai=f"{OUTDIR}/alignment/{{sample}}.dedup.bam.bai",
            metrics=f"{OUTDIR}/alignment/{{sample}}.dedup.metrics.txt"
        threads: THREADS_DEFAULT
        conda:
            "envs/variant.yaml"
        params:
            name_sorted=f"{OUTDIR}/alignment/{{sample}}.name_sorted.bam",
            fixmate=f"{OUTDIR}/alignment/{{sample}}.fixmate.bam",
            fixmate_sorted=f"{OUTDIR}/alignment/{{sample}}.fixmate.sorted.bam"
        shell:
            "samtools sort -n -@ {threads} -o {params.name_sorted} {input.bam} && "
            "samtools fixmate -m -@ {threads} {params.name_sorted} {params.fixmate} && "
            "samtools sort -@ {threads} -o {params.fixmate_sorted} {params.fixmate} && "
            "samtools markdup -@ {threads} -r -s {params.fixmate_sorted} {output.bam} 2> {output.metrics} && "
            "samtools index -@ {threads} {output.bam} && "
            "rm -f {params.name_sorted} {params.fixmate} {params.fixmate_sorted}"


rule alignment_qc:
    input:
        bam=lambda wc: variant_bam_for_sample(wc.sample),
        bai=lambda wc: variant_bai_for_sample(wc.sample)
    output:
        flagstat=f"{OUTDIR}/qc/alignment/{{sample}}.flagstat.txt",
        stats=f"{OUTDIR}/qc/alignment/{{sample}}.stats.txt"
    threads: THREADS_DEFAULT
    conda:
        "envs/variant.yaml"
    shell:
        "mkdir -p {OUTDIR}/qc/alignment && "
        "samtools flagstat -@ {threads} {input.bam} > {output.flagstat} && "
        "samtools stats -@ {threads} {input.bam} > {output.stats}"


rule mosdepth_coverage:
    input:
        bam=lambda wc: variant_bam_for_sample(wc.sample),
        bai=lambda wc: variant_bai_for_sample(wc.sample)
    output:
        summary=f"{OUTDIR}/depth/{{sample}}.mosdepth.summary.txt",
        dist=f"{OUTDIR}/depth/{{sample}}.mosdepth.global.dist.txt"
    threads: THREADS_DEFAULT
    conda:
        "envs/variant.yaml"
    params:
        prefix=lambda wc: f"{OUTDIR}/depth/{wc.sample}"
    shell:
        (
            "if [ '{RUN_ENABLE_DEPTH}' != 'True' ]; then "
            "mkdir -p {OUTDIR}/depth && "
            "printf 'chrom\\tlength\\tbases\\tmean\\n' > {output.summary} && "
            ": > {output.dist}; "
            "else "
            "mkdir -p {OUTDIR}/depth && mosdepth -t {threads} -n --fast-mode {params.prefix} {input.bam}; "
            "fi"
        )


rule depth_per_sample:
    input:
        bam=lambda wc: variant_bam_for_sample(wc.sample),
        bai=lambda wc: variant_bai_for_sample(wc.sample)
    output:
        depth=f"{OUTDIR}/depth/{{sample}}.depth.tsv"
    conda:
        "envs/variant.yaml"
    run:
        Path(output.depth).parent.mkdir(parents=True, exist_ok=True)
        if (not RUN_ENABLE_DEPTH) or (not DEPTH_EMIT_PER_BASE_TSV):
            shell(": > {output.depth}")
        else:
            depth_flag = "-a" if DEPTH_EMIT_ALL_POSITIONS else ""
            shell("samtools depth {depth_flag} {input.bam} > {output.depth}")


rule depth_summary:
    input:
        depths=expand(f"{OUTDIR}/depth/{{sample}}.depth.tsv", sample=SAMPLES),
        mosdepth_summaries=expand(f"{OUTDIR}/depth/{{sample}}.mosdepth.summary.txt", sample=SAMPLES),
        mosdepth_dists=expand(f"{OUTDIR}/depth/{{sample}}.mosdepth.global.dist.txt", sample=SAMPLES)
    output:
        table=f"{OUTDIR}/depth/depth_summary.tsv",
        plot=f"{OUTDIR}/depth/depth_distribution.png"
    conda:
        "envs/python_plot.yaml"
    script:
        "workflow/scripts/depth_summary.py"


rule call_snps:
    input:
        bams=variant_bam_outputs,
        bais=variant_bai_outputs,
        ref=REF_FASTA
    output:
        vcf=f"{OUTDIR}/variants/called.raw.vcf.gz",
        tbi=f"{OUTDIR}/variants/called.raw.vcf.gz.tbi"
    threads: THREADS_DEFAULT
    conda:
        "envs/variant.yaml"
    shell:
        "mkdir -p {OUTDIR}/variants && "
        "bcftools mpileup -Ou -f {input.ref} {input.bams} "
        "| bcftools call -mv -Oz -o {output.vcf} && "
        "tabix -f -p vcf {output.vcf}"


rule filter_called_snps:
    input:
        vcf=f"{OUTDIR}/variants/called.raw.vcf.gz",
        tbi=f"{OUTDIR}/variants/called.raw.vcf.gz.tbi"
    output:
        vcf=f"{OUTDIR}/variants/called.filtered.vcf.gz",
        tbi=f"{OUTDIR}/variants/called.filtered.vcf.gz.tbi"
    conda:
        "envs/variant.yaml"
    run:
        if MIN_VAR_QUAL <= 0:
            expr = "(QUAL>=0 || QUAL=\".\")"
        else:
            expr = f"(QUAL>={MIN_VAR_QUAL} || QUAL=\".\")"
        if MIN_VAR_DP > 0:
            expr = f"{expr} && (INFO/DP>={MIN_VAR_DP} || INFO/DP=\".\")"
        shell(
            "bcftools filter -i '{expr}' -Oz -o {output.vcf} {input.vcf} && "
            "tabix -f -p vcf {output.vcf}"
        )


if EXT_VARIANTS:
    rule normalize_external_variants:
        input:
            EXT_VARIANTS
        output:
            vcf=f"{OUTDIR}/variants/external.normalized.vcf.gz",
            tbi=f"{OUTDIR}/variants/external.normalized.vcf.gz.tbi"
        conda:
            "envs/variant.yaml"
        run:
            Path(f"{OUTDIR}/variants").mkdir(parents=True, exist_ok=True)
            ext = str(input[0]).lower()
            tmp_vcf = f"{OUTDIR}/variants/external.input.vcf"
            tmp_vcfgz = f"{OUTDIR}/variants/external.input.vcf.gz"
            tmp_norm_table = f"{OUTDIR}/variants/external.input.normalized.tsv"
            tmp_norm_report = f"{OUTDIR}/variants/external.input.field_report.tsv"

            if ext.endswith(".csv") or ext.endswith(".tsv"):
                shell(
                    "python workflow/scripts/csv_to_vcf.py "
                    "--input {input[0]} --output {tmp_vcf} "
                    "--ref-fasta {REF_FASTA} "
                    "--normalized-table {tmp_norm_table} "
                    "--report {tmp_norm_report}"
                )
                shell("bgzip -f -c {tmp_vcf} > {tmp_vcfgz}")
                shell("tabix -f -p vcf {tmp_vcfgz}")
            elif ext.endswith(".vcf"):
                shell("bgzip -f -c {input[0]} > {tmp_vcfgz}")
                shell("tabix -f -p vcf {tmp_vcfgz}")
            elif ext.endswith(".vcf.gz"):
                shell("bcftools view -Oz -o {tmp_vcfgz} {input[0]}")
                shell("tabix -f -p vcf {tmp_vcfgz}")
            else:
                raise ValueError("paths.variants_input must end with .csv, .tsv, .vcf, or .vcf.gz")

            shell(
                "if bcftools norm -m -any -f {REF_FASTA} {tmp_vcfgz} -Oz -o {output.vcf}; then "
                "echo 'Reference-based normalization completed.'; "
                "else "
                "echo 'Reference mismatch detected. Falling back to split-only normalization.' >&2; "
                "bcftools norm -m -any {tmp_vcfgz} -Oz -o {output.vcf}; "
                "fi"
            )
            shell("tabix -f -p vcf {output.vcf}")


    rule filter_external_variants:
        input:
            vcf=f"{OUTDIR}/variants/external.normalized.vcf.gz",
            tbi=f"{OUTDIR}/variants/external.normalized.vcf.gz.tbi"
        output:
            vcf=f"{OUTDIR}/variants/external.filtered.vcf.gz",
            tbi=f"{OUTDIR}/variants/external.filtered.vcf.gz.tbi"
        conda:
            "envs/variant.yaml"
        run:
            filter_flag = "-f PASS,." if REQUIRE_FILTER_PASS else ""
            if MIN_VAR_QUAL <= 0:
                expr = "(QUAL>=0 || QUAL=\".\")"
            else:
                expr = f"(QUAL>={MIN_VAR_QUAL} || QUAL=\".\")"
            if MIN_VAR_DP > 0:
                expr = f"{expr} && (INFO/DP>={MIN_VAR_DP} || INFO/DP=\".\")"

            shell(
                "bcftools view {filter_flag} -v snps -i '{expr}' -Oz -o {output.vcf} {input.vcf} && "
                "tabix -f -p vcf {output.vcf}"
            )


def final_variant_source(wildcards=None):
    if MODE == "variant_only":
        return f"{OUTDIR}/variants/external.filtered.vcf.gz"
    if EXT_VARIANTS and PRS_USE_EXTERNAL:
        return f"{OUTDIR}/variants/external.filtered.vcf.gz"
    return f"{OUTDIR}/variants/called.filtered.vcf.gz"


def final_variant_index_source(wildcards=None):
    return f"{final_variant_source()}.tbi"


rule final_variants:
    input:
        vcf=final_variant_source,
        tbi=final_variant_index_source
    output:
        vcf=f"{OUTDIR}/variants/final.filtered.vcf.gz",
        tbi=f"{OUTDIR}/variants/final.filtered.vcf.gz.tbi"
    shell:
        "mkdir -p {OUTDIR}/variants && cp {input.vcf} {output.vcf} && cp {input.tbi} {output.tbi}"


if ANNOTATION_METHOD == "snpeff":
    rule snpeff_annotate_variants:
        input:
            vcf=f"{OUTDIR}/variants/final.filtered.vcf.gz",
            tbi=f"{OUTDIR}/variants/final.filtered.vcf.gz.tbi"
        output:
            vcf=f"{OUTDIR}/variants/final.snpeff.vcf.gz",
            tbi=f"{OUTDIR}/variants/final.snpeff.vcf.gz.tbi",
            html=f"{OUTDIR}/variants/snpeff_summary.html"
        threads: THREADS_DEFAULT
        conda:
            "envs/variant.yaml"
        params:
            db=SNPEFF_DATABASE
        shell:
            "snpEff -canon -stats {output.html} {params.db} {input.vcf} | bgzip -c > {output.vcf} && "
            "tabix -f -p vcf {output.vcf}"


rule vcf_to_variant_bed:
    input:
        vcf=annotation_vcf_source,
        tbi=annotation_vcf_index_source
    output:
        bed=f"{OUTDIR}/variants/intermediate/variants.bed",
        table=f"{OUTDIR}/variants/intermediate/variants.base.tsv"
    conda:
        "envs/variant.yaml"
    shell:
        "mkdir -p {OUTDIR}/variants/intermediate && "
        "python workflow/scripts/vcf_to_variant_bed.py "
        "--vcf {input.vcf} --out-bed {output.bed} --out-table {output.table}"


rule gtf_to_feature_bed:
    input:
        gtf=REF_GTF
    output:
        bed=f"{OUTDIR}/variants/intermediate/reference_features.bed"
    conda:
        "envs/variant.yaml"
    shell:
        "mkdir -p {OUTDIR}/variants/intermediate && "
        "python workflow/scripts/gtf_to_features_bed.py "
        "--gtf {input.gtf} --out-bed {output.bed}"


rule intersect_variant_features:
    input:
        variants=f"{OUTDIR}/variants/intermediate/variants.bed",
        features=f"{OUTDIR}/variants/intermediate/reference_features.bed"
    output:
        overlaps=f"{OUTDIR}/variants/intermediate/variant_feature_overlaps.tsv"
    conda:
        "envs/variant.yaml"
    shell:
        "bedtools intersect -a {input.variants} -b {input.features} -wa -wb > {output.overlaps}"


rule annotate_variants:
    input:
        table=f"{OUTDIR}/variants/intermediate/variants.base.tsv",
        overlaps=f"{OUTDIR}/variants/intermediate/variant_feature_overlaps.tsv"
    output:
        annotated=f"{OUTDIR}/variants/annotated_variants.tsv",
        genes=f"{OUTDIR}/variants/variant_genes.tsv",
        summary=f"{OUTDIR}/variants/variant_annotation_summary.tsv"
    conda:
        "envs/variant.yaml"
    shell:
        "python workflow/scripts/aggregate_variant_annotations.py "
        "--variants {input.table} --overlaps {input.overlaps} "
        "--out-annotated {output.annotated} --out-genes {output.genes} --out-summary {output.summary}"


rule variant_qc_stats:
    input:
        vcf=f"{OUTDIR}/variants/final.filtered.vcf.gz",
        tbi=f"{OUTDIR}/variants/final.filtered.vcf.gz.tbi"
    output:
        metrics=f"{OUTDIR}/variants/variant_qc_metrics.tsv",
        missingness=f"{OUTDIR}/variants/sample_missingness.tsv",
        type_plot=f"{OUTDIR}/variants/variant_type_counts.png",
        qual_plot=f"{OUTDIR}/variants/variant_qual_distribution.png",
        dp_plot=f"{OUTDIR}/variants/variant_dp_distribution.png",
        gq_plot=f"{OUTDIR}/variants/variant_gq_distribution.png",
        ab_plot=f"{OUTDIR}/variants/variant_het_allele_balance.png",
        miss_plot=f"{OUTDIR}/variants/variant_missingness_by_sample.png"
    conda:
        "envs/python_plot.yaml"
    params:
        max_gt_samples=VAR_QC_MAX_GT_SAMPLES
    shell:
        "python workflow/scripts/variant_qc_stats.py "
        "--vcf {input.vcf} "
        "--out-metrics {output.metrics} "
        "--out-missingness {output.missingness} "
        "--out-type-plot {output.type_plot} "
        "--out-qual-plot {output.qual_plot} "
        "--out-dp-plot {output.dp_plot} "
        "--out-gq-plot {output.gq_plot} "
        "--out-ab-plot {output.ab_plot} "
        "--out-missingness-plot {output.miss_plot} "
        "--max-genotype-samples {params.max_gt_samples}"


rule variant_visuals:
    input:
        annotated=f"{OUTDIR}/variants/annotated_variants.tsv",
        vcf=f"{OUTDIR}/variants/final.filtered.vcf.gz",
        tbi=f"{OUTDIR}/variants/final.filtered.vcf.gz.tbi"
    output:
        assoc=f"{OUTDIR}/variants/variant_association_hits.tsv",
        feature_plot=f"{OUTDIR}/variants/variant_functional_classes.png",
        chrom_plot=f"{OUTDIR}/variants/variant_chromosome_counts.png",
        volcano=f"{OUTDIR}/variants/variant_volcano.png",
        manhattan=f"{OUTDIR}/variants/variant_manhattan.png",
        qq=f"{OUTDIR}/variants/variant_qq.png",
        heatmap=f"{OUTDIR}/variants/variant_heatmap_top.png"
    params:
        pval=VAR_ASSOC_PVALUE,
        effect=VAR_ASSOC_ABS_EFFECT,
        topn=VAR_HEATMAP_TOP_N,
        analysis_mode=ACTIVE_VARIANT_DATA_MODE
    conda:
        "envs/python_plot.yaml"
    shell:
        "python workflow/scripts/variant_visualization.py "
        "--annotated {input.annotated} --vcf {input.vcf} "
        "--analysis-mode {params.analysis_mode} "
        "--p-threshold {params.pval} --effect-threshold {params.effect} "
        "--heatmap-top-n {params.topn} "
        "--out-assoc {output.assoc} "
        "--out-feature-plot {output.feature_plot} --out-chrom-plot {output.chrom_plot} "
        "--out-volcano {output.volcano} --out-manhattan {output.manhattan} --out-qq {output.qq} "
        "--out-heatmap {output.heatmap}"


rule variant_enrichment:
    input:
        genes=f"{OUTDIR}/variants/variant_genes.tsv"
    output:
        table=f"{OUTDIR}/variants/variant_enrichment.tsv",
        plot=f"{OUTDIR}/variants/variant_enrichment_dotplot.png"
    params:
        genesets=GENE_SETS,
        qvalue=ENRICH_QVALUE
    conda:
        "envs/r_stats.yaml"
    shell:
        "Rscript workflow/scripts/variant_enrichment.R "
        "--genes {input.genes} --genesets '{params.genesets}' "
        "--qvalue-threshold {params.qvalue} "
        "--out-table {output.table} --out-plot {output.plot}"


if MODE != "variant_only" and ASSAY == "rna":
    rule stringtie_assemble:
        input:
            bam=f"{OUTDIR}/alignment/{{sample}}.sorted.bam",
            bai=f"{OUTDIR}/alignment/{{sample}}.sorted.bam.bai",
            gtf=REF_GTF
        output:
            gtf=f"{OUTDIR}/transcripts/{{sample}}.stringtie.gtf",
            abund=f"{OUTDIR}/transcripts/{{sample}}.gene_abund.tsv"
        threads: THREADS_DEFAULT
        conda:
            "envs/rna.yaml"
        shell:
            "mkdir -p {OUTDIR}/transcripts && "
            "stringtie {input.bam} -p {threads} -G {input.gtf} -o {output.gtf} -A {output.abund}"


    rule stringtie_merge:
        input:
            gtfs=expand(f"{OUTDIR}/transcripts/{{sample}}.stringtie.gtf", sample=SAMPLES),
            gtf=REF_GTF
        output:
            merged=f"{OUTDIR}/transcripts/stringtie_merged.gtf"
        threads: THREADS_DEFAULT
        conda:
            "envs/rna.yaml"
        shell:
            "stringtie --merge -p {threads} -G {input.gtf} -o {output.merged} {input.gtfs}"


    rule quantify_gene_counts:
        input:
            bams=expand(f"{OUTDIR}/alignment/{{sample}}.sorted.bam", sample=SAMPLES),
            bai=expand(f"{OUTDIR}/alignment/{{sample}}.sorted.bam.bai", sample=SAMPLES),
            gtf=f"{OUTDIR}/transcripts/stringtie_merged.gtf"
        output:
            counts=f"{OUTDIR}/transcripts/gene_counts.tsv",
            summary=f"{OUTDIR}/transcripts/gene_counts.summary"
        threads: THREADS_DEFAULT
        conda:
            "envs/rna.yaml"
        shell:
            "mkdir -p {OUTDIR}/transcripts && "
            "featureCounts -T {threads} -a {input.gtf} -o {output.counts} {input.bams} && "
            "cp {output.counts}.summary {output.summary}"


    rule differential_expression:
        input:
            counts=f"{OUTDIR}/transcripts/gene_counts.tsv",
            metadata=METADATA
        output:
            res=f"{OUTDIR}/dge/deseq2_results.tsv",
            norm=f"{OUTDIR}/dge/normalized_counts.tsv",
            volcano=f"{OUTDIR}/dge/volcano.png",
            heatmap=f"{OUTDIR}/dge/heatmap_top50.png"
        params:
            formula=config.get("deseq2", {}).get("design_formula", "~ condition"),
            contrast=config.get("deseq2", {}).get("contrast", "condition,case,control"),
            padj=DE_PADJ,
            lfc=DE_ABS_LOG2FC
        conda:
            "envs/r_stats.yaml"
        shell:
            "mkdir -p {OUTDIR}/dge && "
            "Rscript workflow/scripts/deseq2_analysis.R "
            "--counts {input.counts} --metadata {input.metadata} "
            "--formula '{params.formula}' --contrast '{params.contrast}' "
            "--padj-threshold {params.padj} --lfc-threshold {params.lfc} "
            "--out-results {output.res} --out-normalized {output.norm} "
            "--out-volcano {output.volcano} --out-heatmap {output.heatmap}"


    rule enrichment:
        input:
            de=f"{OUTDIR}/dge/deseq2_results.tsv"
        output:
            table=f"{OUTDIR}/enrichment/enrichment_results.tsv",
            plot=f"{OUTDIR}/enrichment/enrichment_dotplot.png"
        params:
            genesets=GENE_SETS,
            padj=DE_PADJ,
            lfc=DE_ABS_LOG2FC,
            qvalue=ENRICH_QVALUE
        conda:
            "envs/r_stats.yaml"
        shell:
            "mkdir -p {OUTDIR}/enrichment && "
            "Rscript workflow/scripts/enrichment_analysis.R "
            "--de {input.de} --genesets '{params.genesets}' "
            "--padj-threshold {params.padj} --lfc-threshold {params.lfc} "
            "--qvalue-threshold {params.qvalue} "
            "--out-table {output.table} --out-plot {output.plot}"


if PRS_WEIGHTS:
    rule prs_score:
        input:
            vcf=f"{OUTDIR}/variants/final.filtered.vcf.gz",
            tbi=f"{OUTDIR}/variants/final.filtered.vcf.gz.tbi",
            weights=PRS_WEIGHTS
        output:
            sscore=f"{OUTDIR}/prs/plink.sscore",
            table=f"{OUTDIR}/prs/prs_scores.tsv",
            qc=f"{OUTDIR}/prs/prs_qc.tsv",
            plot=f"{OUTDIR}/prs/prs_distribution.png"
        params:
            wfmt=PRS_WEIGHTS_FORMAT,
            prefer_ds="--prefer-ds" if PRS_PREFER_DS else "",
            expected_build=REFERENCE.get("genome_build", ""),
            allow_ambiguous="--allow-ambiguous-strand" if PRS_ALLOW_AMBIGUOUS_STRAND else ""
        conda:
            "envs/prs.yaml"
        shell:
            "mkdir -p {OUTDIR}/prs && "
            "python workflow/scripts/prs_from_vcf.py "
            "--vcf {input.vcf} --weights {input.weights} "
            "--weights-format {params.wfmt} "
            "--expected-build '{params.expected_build}' "
            "--out-sscore {output.sscore} --out-qc {output.qc} {params.prefer_ds} {params.allow_ambiguous} && "
            "python workflow/scripts/prs_report.py "
            "--sscore {output.sscore} --out-table {output.table} --out-plot {output.plot}"


rule run_manifest:
    input:
        core_targets()
    output:
        manifest=f"{OUTDIR}/docs/run_manifest.txt"
    run:
        lines = [
            "Bioinformatics Pipeline Run Manifest",
            f"timestamp: {datetime.utcnow().isoformat()}Z",
            f"mode_raw: {MODE_RAW}",
            f"mode: {MODE}",
            f"variant_data_mode_config: {VARIANT_DATA_MODE}",
            f"variant_data_mode_active: {ACTIVE_VARIANT_DATA_MODE}",
            f"assay: {ASSAY}",
            f"read_type: {READ_TYPE}",
            f"run_enable_depth: {RUN_ENABLE_DEPTH}",
            f"samplesheet: {SAMPLESHEET}",
            f"metadata: {METADATA}",
            f"reference_fasta: {REF_FASTA}",
            f"reference_gtf: {REF_GTF}",
            f"external_variants: {EXT_VARIANTS or 'none'}",
            f"prs_weights: {PRS_WEIGHTS or 'none'}",
            f"prs_weights_format: {PRS_WEIGHTS_FORMAT}",
            f"prs_prefer_ds: {PRS_PREFER_DS}",
            f"prs_allow_ambiguous_strand: {PRS_ALLOW_AMBIGUOUS_STRAND}",
            f"annotation_method: {ANNOTATION_METHOD}",
            f"snpeff_database: {SNPEFF_DATABASE or 'none'}",
            f"dna_short_aligner: {DNA_SHORT_ALIGNER}",
            f"depth_emit_all_positions: {DEPTH_EMIT_ALL_POSITIONS}",
            f"depth_emit_per_base_tsv: {DEPTH_EMIT_PER_BASE_TSV}",
            f"validation_enforce_build_match: {VALIDATE_ENFORCE_BUILD_MATCH}",
            f"validation_fail_on_warning: {VALIDATE_FAIL_ON_WARNING}",
            f"preflight_report: {OUTDIR}/docs/preflight_checks.tsv",
            f"final_variant_source: {final_variant_source()}",
            f"de_padj_threshold: {DE_PADJ}",
            f"de_abs_log2fc_threshold: {DE_ABS_LOG2FC}",
            f"enrichment_qvalue_threshold: {ENRICH_QVALUE}",
            f"min_variant_depth_threshold: {MIN_VAR_DP}",
            f"require_filter_pass: {REQUIRE_FILTER_PASS}",
            f"variant_assoc_pvalue_threshold: {VAR_ASSOC_PVALUE}",
            f"variant_assoc_abs_effect_threshold: {VAR_ASSOC_ABS_EFFECT}",
            f"variant_heatmap_top_n: {VAR_HEATMAP_TOP_N}",
            f"gene_sets: {GENE_SETS or 'none'}",
        ]
        Path(output.manifest).parent.mkdir(parents=True, exist_ok=True)
        with open(output.manifest, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
