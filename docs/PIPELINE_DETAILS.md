# Pipeline Details (Step-by-Step)

## Pipeline Flowcharts

### 1) Preflight + Mode Routing

```mermaid
flowchart TD
%% Nodes
A("Start"):::green
B("Load Config <br> config.yaml + override"):::orange
C("Resolve Runtime Settings <br> mode / assay / read_type / variant_data_mode"):::blue
D("Build Snakemake DAG"):::blue
E("Preflight Validation <br> preflight_validate_resources.py"):::yellow
F{"Preflight Status?"}:::yellow
ZERR("Stop Run <br> preflight_checks.tsv"):::pink
G{"run.mode"}:::yellow

H("FULL MODE branch"):::pink
I("VARIANT_ONLY branch"):::purple
J("Shared Interpretation Layer"):::green
K("Run Manifest + Final Outputs"):::green
L("End"):::green

%% Edges
A --> B --> C --> D --> E --> F
F -- Error --> ZERR
F -- OK / Warning --> G
G -- full --> H
G -- variant_only --> I
H --> J
I --> J
J --> K --> L
E --> K

%% Styling
classDef green fill:#B2DFDB,stroke:#00897B,stroke-width:2px;
classDef orange fill:#FFE0B2,stroke:#FB8C00,stroke-width:2px;
classDef blue fill:#BBDEFB,stroke:#1976D2,stroke-width:2px;
classDef yellow fill:#FFF9C4,stroke:#FBC02D,stroke-width:2px;
classDef pink fill:#F8BBD0,stroke:#C2185B,stroke-width:2px;
classDef purple fill:#E1BEE7,stroke:#8E24AA,stroke-width:2px;
```

### 2) Full Mode (FASTQ to Variant Calling / RNA Branch)

```mermaid
flowchart TD
%% Nodes
A("FULL MODE Inputs <br> samplesheet + metadata + FASTQ + references"):::pink
B("Raw QC <br> FastQC"):::orange
C("Trim / Filter Reads <br> fastp"):::orange
D("Post-trim QC + Aggregation <br> FastQC + MultiQC + qc_compare.py"):::orange
E{"assay + read_type"}:::yellow

F{"DNA short aligner <br> auto / bwa_mem2 / bwa / minimap2_sr"}:::yellow
G("Align DNA short <br> bwa-mem2 mem + samtools sort"):::blue
H("Align DNA short <br> bwa mem + samtools sort"):::blue
I("Align DNA short <br> minimap2 -ax sr + samtools sort"):::blue

J("Align DNA long <br> minimap2 DNA preset"):::blue
K{"STAR index ready?"}:::yellow
L("Build STAR index"):::purple
M("Align RNA short <br> STAR"):::blue
N("Align RNA long <br> minimap2 splice preset"):::blue

O("Sorted BAM"):::green
P("Index BAM <br> samtools index"):::blue
Q{"DNA short mode?"}:::yellow
R("Mark duplicates + index <br> samtools markdup"):::purple
S("Final BAM for downstream"):::green

T("Alignment QC <br> samtools flagstat + stats"):::orange
U{"run.enable_depth?"}:::yellow
V("Depth branch <br> mosdepth + samtools depth"):::purple
W("Depth skipped / placeholders"):::orange

X{"DNA assay?"}:::yellow
Y("Variant Calling <br> bcftools mpileup + call"):::pink
Z("Variant Filtering <br> bcftools filter / view"):::pink
AA("called.filtered.vcf.gz"):::green

AB("RNA branch entry"):::purple
AC("StringTie assembly / merge"):::purple
AD("featureCounts gene quantification"):::purple
AE("DESeq2 differential expression"):::purple
AF("RNA enrichment analysis"):::purple

AG("Output to shared interpretation layer"):::green

%% Edges
A --> B --> C --> D --> E

E -- DNA short --> F
F -- auto + bwa-mem2 --> G
F -- auto fallback / bwa --> H
F -- minimap2_sr --> I

E -- DNA long --> J
E -- RNA short --> K
K -- No --> L --> M
K -- Yes --> M
E -- RNA long --> N

G --> O
H --> O
I --> O
J --> O
M --> O
N --> O

O --> P --> Q
Q -- Yes --> R --> S
Q -- No --> S

S --> T
S --> U
U -- Yes --> V
U -- No --> W

S --> X
X -- Yes --> Y --> Z --> AA --> AG
X -- No (RNA) --> AB --> AC --> AD --> AE --> AF --> AG

%% Styling
classDef green fill:#B2DFDB,stroke:#00897B,stroke-width:2px;
classDef orange fill:#FFE0B2,stroke:#FB8C00,stroke-width:2px;
classDef blue fill:#BBDEFB,stroke:#1976D2,stroke-width:2px;
classDef yellow fill:#FFF9C4,stroke:#FBC02D,stroke-width:2px;
classDef pink fill:#F8BBD0,stroke:#C2185B,stroke-width:2px;
classDef purple fill:#E1BEE7,stroke:#8E24AA,stroke-width:2px;
```

### 3) Variant-Only + Shared Interpretation + PRS

```mermaid
flowchart TD
%% Nodes
A("VARIANT_ONLY Input <br> VCF / VCF.GZ / CSV / TSV"):::purple
B{"Input type"}:::yellow
C("Normalize external variants <br> bcftools norm/view + tabix"):::blue
D("CSV/TSV -> VCF conversion <br> csv_to_vcf.py"):::purple
E("Alias mapping + REF derivation <br> diagnostics reports"):::yellow
F("Filter external variants <br> bcftools view + tabix"):::blue
G("external.filtered.vcf.gz"):::green

H{"Final variants source select"}:::yellow
I("final.filtered.vcf.gz"):::green

J{"annotation.method"}:::yellow
K("Overlap mode <br> use final VCF directly"):::orange
L("snpEff annotation + bgzip + tabix"):::purple
M("final.snpeff.vcf.gz"):::green

N("VCF -> variant BED + base table"):::blue
O("GTF -> feature BED"):::blue
P("Variant-feature overlap <br> bedtools intersect"):::blue
Q("Aggregate variant annotations"):::blue

R("annotated_variants.tsv"):::green
S("variant_genes.tsv"):::green
T("variant_annotation_summary.tsv"):::green

U("Variant QC metrics + plots"):::orange
V{"variant_data_mode"}:::yellow
W("GWAS-style visuals <br> Manhattan / QQ / Volcano / Heatmap"):::pink
X("VCF interpretation visuals"):::pink

Y("Variant enrichment analysis"):::purple
Z("variant_enrichment.tsv + dotplot"):::green

AA{"PRS weights configured?"}:::yellow
AB("Skip PRS branch"):::orange
AC("PRS scoring + harmonization"):::purple
AD("PRS reports + QC + distribution"):::purple
AE{"Genotype-bearing VCF?"}:::yellow
AF("No sample-level PRS <br> summary-style inputs"):::orange
AG("Sample-level PRS scores <br> if loci overlap"):::green

AH("Shared outputs to manifest"):::green

%% Edges
A --> B
B -- VCF/VCF.GZ --> C
B -- CSV/TSV --> D --> E --> C
C --> F --> G

G --> H
H --> I

I --> J
J -- overlap --> K
J -- snpeff --> L --> M --> K

K --> N
O --> P
N --> P --> Q

Q --> R
Q --> S
Q --> T

I --> U
R --> V
I --> V
V -- gwas_summary --> W
V -- vcf_interpretation --> X

S --> Y --> Z

I --> AA
AA -- No --> AB
AA -- Yes --> AC --> AD
AC --> AE
AE -- No --> AF
AE -- Yes --> AG

R --> AH
S --> AH
T --> AH
U --> AH
W --> AH
X --> AH
Z --> AH
AD --> AH
AB --> AH

%% Styling
classDef green fill:#B2DFDB,stroke:#00897B,stroke-width:2px;
classDef orange fill:#FFE0B2,stroke:#FB8C00,stroke-width:2px;
classDef blue fill:#BBDEFB,stroke:#1976D2,stroke-width:2px;
classDef yellow fill:#FFF9C4,stroke:#FBC02D,stroke-width:2px;
classDef pink fill:#F8BBD0,stroke:#C2185B,stroke-width:2px;
classDef purple fill:#E1BEE7,stroke:#8E24AA,stroke-width:2px;
```

## Step 1: Pre-processing (`FastQC`)

- Runs FastQC on raw reads.
- Outputs per-sample HTML/ZIP reports in `results/qc/raw/`.

## Step 2: Quality Control Before/After + Read Depth

- `fastp` performs filtering/trimming.
- FastQC reruns on trimmed reads.
- `workflow/scripts/qc_compare.py` summarizes before/after read retention.
- `samtools depth` computes per-position depth.
- `workflow/scripts/depth_summary.py` creates depth summary and distribution plots.
- `samtools flagstat` and `samtools stats` create per-sample alignment QC summaries.
- `mosdepth` creates per-sample coverage summaries.

Outputs:
- `results/qc/qc_before_after.tsv`
- `results/qc/qc_before_after.png`
- `results/depth/depth_summary.tsv`
- `results/depth/depth_distribution.png`

## Step 3: Filtering

- Read-level filtering by fastp (quality + length thresholds).
- Variant-level filtering by bcftools (`FILTER=PASS`, `QUAL`, optional `INFO/DP` threshold).

## Step 4: Alignment to Reference Genome

- DNA short reads: `bwa-mem2`
- DNA long reads: `minimap2`
- RNA short reads: `STAR`
- RNA long reads: `minimap2` splice mode

Outputs:
- `results/alignment/*.sorted.bam`
- `results/alignment/*.sorted.bam.bai`

## Step 5: Transcript + SNP Identification + Functional Variant Annotation

RNA branch:
- Transcript assembly: `StringTie`
- Merged transcript model: `stringtie --merge`
- Gene-level quantification: `featureCounts`

Variant branch:
- Variant calling: `bcftools mpileup + call`
- Short-read DNA mode marks duplicates with `samtools markdup` before calling.
- Variant filtering: `bcftools filter`
- Biological variant QC summary: SNP/indel counts, Ts/Tv, QUAL/DP/GQ distributions, allele balance, per-sample missingness
- Functional annotation: genomic feature overlap against GTF (`CDS`, `UTR`, `exon`, `intronic_or_gene`, `intergenic`)
- Functional annotation merges overlap-based mapping with reported VCF gene tags when available (`GENEINFO`, `GENE`, etc.)
- Optional primary consequence annotation via `snpEff` is supported (`annotation.method=snpeff`).
- Gene mapping table for downstream biology interpretation
- Variant plots: feature distribution, chromosome counts, volcano, Manhattan, genotype heatmap
  - Sample VCF QC family: SNP/indel counts, Ts/Tv, QUAL/DP/GQ distributions, allele balance, per-sample missingness
  - GWAS family: Manhattan/volcano/QQ when reported association p-values are available
  - If p-value/effect are missing, volcano/Manhattan use proxy prioritization (QUAL/AF and ClinVar-style clinical significance where available).
  - GWAS mode additionally generates a QQ plot (reported p-values only).
  - If genotype columns are absent, heatmap falls back to a variant-feature heatmap.

Outputs:
- `results/transcripts/`
- `results/variants/`

## Step 6: Differential Gene Expression

- `DESeq2` on gene count matrix.
- Produces full results table, normalized counts, volcano plot, top-gene heatmap.

Outputs:
- `results/dge/deseq2_results.tsv`
- `results/dge/normalized_counts.tsv`
- `results/dge/volcano.png`
- `results/dge/heatmap_top50.png`

## Step 7: Enrichment Analysis

- `clusterProfiler::enricher` using provided GMT pathways.
- Generates enrichment table and dotplot.

Outputs:
- `results/enrichment/enrichment_results.tsv`
- `results/enrichment/enrichment_dotplot.png`
- `results/variants/variant_enrichment.tsv`
- `results/variants/variant_enrichment_dotplot.png`

## Step 8: PRS Scoring

- `workflow/scripts/prs_from_vcf.py` computes per-sample PRS from filtered VCF + PRS model file.
- Supports official PGS Catalog scoring files and simple TSV weights.
- Includes allele harmonization and optional dosage (`DS`) parsing.
- PRS QC includes build-compatibility checks, rsID/position match counts, strand flips, and skipped ambiguous SNP counts.
- `workflow/scripts/prs_report.py` produces table + distribution plot.

Outputs:
- `results/prs/plink.sscore`
- `results/prs/prs_scores.tsv`
- `results/prs/prs_qc.tsv`
- `results/prs/prs_distribution.png`

## Destination and Documentation

- Result destination folder: `results/`
- Detailed pipeline documentation: `docs/`
- Run-specific manifest: `results/docs/run_manifest.txt`

## Default Statistical Thresholds

- Differential expression significance: `padj < 0.05` and `|log2FC| >= 1.0`
- Enrichment significance: `q-value (BH-adjusted p) < 0.05`
- Variant association plots/tables: `p-value <= 5e-8` (genome-wide significance), configurable in `config.yaml`
