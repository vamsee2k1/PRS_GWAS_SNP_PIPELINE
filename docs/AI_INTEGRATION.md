# AI Integration Guide

This repository now includes an optional AI layer designed for reproducible, offline-friendly analysis support.

## Modules

1. Results Explainer
- Script: `workflow/scripts/ai_explainer.py`
- Output: `results/docs/ai/ai_explainer.md`
- Purpose: generate a structured narrative summary from QC, association, enrichment, and optional PRS/anomaly/prioritization tables.

2. QC Anomaly Detector
- Script: `workflow/scripts/ai_qc_anomaly.py`
- Outputs:
  - `results/docs/ai/ai_qc_anomalies.tsv`
  - `results/docs/ai/ai_qc_anomaly_report.md`
- Purpose: detect potential quality concerns (Ts/Tv outliers, low retention, low depth, preflight warnings/errors).

3. Variant Prioritization
- Script: `workflow/scripts/ai_variant_prioritization.py`
- Outputs:
  - `results/docs/ai/ai_variant_prioritization.tsv`
  - `results/docs/ai/ai_variant_prioritization_top.png`
  - `results/docs/ai/ai_variant_prioritization_report.md`
- Purpose: rank variants using association strength, effect size, functional class, clinical signal, and priority-gene context.

## Configuration

Configure in `config/config.yaml`:

```yaml
ai:
  enabled: true
  explainer: true
  qc_anomaly: true
  variant_prioritization: true
  top_n: 50
  trait_context: "Alzheimer disease example context"
  priority_genes:
    - APOE
    - TOMM40
    - TREM2
```

## Design Notes

- The AI layer is deterministic and does not require external API calls.
- Core scientific outputs remain unchanged; AI files are additive.
- Run manifest captures AI settings (`ai_*` fields) for reproducibility.

