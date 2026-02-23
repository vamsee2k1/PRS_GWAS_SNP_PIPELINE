# Statistical Threshold Defaults Used in This Pipeline

These defaults are implemented in `config/config.yaml` and are configurable.

## 1. Differential Expression (RNA)

- `thresholds.de_padj: 0.05`
- `thresholds.de_abs_log2fc: 1.0`

Rationale:
- BH-adjusted p-values (`padj`) control false discovery rate.
- `|log2FC| >= 1` is a common practical effect-size filter for biological interpretability.

References:
- Love MI, Huber W, Anders S. DESeq2. Genome Biology (2014). https://doi.org/10.1186/s13059-014-0550-8
- Benjamini Y, Hochberg Y. FDR control. JRSS B (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

## 2. Enrichment

- `thresholds.enrichment_qvalue: 0.05`

Rationale:
- Enrichment uses BH-adjusted p-values and reports terms below q-value/FDR threshold.

References:
- Benjamini Y, Hochberg Y. JRSS B (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Wu T et al. clusterProfiler 4.0. Innovation (2021). https://doi.org/10.1016/j.xinn.2021.100141

## 3. Variant Association Visuals and Hit Table

- `thresholds.variant_assoc_pvalue: 5e-8`
- `thresholds.variant_assoc_abs_effect: 0.0` (set >0 to enforce effect-size floor)

Rationale:
- `5e-8` is the conventional genome-wide significance threshold for GWAS-scale testing.
- Effect-size threshold is left open by default because acceptable minimum effect depends on trait architecture and study design.

References:
- Jannot AS et al. J Clin Epidemiol (2015). https://doi.org/10.1016/j.jclinepi.2015.01.001
- Cheruiyot EK et al. Genetics (2025). https://doi.org/10.1093/genetics/iyaf056

## Notes

- Thresholds are not universal constants. Adjust for cohort size, phenotype model, and multiple-testing burden.
- For very large GWAS cohorts, stricter significance thresholds can be justified.
