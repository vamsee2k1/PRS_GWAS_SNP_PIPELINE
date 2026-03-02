# Biological Significance Deep Dive: PRS, Differential Expression, and Integrated Multi-Mode Genomics

## 1) Why This Pipeline Is Biologically Meaningful

Most disease genomics workflows answer only one layer of biology:

- GWAS/variant workflows estimate inherited risk architecture.
- RNA differential-expression workflows estimate dynamic molecular state.
- Pathway workflows estimate systems-level convergence.

This pipeline is biologically significant because it lets users connect those layers in one reproducible framework:

- inherited susceptibility (`variant_only`, GWAS + VCF interpretation),
- molecular activity (`full` RNA/DNA branches with optional DE),
- pathway-level interpretation (enrichment from variant-mapped genes and/or DE genes),
- optional PRS stratification to contextualize polygenic burden.

That integration is important because complex disease is both polygenic and regulatory, not only driven by single high-effect loci [1,2,3].

## 2) PRS: What It Adds Biologically

### 2.1 Biological Concept

Polygenic risk scores (PRS) model cumulative inherited burden across many common variants. Conceptually:

`PRS = Σ (effect allele dosage_i × GWAS effect size_i)`

This enables a shift from single-variant significance to genome-wide genetic liability estimates, which is essential for common complex traits where no single SNP explains most risk [4,5,6].

### 2.2 Why PRS Matters in Disease Biology

PRS can:

- stratify individuals into relative risk distributions for a trait,
- identify tails of liability distribution that may be relevant for screening or prevention studies,
- support research questions about interaction between inherited burden and molecular phenotypes.

Empirically, high PRS can identify subsets with risk comparable to some monogenic-risk strata for specific diseases [7], while still representing probabilistic risk and not diagnosis [6,8].

### 2.3 Alzheimer Context in This Project

This pipeline uses PGS Catalog-compatible PRS models and supports AD-oriented use cases, including scores derived from large AD GWAS (for example, Bellenguez et al.) [9,10,11].

Biological value in AD context:

- captures distributed polygenic architecture beyond APOE-centered interpretation alone,
- supports cohort-level burden interpretation when genotype coverage is adequate,
- allows transparent reporting of coverage limits when matched model loci are low.

### 2.4 Why PRS QC Is Biologically Critical

Without strict QC, PRS can be biologically misleading. Key checks implemented in the pipeline are therefore biologically necessary:

- genome build consistency (e.g., GRCh38 alignment),
- allele harmonization and strand handling,
- matched-model-variant counts,
- explicit reporting when coverage is insufficient for interpretation.

These controls align with PRS reporting and reproducibility recommendations [8,11,12].

### 2.5 PRS Limitations (Scientific, Not Optional)

Important interpretation boundaries:

- PRS portability across ancestries is still uneven due to Eurocentric GWAS composition [13].
- PRS does not encode environmental exposure, lifestyle, or disease trajectory.
- Low overlap between callable loci and model weights reduces interpretability.

Accordingly, PRS in this pipeline should be used as research-grade stratification evidence unless clinical validation context is explicitly provided.

## 3) Differential Expression: What It Adds Biologically

### 3.1 Biological Concept

Differential expression (DE) quantifies condition-associated transcript abundance changes, capturing dynamic biology not inferable from germline variation alone.

If PRS addresses inherited predisposition, DE addresses active molecular response.

### 3.2 Why DE Is Complementary to PRS/Variants

Variant and GWAS analyses can identify susceptibility regions, but they do not directly show:

- which genes are transcriptionally perturbed in the measured context,
- whether pathway activity is currently shifted,
- whether molecular changes align with variant-prioritized biology.

DE fills that gap by providing expression effect direction and magnitude per gene and by enabling functional interpretation with pathway tools [14,15,16,17].

### 3.3 Statistical Relevance

The pipeline uses established count-based DE methodology (DESeq2 in RNA branch), which is biologically important because it stabilizes inference under:

- low count regimes,
- gene-wise dispersion heterogeneity,
- limited replicate designs.

Shrinkage-based modeling improves effect-size interpretability and controls false positives in realistic RNA-seq settings [14].

### 3.4 Biological Interpretation Caveats

DE associations are context-dependent:

- tissue composition, cell mixture effects, and batch structure can change observed signals,
- DE does not by itself prove causality,
- pathway enrichment should be interpreted with appropriate multiple-testing control (FDR).

This is why reproducible QC + transparent thresholds are as biologically important as the DE test itself [17,18].

## 4) Why Integrated Pipelines Like This One Matter

### 4.1 From Association to Mechanism Hypotheses

A reproducible multi-mode pipeline enables evidence stacking:

1. GWAS/variant branch: where risk-associated loci are.
2. Annotation branch: which genes/features they may affect.
3. DE branch (when RNA is available): whether related genes/pathways are transcriptionally altered in context.
4. Enrichment branch: whether coherent biological systems emerge.
5. PRS branch: whether aggregate inherited burden supports cohort-level stratification.

Convergence across these layers is stronger biologically than any single metric alone.

### 4.2 Translational Utility

For research and translational preclinical pipelines, this architecture supports:

- prioritizing candidate pathways for follow-up experiments,
- selecting genes/loci for functional validation,
- triaging samples/cohorts for deeper phenotyping based on polygenic burden and molecular profiles,
- generating publication-ready, auditable evidence trails.

### 4.3 Reproducibility as Biological Validity Infrastructure

In computational biology, reproducibility is not only engineering hygiene; it is part of biological validity. If the same inputs cannot regenerate the same results, downstream biological claims are weak.

This pipeline’s Snakemake orchestration, environment isolation, preflight checks, and run manifests support FAIR-aligned reproducible interpretation [19,20,21].

## 5) Practical Interpretation Framework for Users

When interpreting outputs from this project:

1. Start with preflight and QC artifacts (`preflight_checks.tsv`, run manifest, variant/RNA QC).
2. Evaluate PRS coverage metrics before interpreting PRS magnitude.
3. Use DE + enrichment as hypothesis generation unless validated externally.
4. Prioritize biology where variant/GWAS evidence and expression/pathway evidence converge.
5. Report ancestry/build/resource versions and threshold settings with every result.

## 6) References

1. Visscher PM, Wray NR, Zhang Q, et al. 10 Years of GWAS Discovery: Biology, Function, and Translation. *Am J Hum Genet* (2017). https://doi.org/10.1016/j.ajhg.2017.06.005  
2. Cerezo M, Buniello A, et al. The NHGRI-EBI GWAS Catalog: standards for reusability, sustainability and diversity. *Nucleic Acids Res* (2025). https://academic.oup.com/nar/article/53/D1/D998/7893318  
3. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* (2020). https://doi.org/10.1126/science.aaz1776  
4. Choi SW, Mak TSH, O'Reilly PF. Tutorial: a guide to performing polygenic risk score analyses. *Nat Protoc* (2020). https://doi.org/10.1038/s41596-020-0353-1  
5. Torkamani A, Wineinger NE, Topol EJ. The personal and clinical utility of polygenic risk scores. *Nat Rev Genet* (2018). https://doi.org/10.1038/s41576-018-0018-x  
6. Lewis CM, Vassos E. Polygenic risk scores: from research tools to clinical instruments. *Genome Med* (2020). https://doi.org/10.1186/s13073-020-00742-5  
7. Khera AV, Chaffin M, Aragam KG, et al. Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. *Nat Genet* (2018). https://doi.org/10.1038/s41588-018-0183-z  
8. Wand H, Lambert SA, Tamburro C, et al. Improving reporting standards for polygenic scores in risk prediction studies. *Nature* (2021). https://doi.org/10.1038/s41586-021-03243-6  
9. Bellenguez C, Küçükali F, Jansen IE, et al. New insights into the genetic etiology of Alzheimer’s disease and related dementias. *Nat Genet* (2022). https://doi.org/10.1038/s41588-022-01024-z  
10. Lambert SA, Gil L, Jupp S, et al. The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation. *Nat Genet* (2021). https://doi.org/10.1038/s41588-021-00783-5  
11. Lambert SA, Abraham G, Inouye M. Towards clinical utility of polygenic risk scores. *Hum Mol Genet* (2019). https://doi.org/10.1093/hmg/ddz187  
12. PGS Catalog. Score PGS002280 (Alzheimer context). https://www.pgscatalog.org/score/PGS002280/  
13. Martin AR, Kanai M, Kamatani Y, et al. Clinical use of current polygenic risk scores may exacerbate health disparities. *Nat Genet* (2019). https://doi.org/10.1038/s41588-019-0379-x  
14. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* (2014). https://doi.org/10.1186/s13059-014-0550-8  
15. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* (2010). https://doi.org/10.1093/bioinformatics/btp616  
16. Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol* (2014). https://doi.org/10.1186/gb-2014-15-2-r29  
17. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* (2005). https://doi.org/10.1073/pnas.0506580102  
18. Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* (2021). https://doi.org/10.1016/j.xinn.2021.100141  
19. Köster J, Rahmann S. Snakemake: a scalable bioinformatics workflow engine. *Bioinformatics* (2012). https://doi.org/10.1093/bioinformatics/bts480  
20. Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. *F1000Research* (2021). https://doi.org/10.12688/f1000research.29032.2  
21. Wilkinson MD, Dumontier M, Aalbersberg IJJ, et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data* (2016). https://doi.org/10.1038/sdata.2016.18  
