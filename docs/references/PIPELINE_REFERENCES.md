# Pipeline References

This folder lists the primary papers and official resources used to design and implement this pipeline.

## Workflow and Reproducibility

1. Snakemake (original)
- Koester J, Rahmann S. Snakemake - a scalable bioinformatics workflow engine. Bioinformatics (2012).
- https://doi.org/10.1093/bioinformatics/bts480

2. Snakemake (current architecture)
- Moelder F et al. Sustainable data analysis with Snakemake. F1000Research (2021).
- https://doi.org/10.12688/f1000research.29032.2

## QC and Preprocessing

3. FastQC
- Official project page: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

4. fastp
- Chen S et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics (2018).
- https://doi.org/10.1093/bioinformatics/bty560

5. MultiQC
- Ewels P et al. MultiQC: summarize analysis results for multiple tools. Bioinformatics (2016).
- https://doi.org/10.1093/bioinformatics/btw354

## Alignment and Variant Processing

6. BWA-MEM2
- Official repository (used for implementation details): https://github.com/bwa-mem2/bwa-mem2

7. STAR
- Dobin A et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics (2013).
- https://doi.org/10.1093/bioinformatics/bts635

8. Minimap2
- Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics (2018).
- https://doi.org/10.1093/bioinformatics/bty191

9. SAMtools / BCFtools
- Li H et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics (2009).
- https://doi.org/10.1093/bioinformatics/btp352
- Danecek P et al. Twelve years of SAMtools and BCFtools. GigaScience (2021).
- https://doi.org/10.1093/gigascience/giab008

10. BEDTools (variant-feature overlap annotation)
- Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics (2010).
- https://doi.org/10.1093/bioinformatics/btq033

## RNA Quantification and Differential Expression

11. StringTie
- Pertea M et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol (2015).
- https://doi.org/10.1038/nbt.3122

12. featureCounts (Subread)
- Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads. Bioinformatics (2014).
- https://doi.org/10.1093/bioinformatics/btt656

13. DESeq2
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol (2014).
- https://doi.org/10.1186/s13059-014-0550-8

14. clusterProfiler
- Wu T et al. clusterProfiler 4.0: A universal enrichment tool for omics data. Innovation (2021).
- https://doi.org/10.1016/j.xinn.2021.100141

## Statistical Thresholding References

15. False discovery rate (BH procedure)
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. JRSS B (1995).
- https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

16. GWAS genome-wide significance threshold
- Jannot AS et al. P < 5 x 10^-8 has emerged as a standard of statistical significance for genome-wide association studies. J Clin Epidemiol (2015).
- https://doi.org/10.1016/j.jclinepi.2015.01.001
- Cheruiyot EK et al. GWAS significance thresholds in large cohorts of European ancestry. Genetics (2025).
- https://doi.org/10.1093/genetics/iyaf056

## PRS and Alzheimer Model

17. PGS Catalog (resource)
- Lambert SA et al. The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation. Nucleic Acids Res (2021).
- https://doi.org/10.1093/nar/gkaa886
- https://www.pgscatalog.org/

18. Alzheimer PRS model used in this pipeline
- PGS ID: PGS002280
- Score page: https://www.pgscatalog.org/score/PGS002280
- Scoring file: https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002280/ScoringFiles/PGS002280.txt.gz
- Source GWAS paper: Bellenguez C et al. Nat Genet (2022).
- https://doi.org/10.1038/s41588-022-01024-z

## Alzheimer Variant Source Used for Real Disease-Annotated VCF

19. ClinVar
- Landrum MJ et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res (2018).
- https://doi.org/10.1093/nar/gkx1153
- Download used here (GRCh38 VCF):
  - https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
