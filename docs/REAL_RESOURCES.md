# Real Resource Files (Downloaded)

This project now includes real public resources for enrichment, annotation, and AD GWAS context.

## Enrichment Gene Sets

Local directory: `resources/enrichment/real/`

- `GO_Biological_Process_2023.gmt`
  - Source: `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023`
- `KEGG_2021_Human.gmt`
  - Source: `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human`
- `Reactome_2022.gmt`
  - Source: `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022`
- `Human_AllPathways_noPFOCR_February_03_2026_symbol.gmt`
  - Source: `https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/`

Default config now points to:
- `resources/enrichment/real/Reactome_2022.gmt`

## ClinVar Annotation Resources

Local directory: `resources/annotation/clinvar/`

- `clinvar.vcf.gz` and `clinvar.vcf.gz.tbi`
  - Source: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`
- `clinvar_alzheimers_sites.vcf.gz` and index
  - Built from ClinVar by filtering disease names containing `Alzheimer`.

## Alzheimer GWAS Catalog Resources (Public)

Local directory: `resources/gwas/alzheimers/`

- `gwas_catalog_alzheimers_associations_2026-02-17.tsv`
  - Built from: `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip`
  - Filter: rows containing `alzheimer` (case-insensitive)
- `gwas_catalog_alzheimers_studies_2026-02-17.tsv`
  - Built from: `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-studies.tsv`
  - Filter: rows containing `alzheimer` (case-insensitive)

## Notes

- `resources/reference/star_grch38_index/` is currently empty and still needs STAR index generation for RNA short-read alignment.
- snpEff database files are still optional and only needed if `annotation.method: snpeff` is enabled.
