# Running in Different Environments (Reproducible)

## 1. Local Workstation (Conda/Mamba)

```bash
snakemake --use-conda --cores 8
```

Recommended:

```bash
snakemake --use-conda --conda-frontend mamba --cores 8 --rerun-incomplete --printshellcmds
```

Optional UX wrapper with spinner/loading indicator during quiet intervals:

```bash
./run_pipeline.sh --use-conda --cores 8 --conda-frontend mamba --rerun-incomplete --printshellcmds
```

Preflight-only validation before full execution:

```bash
snakemake --use-conda --conda-frontend mamba --cores 1 --until preflight_resources
```

## 2. Shared Server / HPC

Use Snakemake profiles and run jobs with a scheduler (SLURM/SGE/PBS).

Example (SLURM profile already configured in your environment):

```bash
snakemake --profile slurm --use-conda --conda-frontend mamba
```

## 3. Containerized Execution

If you maintain container images for each step, run with:

```bash
snakemake --use-singularity --cores 8
```

or

```bash
snakemake --software-deployment-method conda apptainer --cores 8
```

## 4. Cloud / Portable Runs

- Keep `config/config.yaml`, `config/samples.tsv`, and `config/metadata.tsv` under version control.
- Pin references and annotation versions (GRCh38 + matching GTF).
- Store a copy of `results/docs/run_manifest.txt` for every run.

## 5. Reproducibility Checklist

1. Same genome build across all resources.
2. Same pipeline commit hash.
3. Same config and metadata files.
4. Same PRS weights version.
5. Keep generated result manifest and MultiQC report.
6. Keep generated preflight validation report (`results/docs/preflight_checks.tsv`).
