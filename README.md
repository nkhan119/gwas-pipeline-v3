# GWAS Pipeline — CDC 1.0.0
**Author:** Nadeem Khan, Bioinformatician, INRS-CAFSB

---

## Directory structure expected in `--raw_dir`

```
{raw_dir}/
  Cohort_A/
    Cohort_A.bed / .bim / .fam     ← genotypes
    phenotypes.tsv                  ← if --use_existing_pheno true
    covariates.tsv                  ← if --use_existing_pheno true
  Cohort_B/
    Cohort_B.bed / .bim / .fam
    phenotypes.tsv
    covariates.tsv
```

## Your specific paths

```
raw_dir = ~/Downloads/1000G_phase3_common_norel/GWAS_analysis
```

Each cohort folder already contains `phenotypes.tsv` and `covariates.tsv`
from `gwas_step0.sh` — the pipeline reads these directly when
`--use_existing_pheno true` (default).

---

## All outputs go to `--out_dir`

Default: `{pipeline_dir}/results/`

```
results/
  GWAS_Report.html              ← interactive report
  qc/Cohort_A/                  ← QC-passed genotypes
  qc/Cohort_B/
  pca/                          ← eigenvec, eigenval, prune.in
  phenotypes/Cohort_A/          ← only if simulating; otherwise reads from raw_dir
  gwas/Cohort_A/sumstats/       ← harmonised sumstats (×7 traits)
  gwas/Cohort_B/sumstats/
  logs/                         
```

---

## Run commands

### Scenario 1 — You already have everything, just want the HTML report
*(Your current situation: QC, PCA, phenotypes, and REGENIE are all done)*

```bash
cd ~/Downloads/gwas-pipeline-v3

nextflow run main.nf \
    -profile local,conda \
    --raw_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --out_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --skip_qc \
    --skip_pca \
    --skip_regenie \
    -resume
```

The report reads sumstats from:
`{out_dir}/gwas/**/sumstats/*_sumstats.tsv.gz`

### Scenario 2 — Run from scratch using your existing phenotypes
*(Full pipeline but use the real phenotypes.tsv already in Cohort_A/ Cohort_B/)*

```bash
nextflow run main.nf \
    -profile local,conda \
    --raw_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --out_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --use_existing_pheno true \
    -resume
```

### Scenario 3 — Run from scratch with simulated phenotypes
*(Regenerates phenotypes from PCA — useful when starting fresh on new data)*

```bash
nextflow run main.nf \
    -profile local,conda \
    --raw_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --out_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --use_existing_pheno false \
    -resume
```

### Scenario 4 — Run on Narval from scratch

```bash
# On Narval
nextflow run main.nf \
    -profile narval,conda \
    --raw_dir ~/projects/def-fveyrier/1000G/GWAS_analysis \
    --out_dir ~/projects/def-fveyrier/1000G/GWAS_results \
    --use_existing_pheno true \
    -resume
```

### Scenario 5 — Custom output directory

```bash
nextflow run main.nf \
    -profile local,conda \
    --raw_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --out_dir ~/Downloads/my_gwas_results \
    --skip_qc --skip_pca --skip_regenie \
    -resume
```

---

## Skip flags explained

| Flag | Skips | Requires in out_dir |
|---|---|---|
| `--skip_qc` | QC + merge + LD prune (A1–A3) | `qc/{cohort}/{cohort}.bed` and `pca/pruned_global.prune.in` |
| `--skip_pca` | PCA (A4) — only works with `--skip_qc` | `pca/1000G_pca.eigenvec` |
| `--skip_regenie` | REGENIE Step1+2 + harmonise (A8–A10) | `gwas/**/sumstats/*_sumstats.tsv.gz` |
| `--skip_report` | HTML report (A11) | — |

---

## `use_existing_pheno` flag

| Value | Behaviour |
|---|---|
| `true` (default) | Reads `phenotypes.tsv` + `covariates.tsv` from `raw_dir/Cohort_*/` |
| `false` | Simulates 4 QT + 3 BT phenotypes from PCA eigenvectors (seed=42) |

---

## Conda env (install once)

```bash
conda env create -f envs/gwas_env.yaml
conda activate cdc_gwas
```
