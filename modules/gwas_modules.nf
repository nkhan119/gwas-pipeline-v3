// ================================================================
//  modules/gwas_modules.nf
// ================================================================

// ── Simulate phenotypes from PCA eigenvectors ──────────────
// Only used when --use_existing_pheno false
process SIMULATE_PHENOTYPES {
    label 'small'
    tag   "${cohort}"
    publishDir "${params.out_dir}/phenotypes/${cohort}", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple val(cohort), path(bed), path(bim), path(fam)
    path  eigenvec

    output:
    tuple val(cohort), path("${cohort}_phenotypes.tsv"), emit: phenotypes

    script:
    """
    python3 << 'PYEOF'
import numpy as np, pandas as pd
from scipy import stats
np.random.seed(42)

cohort = "${cohort}"
fam    = pd.read_csv("${fam}", sep=r"\\s+", header=None,
                     names=["FID","IID","PAT","MAT","SEX","PHENO"])
pca    = pd.read_csv("${eigenvec}", sep=r"\\s+")
if "#IID" in pca.columns: pca=pca.rename(columns={"#IID":"IID"})
elif pca.columns[0].startswith("#"):
    pca.columns = ["FID","IID"] + [f"PC{i}" for i in range(1,len(pca.columns)-1)]
if "FID" not in pca.columns: pca.insert(0,"FID",pca["IID"])

pc_cols = [f"PC{i}" for i in range(1,11)]
merged  = fam[["FID","IID","SEX"]].merge(pca[["IID"]+pc_cols], on="IID", how="left")
n       = len(merged)
pc_mat  = merged[pc_cols].fillna(0).values

def pce(w): e=pc_mat@np.array(w); return (e-e.mean())/(e.std()+1e-8)

QT = {
    "ldl_cholesterol":{"mean":130,"sd":35,"h2":.50,"w":[ 8,-4, 2,-1, 1,0,0,0,0,0]},
    "bmi":            {"mean": 27,"sd": 5,"h2":.40,"w":[ 3, 2,-3, 1,-2,0,0,0,0,0]},
    "crp_log":        {"mean":  1,"sd": 1,"h2":.25,"w":[ 2,-1, 1, 2,-1,0,0,0,0,0]},
    "height_cm":      {"mean":170,"sd": 9,"h2":.80,"w":[ 5, 3,-2, 0, 1,0,0,0,0,0]},
}
BT = {
    "cad":          {"prev":.08, "h2":.40,"w":[ 3,-2, 1,0, 1,0,0,0,0,0]},
    "t2d":          {"prev":.12, "h2":.35,"w":[ 2, 3,-1,1,-2,0,0,0,0,0]},
    "hypertension": {"prev":.30, "h2":.30,"w":[ 4, 2,-2,1, 0,0,0,0,0,0]},
}
B_SHIFT={"Cohort_B":{"ldl_cholesterol":-10,"bmi":1.5,"crp_log":.2,"height_cm":-2}}
B_PREV= {"Cohort_B":{"cad":.02,"t2d":.03,"hypertension":.05}}

df = merged[["FID","IID"]].copy()
for t,p in QT.items():
    g = np.sqrt(p["h2"])*pce(p["w"])*p["sd"]
    e = np.random.normal(0,np.sqrt(1-p["h2"])*p["sd"],n)
    df[t] = np.round(p["mean"]+g+e+B_SHIFT.get(cohort,{}).get(t,0),4)
for t,p in BT.items():
    prev = p["prev"]+B_PREV.get(cohort,{}).get(t,0)
    lia  = np.sqrt(p["h2"])*pce(p["w"])+np.random.normal(0,np.sqrt(1-p["h2"]),n)
    df[t] = np.where(lia>=stats.norm.ppf(1-prev),2,1).astype(int)
age = 45+10*pce([1,0,0,0,0,0,0,0,0,0])+np.random.normal(0,8,n)
df["Age"] = np.clip(np.round(age).astype(int),20,80)
df["Sex"] = merged["SEX"].values
df.to_csv("${cohort}_phenotypes.tsv", sep="\\t", index=False)
print(f"[PHENO] Simulated {n} samples for {cohort}")
for t in list(QT): print(f"  {t:<22} mean={df[t].mean():.2f}  SD={df[t].std():.2f}")
for t in list(BT):
    nc=(df[t]==2).sum(); print(f"  {t:<22} prev={nc/n:.1%}  n_cases={nc:,}")
PYEOF
    """
}

// ─ Build covariates (PC1-10 + Age + Sex, z-scored) ────────
process BUILD_COVARIATES {
    label 'small'
    tag   "${cohort}"
    publishDir "${params.out_dir}/phenotypes/${cohort}", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple val(cohort), path(phenotypes)
    path  eigenvec

    output:
    tuple val(cohort), path("${cohort}_covariates.tsv"), emit: covariates

    script:
    """
    python3 << 'PYEOF'
import pandas as pd, numpy as np
pca = pd.read_csv("${eigenvec}", sep=r"\\s+")
if "#IID" in pca.columns: pca=pca.rename(columns={"#IID":"IID"})
elif pca.columns[0].startswith("#"):
    pca.columns=["FID","IID"]+[f"PC{i}" for i in range(1,len(pca.columns)-1)]
if "FID" not in pca.columns: pca.insert(0,"FID",pca["IID"])
pc_cols=[f"PC{i}" for i in range(1,11)]
pheno=pd.read_csv("${phenotypes}",sep="\\t")
merged=pheno[["FID","IID","Age","Sex"]].merge(pca[["IID"]+pc_cols],on="IID",how="left")
for c in pc_cols+["Age"]:
    merged[c]=(merged[c]-merged[c].mean())/(merged[c].std()+1e-8)
merged[["FID","IID","Age","Sex"]+pc_cols].to_csv(
    "${cohort}_covariates.tsv",sep="\\t",index=False,float_format="%.6f")
print(f"[COV] {cohort}: {len(merged)} samples  Age+Sex+{len(pc_cols)} PCs")
PYEOF
    """
}

// ── Recode binary phenotypes PLINK 1/2 → REGENIE 0/1/NA ───
process RECODE_BINARY {
    label 'tiny'
    tag   "${cohort}"
    publishDir "${params.out_dir}/phenotypes/${cohort}", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple val(cohort), path(phenotypes)

    output:
    tuple val(cohort), path("${cohort}_phenotypes_bt.tsv"), emit: bt_phenotypes

    script:
    """
    python3 << 'PYEOF'
import pandas as pd, sys
bt_cols = [c.strip() for c in "${params.binary_traits}".split(",")]
df      = pd.read_csv("${phenotypes}", sep="\\t")
missing = [c for c in bt_cols if c not in df.columns]
if missing:
    print(f"ERROR: BT columns not found: {missing}\\nAvailable: {df.columns.tolist()}", file=sys.stderr)
    sys.exit(1)
out = df[["FID","IID"]+bt_cols].copy()
for col in bt_cols:
    s = pd.to_numeric(out[col], errors="coerce").replace({-9.0: float("nan")})
    ctrl=(s==1).sum(); case=(s==2).sum()
    out[col] = s.map({1.0:0, 2.0:1})
    print(f"  {col:<22}  controls={ctrl:>5,}  cases={case:>5,}")
out.to_csv("${cohort}_phenotypes_bt.tsv", sep="\\t", index=False, na_rep="NA")
print(f"[RECODE] ${cohort} → ${cohort}_phenotypes_bt.tsv")
PYEOF
    """
}

// ── REGENIE Step 1 (QT and BT separately) ──────────────────
process REGENIE_STEP1 {
    label 'xlarge'
    tag   "${cohort}_${trait_type}"
    publishDir "${params.out_dir}/gwas/${cohort}/step1", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple val(cohort),
          path(bed), path(bim), path(fam),
          path(phenotypes), path(covariates),
          path(prune_in)
    val   trait_type   // "qt" or "bt"

    output:
    tuple val(cohort), val(trait_type),
          path("step1_${trait_type}_pred.list"),
          path("step1_${trait_type}*.loco.gz"),
          emit: step1_pred

    script:
    def traits = trait_type == "qt" ? params.quant_traits : params.binary_traits
    def flag   = trait_type == "qt" ? "--qt" : "--bt"
    """
    set -euo pipefail
    STEM=\$(basename ${bed} .bed)
    echo "[REGENIE S1] ${cohort} ${trait_type}: ${traits}"

    # Extract pruned SNPs for Step 1
    plink2 --bfile "\${STEM}" --extract ${prune_in} \\
           --make-bed --out step1_input --threads ${params.threads}
    echo "[REGENIE S1] Pruned SNPs for Step1: \$(wc -l < step1_input.bim)"

    regenie --step 1 \\
        --bed          step1_input \\
        --phenoFile    ${phenotypes} \\
        --covarFile    ${covariates} \\
        --phenoColList ${traits} \\
        ${flag} \\
        --bsize        ${params.step1_bsize} \\
        --lowmem \\
        --lowmem-prefix tmp_${trait_type} \\
        --out          step1_${trait_type} \\
        --threads      ${params.threads} \\
        --gz

    echo "[REGENIE S1] ${cohort} ${trait_type} complete"
    """
}

// ── A9: REGENIE Step 2 (per-variant association) ───────────────
process REGENIE_STEP2 {
    label 'xlarge'
    tag   "${cohort}_${trait_type_val}"
    publishDir "${params.out_dir}/gwas/${cohort}/step2", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple val(cohort),
          path(bed), path(bim), path(fam),
          path(phenotypes), path(covariates),
          path(prune_in),
          val(trait_type), path(pred_list), path(loco_files)
    val   trait_type_val

    output:
    tuple val(cohort), val(trait_type_val),
          path("${cohort}_${trait_type_val}_*.regenie.gz"),
          emit: step2_out

    script:
    def traits = trait_type_val == "qt" ? params.quant_traits : params.binary_traits
    def flag   = trait_type_val == "qt" ? "--qt" : "--bt"
    def firth  = trait_type_val == "bt"
                 ? "--firth --approx --pThresh ${params.firth_pthresh}" : ""
    """
    STEM=\$(basename ${bed} .bed)
    echo "[REGENIE S2] ${cohort} ${trait_type_val}: ${traits}"
    regenie --step 2 --bed "\${STEM}" \\
        --phenoFile    ${phenotypes} \\
        --covarFile    ${covariates} \\
        --phenoColList ${traits} \\
        ${flag} \\
        --pred  ${pred_list} \\
        ${firth} \\
        --bsize ${params.step2_bsize} \\
        --out   ${cohort}_${trait_type_val} \\
        --threads ${params.threads} \\
        --gz
    echo "[REGENIE S2] ${cohort} ${trait_type_val} complete"
    """
}

// ── Harmonise to MR-ready sumstats ────────────────────────
// REGENIE: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P
// Output:  CHR BP SNP A1 A2 BETA SE P LOG10P N EAF TRAIT COHORT
process HARMONISE_SUMSTATS {
    label 'small'
    tag   "${cohort}_${trait_type}"
    publishDir "${params.out_dir}/gwas/${cohort}/sumstats", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple val(cohort), val(trait_type), path(regenie_files)

    output:
    tuple val(cohort), path("${cohort}_*_sumstats.tsv.gz"), emit: sumstats
    path "${cohort}_${trait_type}_summary.tsv",             emit: gwas_summary

    script:
    def traits = (trait_type == "qt" ? params.quant_traits : params.binary_traits)
                  .split(",").collect{it.trim()}.join(",")
    """
    python3 << 'PYEOF'
import pandas as pd, numpy as np, glob
from pathlib import Path

cohort     = "${cohort}"
trait_type = "${trait_type}"
traits     = [t.strip() for t in "${traits}".split(",")]

COL_MAP = {
    "CHROM":"CHR", "GENPOS":"BP", "ID":"SNP",
    "ALLELE0":"A2",   # REGENIE: non-effect allele
    "ALLELE1":"A1",   # REGENIE: effect allele
    "A1FREQ":"EAF", "N":"N", "BETA":"BETA", "SE":"SE", "LOG10P":"LOG10P"
}

summary = []
for trait in traits:
    df = None
    for ext in [f"${cohort}_{trait_type}_{trait}.regenie.gz",
                f"${cohort}_{trait_type}_{trait}.regenie"]:
        p = Path(ext)
        if p.exists():
            df = pd.read_csv(p, sep=r"\\s+", dtype=str,
                             compression="gzip" if ext.endswith(".gz") else None)
            break
    if df is None:
        print(f"  [WARN] Not found: ${cohort}_{trait_type}_{trait}.regenie(.gz)")
        continue

    df = df.rename(columns=COL_MAP)
    keep = [c for c in ["CHR","BP","SNP","A1","A2","BETA","SE","LOG10P","EAF","N"]
            if c in df.columns]
    df = df[keep]
    for c in ["BP","BETA","SE","LOG10P","EAF"]:
        if c in df: df[c] = pd.to_numeric(df[c], errors="coerce")
    if "N"   in df: df["N"]   = pd.to_numeric(df["N"],   errors="coerce").astype("Int64")
    if "CHR" in df: df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce").astype("Int64")
    df["P"]      = 10 ** (-df["LOG10P"])
    df["TRAIT"]  = trait
    df["COHORT"] = cohort

    n0 = len(df)
    df = df.dropna(subset=["BETA","SE","P"])
    n_drop = n0 - len(df)

    out = f"${cohort}_{trait}_sumstats.tsv.gz"
    df.to_csv(out, sep="\\t", index=False, compression="gzip")

    z2  = (df["BETA"]/df["SE"])**2
    lgc = round(float(np.median(z2)/0.4549), 4) if len(df) else None
    n_sig = int((df["P"]<5e-8).sum())
    summary.append({
        "cohort":cohort,"trait":trait,"trait_type":trait_type,
        "n_variants":len(df),"n_dropped":n_drop,
        "n_gw_sig":n_sig,"min_p":float(df["P"].min()) if len(df) else None,
        "lambda_gc":lgc
    })
    print(f"  ✓ {trait:<22} {len(df):>9,} variants  GW-sig={n_sig:>5,}  λ={lgc}  min_P={df['P'].min():.2e}")

pd.DataFrame(summary).to_csv(
    f"${cohort}_{trait_type}_summary.tsv", sep="\\t", index=False)
print(f"[HARMONISE] Summary → ${cohort}_{trait_type}_summary.tsv")
PYEOF
    """
}
