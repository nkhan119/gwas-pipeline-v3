// ================================================================
//  modules/genotype_qc.nf
// ================================================================
process GENOTYPE_QC {
    label 'large'
    tag   "${cohort}"
    publishDir "${params.out_dir}/qc/${cohort}", mode: 'copy',
               pattern: "${cohort}.{bed,bim,fam}"
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple val(cohort), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort),
          path("${cohort}.bed"), path("${cohort}.bim"), path("${cohort}.fam"),
          emit: qc_beds
    path "${cohort}_qc_metrics.tsv", emit: qc_metrics
    path "${cohort}_qc.log",         emit: qc_log

    script:
    """
    set -euo pipefail
    NS=\$(wc -l < ${fam}); NV=\$(wc -l < ${bim})
    echo "[QC] ${cohort}: \${NS} samples  \${NV} variants" | tee ${cohort}_qc.log

    # A1a variant filters
    plink2 --bed ${bed} --bim ${bim} --fam ${fam} \\
        --maf ${params.maf} --hwe ${params.hwe} midp \\
        --geno ${params.geno} \\
        --make-bed --out ${cohort}_qc1 --threads ${params.threads} \\
        2>&1 | tee -a ${cohort}_qc.log
    echo "[QC] After variant QC: \$(wc -l < ${cohort}_qc1.bim) variants" | tee -a ${cohort}_qc.log

    # A1b sample missingness
    plink2 --bfile ${cohort}_qc1 --mind ${params.mind} \\
        --make-bed --out ${cohort}_qc2 --threads ${params.threads} \\
        2>&1 | tee -a ${cohort}_qc.log
    echo "[QC] After sample QC: \$(wc -l < ${cohort}_qc2.fam) samples" | tee -a ${cohort}_qc.log

    # A1c heterozygosity outliers
    plink2 --bfile ${cohort}_qc2 --het --out ${cohort}_het \\
        --threads ${params.threads} 2>/dev/null || true

    python3 << 'PYEOF'
import pandas as pd, numpy as np
from pathlib import Path
het_f = Path("${cohort}_het.het")
out_f = Path("${cohort}_het_outliers.txt")
if not het_f.exists():
    out_f.write_text(""); print("[QC] No .het file — skipping het filter")
else:
    df  = pd.read_csv(het_f, sep=r"\\s+")
    col = next((c for c in ["F","HET_RATE"] if c in df.columns), None)
    if col is None:
        out_f.write_text("")
        print(f"[QC] Unexpected het columns: {df.columns.tolist()}")
    else:
        mu, sd = df[col].mean(), df[col].std()
        lo, hi = mu - ${params.het_sd}*sd, mu + ${params.het_sd}*sd
        out    = df[(df[col]<lo)|(df[col]>hi)]
        iid    = "#IID" if "#IID" in df.columns else df.columns[0]
        out[[iid,iid]].to_csv(out_f, sep="\\t", index=False, header=False)
        print(f"[QC] Het: mean={mu:.4f} SD={sd:.4f} [{lo:.4f},{hi:.4f}]  outliers={len(out)}")
PYEOF

    plink2 --bfile ${cohort}_qc2 --remove ${cohort}_het_outliers.txt \\
        --make-bed --out ${cohort} --threads ${params.threads} \\
        2>&1 | tee -a ${cohort}_qc.log

    NFS=\$(wc -l < ${cohort}.fam); NFV=\$(wc -l < ${cohort}.bim)
    echo "[QC] Final: \${NFS} samples  \${NFV} variants" | tee -a ${cohort}_qc.log

    python3 << 'PYEOF'
import csv, subprocess
def wc(p):
    try: return int(subprocess.check_output(["wc","-l",p]).split()[0])
    except: return "N/A"
with open("${cohort}_qc_metrics.tsv","w",newline="") as f:
    w = csv.writer(f, delimiter="\\t")
    w.writerow(["metric","value","cohort"])
    for m,v in [
        ("raw_samples",   wc("${fam}")),
        ("raw_variants",  wc("${bim}")),
        ("final_samples", wc("${cohort}.fam")),
        ("final_variants",wc("${cohort}.bim")),
        ("maf",           "${params.maf}"),
        ("hwe",           "${params.hwe}"),
        ("geno",          "${params.geno}"),
        ("mind",          "${params.mind}"),
        ("het_sd",        "${params.het_sd}"),
    ]:
        w.writerow([m, v, "${cohort}"])
PYEOF
    echo "[QC] ${cohort} complete" | tee -a ${cohort}_qc.log
    """
}
