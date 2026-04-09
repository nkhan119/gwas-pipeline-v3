// ================================================================
//  modules/merge_ld_pca.nf
// ================================================================

process MERGE_COHORTS {
    label 'large'
    publishDir "${params.out_dir}/pca", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    path qc_files   // collected: Cohort_A.bed/.bim/.fam  Cohort_B.bed/.bim/.fam

    output:
    tuple path("1000G_dedup.bed"), path("1000G_dedup.bim"), path("1000G_dedup.fam"),
          emit: dedup_plink

    script:
    """
    set -euo pipefail

    # ── Identify the two cohort stems ────────────────────────────
    ls *.bim | sed 's/\\.bim//' | sort -u > stems.txt
    mapfile -t STEMS < stems.txt
    A="\${STEMS[0]}"; B="\${STEMS[1]}"
    echo "[MERGE] Cohorts: \${A}  \${B}"

    # ── Extract common SNPs ──────────────────────────────────────
    awk '{print \$2}' "\${A}.bim" | sort > snps_a.txt
    awk '{print \$2}' "\${B}.bim" | sort > snps_b.txt
    comm -12 snps_a.txt snps_b.txt > common_snps.txt
    N_COMMON=\$(wc -l < common_snps.txt)
    echo "[MERGE] Common SNPs: \${N_COMMON}"

    if [ "\${N_COMMON}" -eq 0 ]; then
        echo "[ERROR] No common SNPs between cohorts — check variant IDs"; exit 1
    fi

    # ── Filter each cohort to common SNPs ───────────────────────
    for S in "\${A}" "\${B}"; do
        plink2 --bfile "\${S}" --extract common_snps.txt \\
               --make-bed --out "\${S}_common" --threads ${params.threads}
    done

    # ── Merge with plink1 (robust, handles strand flips) ────────
    # plink2 --pmerge-list only scans files in the first pass and does
    # not write .bed — plink1 --bmerge is the reliable choice here.
    plink --bfile "\${A}_common" \\
          --bmerge "\${B}_common" \\
          --allow-no-sex \\
          --make-bed --out merged \\
          --threads ${params.threads} 2>&1 | tee merge.log || true

    # ── If mis-SNP file produced, remove problem SNPs and retry ─
    if [ -f merged-merge.missnp ]; then
        N_MIS=\$(wc -l < merged-merge.missnp)
        echo "[MERGE] \${N_MIS} mis-SNPs detected — excluding and retrying..."
        for S in "\${A}_common" "\${B}_common"; do
            plink --bfile "\${S}" --exclude merged-merge.missnp \\
                  --make-bed --out "\${S}_clean" --allow-no-sex \\
                  --threads ${params.threads}
        done
        plink --bfile "\${A}_common_clean" \\
              --bmerge "\${B}_common_clean" \\
              --allow-no-sex \\
              --make-bed --out merged \\
              --threads ${params.threads}
    fi

    if [ ! -f merged.bed ]; then
        echo "[ERROR] Merge failed — check merge.log"
        cat merge.log
        exit 1
    fi
    echo "[MERGE] Merged: \$(wc -l < merged.fam) samples  \$(wc -l < merged.bim) variants"

    # ── Deduplicate positions ────────────────────────────────────
    plink2 --bfile merged --rm-dup exclude-all \\
           --make-bed --out 1000G_dedup --threads ${params.threads}
    echo "[MERGE] Final: \$(wc -l < 1000G_dedup.fam) samples  \$(wc -l < 1000G_dedup.bim) variants"
    """
}

// ================================================================
process LD_PRUNE {
    label 'medium'
    publishDir "${params.out_dir}/pca", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "pruned_global.prune.in",  emit: prune_in
    path "pruned_global.prune.out", emit: prune_out

    script:
    """
    STEM=\$(basename ${bed} .bed)
    echo "[PRUNE] win=${params.prune_win}kb  step=${params.prune_step}  r2<${params.prune_r2}"
    plink2 --bfile "\${STEM}" \\
           --indep-pairwise ${params.prune_win} ${params.prune_step} ${params.prune_r2} \\
           --out pruned_global --threads ${params.threads}
    echo "[PRUNE] In: \$(wc -l < pruned_global.prune.in)  Out: \$(wc -l < pruned_global.prune.out)"
    """
}

// ================================================================
process RUN_PCA {
    label 'large'
    publishDir "${params.out_dir}/pca", mode: 'copy'
    conda "${projectDir}/envs/gwas_env.yaml"

    input:
    tuple path(bed), path(bim), path(fam)
    path  prune_in

    output:
    path "1000G_pca.eigenvec",  emit: eigenvec
    path "1000G_pca.eigenval",  emit: eigenval
    path "pca_variance.tsv",    emit: pca_variance

    script:
    """
    STEM=\$(basename ${bed} .bed)
    echo "[PCA] Running ${params.n_pcs} PCs on \$(wc -l < ${prune_in}) pruned variants"
    plink2 --bfile "\${STEM}" --extract ${prune_in} \\
           --pca approx ${params.n_pcs} --out 1000G_pca \\
           --threads ${params.threads}

    python3 << 'PYEOF'
import csv
vals = [float(l) for l in open("1000G_pca.eigenval").read().splitlines() if l.strip()]
tot  = sum(vals); cum = 0
with open("pca_variance.tsv","w",newline="") as f:
    w = csv.writer(f, delimiter="\\t")
    w.writerow(["PC","eigenvalue","pct_variance","cumulative_pct"])
    for i,v in enumerate(vals, 1):
        p = 100*v/tot; cum += p
        w.writerow([f"PC{i}", f"{v:.4f}", f"{p:.4f}", f"{cum:.4f}"])
        if i<=10: print(f"  PC{i:>2}: {p:.2f}%")
print(f"  PC1-10 cumulative: {sum(100*v/tot for v in vals[:10]):.1f}%")
PYEOF
    """
}

// ================================================================
//  Sub-workflow — called from main.nf as MERGE_LD_PCA(collected_beds)
//  Emits: eigenvec, eigenval, pca_variance, prune_in
// ================================================================
workflow MERGE_LD_PCA {
    take:
    qc_files

    main:
    MERGE_COHORTS(qc_files)
    LD_PRUNE(MERGE_COHORTS.out.dedup_plink)
    RUN_PCA(MERGE_COHORTS.out.dedup_plink, LD_PRUNE.out.prune_in)

    emit:
    eigenvec     = RUN_PCA.out.eigenvec
    eigenval     = RUN_PCA.out.eigenval
    pca_variance = RUN_PCA.out.pca_variance
    prune_in     = LD_PRUNE.out.prune_in
}
