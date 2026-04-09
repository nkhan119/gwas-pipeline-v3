// ================================================================
//  modules/generate_figures.nf
// ================================================================
process GENERATE_FIGURES {
    label 'medium'
    publishDir "${params.out_dir}/figures", mode: 'copy'
    conda      "${projectDir}/envs/gwas_env.yaml"

    input:
    path  eigenvec
    path  eigenval
    path  sumstat_files        // collected gz files
    path  summary_files        // collected summary TSVs
    val   author
    val   institute
    val   affiliation

    output:
    path "F1_miami_plot.png",       emit: miami_png,    optional: true
    path "F1_miami_plot.pdf",       emit: miami_pdf,    optional: true
    path "F2_qq_plots.png",         emit: qq_png,       optional: true
    path "F2_qq_plots.pdf",         emit: qq_pdf,       optional: true
    path "F3a_pca_combined.png",    emit: pca_all_png,  optional: true
    path "F3a_pca_combined.pdf",    emit: pca_all_pdf,  optional: true
    path "F3b_pca_cohort_A.png",    emit: pca_a_png,    optional: true
    path "F3b_pca_cohort_A.pdf",    emit: pca_a_pdf,    optional: true
    path "F3c_pca_cohort_B.png",    emit: pca_b_png,    optional: true
    path "F3c_pca_cohort_B.pdf",    emit: pca_b_pdf,    optional: true
    path "F4_lgc_summary.png",      emit: lgc_png,      optional: true
    path "F4_lgc_summary.pdf",      emit: lgc_pdf,      optional: true
    path "F5_locus_zoom.png",       emit: locus_png,    optional: true
    path "F5_locus_zoom.pdf",       emit: locus_pdf,    optional: true
    path "F6_forest_plot.png",      emit: forest_png,   optional: true
    path "F6_forest_plot.pdf",      emit: forest_pdf,   optional: true

    script:
    def pop_arg  = params.population_file
                   ? "--population_file ${params.population_file}" : ""
    def ref_arg  = params.ref_bfile
                   ? "--ref_bfile ${params.ref_bfile}" : ""

    // Per-cohort covariate files for PCA F3b / F3c
    def cov_a = "${params.pheno_dir}/Cohort_A/covariates.tsv"
    def cov_b = "${params.pheno_dir}/Cohort_B/covariates.tsv"
    """
    set -euo pipefail

    # ── Stage sumstats into expected directory layout ────────────
    mkdir -p sumstats_staged
    for gz in *_sumstats.tsv.gz; do
        cohort=\$(echo "\$gz" | grep -oP '^Cohort_[AB]')
        mkdir -p "sumstats_staged/\${cohort}/sumstats"
        ln -sf "\$(pwd)/\${gz}" "sumstats_staged/\${cohort}/sumstats/\${gz}" 2>/dev/null || \
        cp "\${gz}" "sumstats_staged/\${cohort}/sumstats/\${gz}"
    done

    # ── Collect summary TSVs ─────────────────────────────────────
    SUMMARY_ARGS=\$(ls *_summary.tsv 2>/dev/null | tr '\\n' ' ' || echo "")
    [ -z "\$SUMMARY_ARGS" ] && { touch _dummy_summary.tsv; SUMMARY_ARGS="_dummy_summary.tsv"; }

    # ── Run figure generator ─────────────────────────────────────
    python3 ${projectDir}/bin/generate_figures.py \\
        --sumstat_dir     sumstats_staged \\
        --eigenvec        ${eigenvec} \\
        --eigenval        ${eigenval} \\
        --cohort_a_covar  "${cov_a}" \\
        --cohort_b_covar  "${cov_b}" \\
        --summary_tsvs    \$SUMMARY_ARGS \\
        --out_dir         . \\
        --cohorts         "${params.cohorts}" \\
        --quant_traits    "${params.quant_traits}" \\
        --binary_traits   "${params.binary_traits}" \\
        --author          "${author}" \\
        --institute       "${institute}" \\
        ${pop_arg} \\
        ${ref_arg}

    echo "[FIGURES] Complete:"
    ls -lh F*.png F*.pdf 2>/dev/null || echo "  No figure files produced"
    """
}
