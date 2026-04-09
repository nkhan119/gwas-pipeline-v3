// ================================================================
//  modules/gwas_report.nf
//  Interactive Plotly HTML report
// ================================================================
process GWAS_REPORT {
    label 'small'
    publishDir "${params.out_dir}", mode: 'copy'
    conda      "${projectDir}/envs/gwas_env.yaml"

    input:
    path eigenvec
    path eigenval
    path pca_variance
    path qc_metrics_files
    path sumstat_files
    path gwas_summary_files
    val  author
    val  affiliation

    output:
    path "GWAS_Report.html", emit: report

    script:
    def institute = params.institute      ?: "INRS-Centre Armand-Frappier Sante-Biotechnologie"
    def github    = params.github         ?: "github.com/nkhan119"
    def pop_file  = params.population_file?: ""
    """
    python3 ${projectDir}/bin/render_report.py \
        --eigenvec      ${eigenvec} \
        --pca_variance  ${pca_variance} \
        --template      ${projectDir}/assets/gwas_report_template.html \
        --author        "${author}" \
        --affiliation   "${affiliation}" \
        --institute     "${institute}" \
        --github        "${github}" \
        --population_file "${pop_file}" \
        --maf      ${params.maf} \
        --hwe      ${params.hwe} \
        --geno     ${params.geno} \
        --mind     ${params.mind} \
        --win      ${params.prune_win} \
        --r2       ${params.prune_r2} \
        --npcs     ${params.n_pcs} \
        --s1bsize  ${params.step1_bsize} \
        --s2bsize  ${params.step2_bsize}
    """
}
