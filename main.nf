#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GENOTYPE_QC                       } from './modules/genotype_qc.nf'
include { MERGE_LD_PCA                      } from './modules/merge_ld_pca.nf'
include { RECODE_BINARY                     } from './modules/gwas_modules.nf'
include { REGENIE_STEP1 as REGENIE_STEP1_QT } from './modules/gwas_modules.nf'
include { REGENIE_STEP1 as REGENIE_STEP1_BT } from './modules/gwas_modules.nf'
include { REGENIE_STEP2 as REGENIE_STEP2_QT } from './modules/gwas_modules.nf'
include { REGENIE_STEP2 as REGENIE_STEP2_BT } from './modules/gwas_modules.nf'
include { HARMONISE_SUMSTATS                } from './modules/gwas_modules.nf'
include { GWAS_REPORT                       } from './modules/gwas_report.nf'
include { GENERATE_FIGURES                  } from './modules/generate_figures.nf'

new File("${params.out_dir}/logs").mkdirs()
new File("${params.out_dir}/figures").mkdirs()

log.info """
╔══════════════════════════════════════════════════════════════════╗
║   GWAS Pipeline · CDC 1.0.0 · Phase A                          ║
║   Author    : ${params.author.padRight(52)}║
║   Institute : ${(params.institute ?: "").take(52).padRight(52)}║
╠══════════════════════════════════════════════════════════════════╣
║   raw_dir   : ${params.raw_dir.padRight(52)}║
║   pca_dir   : ${(params.pca_dir ?: params.out_dir+"/pca").padRight(52)}║
║   pheno_dir : ${params.pheno_dir.padRight(52)}║
║   out_dir   : ${params.out_dir.padRight(52)}║
╠══════════════════════════════════════════════════════════════════╣
║   skip_qc       = ${params.skip_qc.toString().padRight(47)}║
║   skip_regenie  = ${params.skip_regenie.toString().padRight(47)}║
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()

workflow {

    def cohort_list = params.cohorts instanceof List
        ? params.cohorts
        : params.cohorts.toString().split(",").collect { it.trim() }

    // ── A1  QC + merge + prune ────────────────────────────────────
    if (!params.skip_qc) {
        ch_raw = Channel.fromList(cohort_list).map { cohort ->
            def base = "${params.raw_dir}/${cohort}/${cohort}"
            tuple(cohort, file("${base}.bed"), file("${base}.bim"), file("${base}.fam"))
        }
        GENOTYPE_QC(ch_raw)
        MERGE_LD_PCA(
            GENOTYPE_QC.out.qc_beds
                .map { cohort, bed, bim, fam -> [bed, bim, fam] }
                .flatten().collect()
        )
        ch_qc_beds  = GENOTYPE_QC.out.qc_beds
        ch_eigenvec = MERGE_LD_PCA.out.eigenvec
        ch_eigenval = MERGE_LD_PCA.out.eigenval
        ch_pca_var  = MERGE_LD_PCA.out.pca_variance
        ch_prune_in = MERGE_LD_PCA.out.prune_in
        ch_qc_files = GENOTYPE_QC.out.qc_metrics.collect()

    } else {
        log.info "⏭  skip_qc — loading existing QC beds and combined PCA from disk"

        // Use raw_dir beds directly — QC was done externally (e.g. gwas_step0.sh)
        ch_qc_beds = Channel.fromList(cohort_list).map { cohort ->
            def bed = file("${params.raw_dir}/${cohort}/${cohort}.bed")
            def bim = file("${params.raw_dir}/${cohort}/${cohort}.bim")
            def fam = file("${params.raw_dir}/${cohort}/${cohort}.fam")
            if (!bed.exists()) error "skip_qc=true but BED not found: ${bed}"
            tuple(cohort, bed, bim, fam)
        }

        // Combined PCA — use params.pca_dir (points to GWAS_analysis/pca/)
        def pca_base = params.pca_dir ?: "${params.out_dir}/pca"
        def evc_file = file("${pca_base}/1000G_pca.eigenvec")
        def eva_file = file("${pca_base}/1000G_pca.eigenval")
        def pcv_file = file("${pca_base}/pca_variance.tsv")
        if (!evc_file.exists()) error "Combined eigenvec not found: ${evc_file}\nSet params.pca_dir in nextflow.config"

        ch_eigenvec = Channel.value(evc_file)
        ch_eigenval = Channel.value(eva_file)
        ch_pca_var  = Channel.value(pcv_file)

        ch_prune_in = Channel.value(
            file("${pca_base}/pruned_global.prune.in").exists()
                ? file("${pca_base}/pruned_global.prune.in")
                : file("${params.raw_dir}/Cohort_A/pruned_Cohort_A.prune.in")  // fallback to gwas_step0 output
        )

        // When skip_qc=true, generate QC metrics on the fly from the raw fam/bim files.
        // This counts samples and variants from the original input data so the
        // QC tab in the HTML report is populated correctly even for external QC runs.
        ch_qc_files = Channel.fromList(cohort_list).map { cohort ->
            def fam = file("${params.raw_dir}/${cohort}/${cohort}.fam")
            def bim = file("${params.raw_dir}/${cohort}/${cohort}.bim")
            if (!fam.exists()) {
                log.warn "QC metrics: fam not found at ${fam} — QC tab will be empty for ${cohort}"
                return null
            }
            def n_samples  = fam.readLines().size()
            def n_variants = bim.readLines().size()
            def out = file("${workDir}/qc_metrics_${cohort}.tsv")
            out.text = "metric\tvalue\tcohort\n" +
                "raw_samples\t${n_samples}\t${cohort}\n" +
                "raw_variants\t${n_variants}\t${cohort}\n" +
                "final_samples\t${n_samples}\t${cohort}\n" +
                "final_variants\t${n_variants}\t${cohort}\n" +
                "maf\t${params.maf}\t${cohort}\n" +
                "hwe\t${params.hwe}\t${cohort}\n" +
                "geno\t${params.geno}\t${cohort}\n" +
                "mind\t${params.mind}\t${cohort}\n" +
                "het_sd\t${params.het_sd}\t${cohort}\n"
            return out
        }.filter { it != null }.collect()
    }

    // ── Phenotypes + covariates (existing files) ──────────────────
    ch_phenotypes = Channel.fromList(cohort_list).map { cohort ->
        def f = file("${params.pheno_dir}/${cohort}/phenotypes.tsv")
        if (!f.exists()) error "phenotypes.tsv not found: ${f}"
        tuple(cohort, f)
    }
    ch_covariates = Channel.fromList(cohort_list).map { cohort ->
        def f = file("${params.pheno_dir}/${cohort}/covariates.tsv")
        if (!f.exists()) error "covariates.tsv not found: ${f}"
        tuple(cohort, f)
    }

    // ── A7  Recode binary phenotypes ──────────────────────────────
    RECODE_BINARY(ch_phenotypes)

    // ── A8–A10  REGENIE + harmonise ───────────────────────────────
    if (!params.skip_regenie) {

        ch_input_qt = ch_qc_beds
            .join(ch_phenotypes, by: 0)
            .join(ch_covariates, by: 0)
            .combine(ch_prune_in)

        ch_input_bt = ch_qc_beds
            .join(RECODE_BINARY.out.bt_phenotypes, by: 0)
            .join(ch_covariates, by: 0)
            .combine(ch_prune_in)

        REGENIE_STEP1_QT(ch_input_qt, 'qt')
        REGENIE_STEP1_BT(ch_input_bt, 'bt')

        REGENIE_STEP2_QT(
            ch_input_qt.join(REGENIE_STEP1_QT.out.step1_pred, by: 0), 'qt')
        REGENIE_STEP2_BT(
            ch_input_bt.join(REGENIE_STEP1_BT.out.step1_pred, by: 0), 'bt')

        HARMONISE_SUMSTATS(
            REGENIE_STEP2_QT.out.step2_out.mix(REGENIE_STEP2_BT.out.step2_out)
        )

        ch_sumstat_files = HARMONISE_SUMSTATS.out.sumstats
            .map { cohort, files -> files }.flatten().collect()
        ch_summary_files = HARMONISE_SUMSTATS.out.gwas_summary.collect()

    } else {
        log.info "⏭  skip_regenie — loading existing sumstats from disk"

        // Look in out_dir first, then in original gwas location
        def ss_pattern = "${params.out_dir}/gwas/**/sumstats/*_sumstats.tsv.gz"
        def ss_pattern2 = "${params.raw_dir}/gwas/**/sumstats/*_sumstats.tsv.gz"

        ch_sumstat_files = Channel.fromPath(ss_pattern)
            .mix(Channel.fromPath(ss_pattern2))
            .unique()
            .collect()
            .ifEmpty { error "No sumstat files found in ${ss_pattern} or ${ss_pattern2}" }

        ch_summary_files = Channel.fromPath("${params.out_dir}/gwas/**/sumstats/*_summary.tsv")
            .mix(Channel.fromPath("${params.raw_dir}/gwas/**/sumstats/*_summary.tsv"))
            .unique().collect().ifEmpty([])
    }

    // ── A11  Interactive HTML report ──────────────────────────────
    if (!params.skip_report) {
        GWAS_REPORT(
            ch_eigenvec,
            ch_eigenval,
            ch_pca_var,
            ch_qc_files,
            ch_sumstat_files,
            ch_summary_files,
            params.author,
            params.affiliation
        )
    }

    // ── A12  Static publication figures ───────────────────────────
    if (!params.skip_report) {
        GENERATE_FIGURES(
            ch_eigenvec,
            ch_eigenval,
            ch_sumstat_files,
            ch_summary_files,
            params.author,
            params.institute ?: "INRS-Centre Armand-Frappier Sante-Biotechnologie",
            params.affiliation
        )
    }
}

workflow.onComplete {
    def status = workflow.success ? "SUCCESS ✓" : "FAILED  ✗"
    log.info """
╔══════════════════════════════════════════════════════════════════╗
║  GWAS Pipeline  ${status.padRight(47)}║
║  Duration : ${workflow.duration.toString().padRight(51)}║
║  HTML     : ${(params.out_dir + '/GWAS_Report.html').padRight(51)}║
║  Figures  : ${(params.out_dir + '/figures/').padRight(51)}║
╚══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
