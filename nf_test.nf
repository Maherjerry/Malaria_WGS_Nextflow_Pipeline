#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Version check
if( !nextflow.version.matches('>= 24.10') ) {
    error "This workflow requires Nextflow version 24.10 or greater"
}

// Input validation
def validateInputs() {
    // Check input directory
    def readDir = file(params.readDIR)
    if (!readDir.exists()) {
        error "Input directory not found: ${params.readDIR}"
    }

    // Check reference file
    def refFile = file(params.fasta)
    if (!refFile.exists()) {
        error "Reference file not found: ${params.fasta}"
    }

    // Check known sites directory
    def knownSitesDir = file(params.known_sites)
    if (!knownSitesDir.exists()) {
        error "Known sites directory not found: ${params.known_sites}"
    }

    log.info "Input validation passed"
}

// Print pipeline header
log.info """
===========================================
 🧬  WGS ANALYSIS PIPELINE  🧬
===========================================

Pipeline Parameters:
-------------------------------------------
📁 Input Directory : ${params.readDIR}
📂 Output Directory: ${params.outDIR}
📂 Reference File  : ${params.fasta}
📂 Known Sites Dir : ${params.known_sites}
💻 CPUs           : ${params.cpus}
-------------------------------------------
"""

// Import modules
include { FASTQC_RAW; FASTQC_TRIMMED } from './modules/fastqc'
include { MULTIQC_RAW; MULTIQC_TRIMMED; ALIGNMENT_MULTIQC } from './modules/multiqc'
include { TRIM_GALORE } from './modules/trimming'
include { FASTQ_TO_UBAM; PREPARE_REFERENCE; ALIGN_READS;
         MARK_DUPLICATES; BASE_RECALIBRATOR; COLLECT_ALIGNMENT_METRICS } from './modules/alignment'
include { HAPLOTYPE_CALLER; COMBINE_GVCFS; GENOTYPE_GVCFS;
         VARIANT_RECALIBRATOR; APPLY_VQSR } from './modules/variant_calling'
include { VCF_ANALYSIS; VCF_ANALYSIS_REPORT } from './modules/vcf_analysis'

// Main workflow
workflow {
    // Validate inputs
    validateInputs()

    // Input channel for FASTQ files
    Channel
        .fromFilePairs("${params.readDIR}/*_R{1,2}.fastq.gz", flat: true)
        .ifEmpty { error "No FASTQ file pairs found in ${params.readDIR}" }
        .set { FASTQ_FILES }

    // Reference channel
    reference_ch = Channel.fromPath(params.fasta, checkIfExists: true)

    // Known sites channel
    known_sites_ch = Channel
        .fromPath("${params.known_sites}/*.vcf.gz")
        .collect()

    // QC and Trimming
    FASTQC_RAW(FASTQ_FILES.map { it -> [it[0], [it[1], it[2]]] })
    MULTIQC_RAW(FASTQC_RAW.out.fastqc_results.collect())

    TRIM_GALORE(FASTQ_FILES.map { it -> [it[0], [it[1], it[2]]] })
    
    // Reshape the trimmed reads channel to convert [sample_id, read1, read2] to [sample_id, [read1, read2]]
    trimmed_reads_ch = TRIM_GALORE.out.trimmed_reads.map { 
        sample_id, read1, read2 -> [sample_id, [read1, read2]] 
    }
    
    FASTQC_TRIMMED(trimmed_reads_ch)
    MULTIQC_TRIMMED(FASTQC_TRIMMED.out.fastqc_results.collect())

    // Alignment steps
    PREPARE_REFERENCE(reference_ch)
    FASTQ_TO_UBAM(trimmed_reads_ch)

    ALIGN_READS(
        FASTQ_TO_UBAM.out.ubam,
        PREPARE_REFERENCE.out.reference_bundle.collect()
    )

    MARK_DUPLICATES(ALIGN_READS.out.aligned_bam)

    BASE_RECALIBRATOR(
        MARK_DUPLICATES.out.markdup_bam,
        PREPARE_REFERENCE.out.reference_bundle.collect()
    )

    COLLECT_ALIGNMENT_METRICS(
        BASE_RECALIBRATOR.out.recal_bam,
        PREPARE_REFERENCE.out.reference_bundle.collect()
    )

    // Alignment QC
    ALIGNMENT_MULTIQC(
        COLLECT_ALIGNMENT_METRICS.out.alignment_metrics
            .mix(COLLECT_ALIGNMENT_METRICS.out.insert_metrics)
            .mix(MARK_DUPLICATES.out.metrics)
            .mix(COLLECT_ALIGNMENT_METRICS.out.coverage_metrics)
            .collect()
    )

    // Variant Calling steps
    HAPLOTYPE_CALLER(
        BASE_RECALIBRATOR.out.recal_bam,
        PREPARE_REFERENCE.out.reference_bundle.collect()
    )

    COMBINE_GVCFS(
        HAPLOTYPE_CALLER.out.gvcf.map { it[1] }.collect(),
        HAPLOTYPE_CALLER.out.gvcf.map { it[2] }.collect(),
        PREPARE_REFERENCE.out.reference_bundle.collect()
    )

    GENOTYPE_GVCFS(
        COMBINE_GVCFS.out.combined_gvcf,
        PREPARE_REFERENCE.out.reference_bundle.collect()
    )

    VARIANT_RECALIBRATOR(
        GENOTYPE_GVCFS.out.final_vcf,
        PREPARE_REFERENCE.out.reference_bundle.collect()
    )

    APPLY_VQSR(
        GENOTYPE_GVCFS.out.final_vcf,
        PREPARE_REFERENCE.out.reference_bundle.collect(),
        VARIANT_RECALIBRATOR.out.recal_table,
        VARIANT_RECALIBRATOR.out.tranches
    )

    // After APPLY_VQSR process
    VCF_ANALYSIS(APPLY_VQSR.out.filtered_vcf)
    VCF_ANALYSIS_REPORT(
        VCF_ANALYSIS.out.genotypes,
        VCF_ANALYSIS.out.allelic_depth,
        VCF_ANALYSIS.out.site_depth,
        VCF_ANALYSIS.out.site_qualities,
        VCF_ANALYSIS.out.sample_names,
        VCF_ANALYSIS.out.allele_freq,
        VCF_ANALYSIS.out.indv_depth,
        VCF_ANALYSIS.out.mean_site_depth,
        VCF_ANALYSIS.out.missing_indv,
        VCF_ANALYSIS.out.missing_sites,
        VCF_ANALYSIS.out.heterozygosity
    )
}
