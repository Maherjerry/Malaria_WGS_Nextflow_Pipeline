// MultiQC module for QC reports
process MULTIQC_RAW {
    label 'qc'
    publishDir "${params.outDIR}/00_qc_raw_reads", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", optional: true

    script:
    """
    multiqc . -f
    """
}

process MULTIQC_TRIMMED {
    label 'qc'
    publishDir "${params.outDIR}/03_qc_trimmed_reads", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", optional: true

    script:
    """
    multiqc . -f
    """
}

// Process to run MultiQC on alignment metrics
process ALIGNMENT_MULTIQC {
    label 'qc'
    publishDir "${params.outDIR}/05_qc_aligned_reads", mode: 'copy'

    input:
    path metrics_files

    output:
    path "multiqc_report.html"
    path "multiqc_data", optional: true

    script:
    """
    multiqc . -f
    """
}

process VARIANT_MULTIQC {
    label 'qc'
    publishDir "${params.outDIR}/multiqc/variants", mode: 'copy'

    input:
    tuple path(combined_vcf), path(combined_vcf_idx), path(filtered_vcf), path(filtered_vcf_idx)

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc . -f
    """
}
