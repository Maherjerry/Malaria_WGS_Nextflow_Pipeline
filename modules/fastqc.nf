// FastQC module for quality control analysis
process FASTQC_RAW {
    tag "${sample_id}"
    label 'fastqc'
    publishDir "${params.outDIR}/00_qc_raw_reads/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results
    
    script:
    """
    fastqc ${reads}
    """
}

process FASTQC_TRIMMED {
    tag "${sample_id}"
    label 'fastqc'
    publishDir "${params.outDIR}/03_qc_trimmed_reads/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results
    
    script:
    """
    fastqc ${reads}
    """
}