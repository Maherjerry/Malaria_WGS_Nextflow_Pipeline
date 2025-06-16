// Module for adapter trimming and quality filtering
process TRIM_GALORE {
    tag "${sample_id}"
    label 'trimming'
    publishDir "${params.outDIR}/01_trimmed_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trim_log

    script:
    """
    trim_galore \\
        --cores ${params.trim_cores} \\
        --paired \\
        --quality ${params.trim_quality} \\
        --phred33 \\
        --stringency ${params.trim_stringency} \\
        --length ${params.trim_length} \\
        --fastqc \\
        --gzip \\
        ${reads[0]} \\
        ${reads[1]}
    """
} 