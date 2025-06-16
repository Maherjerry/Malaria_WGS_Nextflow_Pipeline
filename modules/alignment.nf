// Module for DNA sequence alignment and processing
process FASTQ_TO_UBAM {
    tag "${sample_id}"
    label 'alignment'
    publishDir "${params.outDIR}/04_alignment/01_unmapped_bam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.unmapped.bam"), emit: ubam

    script:
    def rgid = "${sample_id}_readgroup1"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" FastqToSam \\
        -F1 ${reads[0]} \\
        -F2 ${reads[1]} \\
        -O ${sample_id}.unmapped.bam \\
        -SM ${sample_id} \\
        -RG ${rgid} \\
        -PL ILLUMINA \\
        -LB lib1 \\
        -PU unit1 \\
        --SORT_ORDER queryname \\
        --CREATE_INDEX false \\
        --VALIDATION_STRINGENCY LENIENT
    """
}

process PREPARE_REFERENCE {
    label 'reference'
    publishDir "${projectDir}/data/reference/genome", mode: 'copy'

    input:
    path(ref_fasta)

    output:
    tuple path("${ref_fasta}"),
          path("${ref_fasta}.fai"),
          path("${ref_fasta}.amb"),
          path("${ref_fasta}.ann"),
          path("${ref_fasta}.bwt"),
          path("${ref_fasta}.pac"),
          path("${ref_fasta}.sa"),
          path("${ref_fasta.baseName}.dict"), emit: reference_bundle

    script:
    """
    # Create BWA index
    bwa index ${ref_fasta}

    # Create FASTA index
    samtools faidx ${ref_fasta}

    # Create sequence dictionary
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CreateSequenceDictionary \\
        -R ${ref_fasta}
    """
}

process ALIGN_READS {
    tag "${sample_id}"
    label 'alignment_high_mem'
    publishDir "${params.outDIR}/04_alignment/02_aligned_bam", mode: 'copy'

    input:
    tuple val(sample_id), path(ubam)
    tuple path(ref_fasta), path(ref_fai), path(ref_amb), path(ref_ann),
          path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"),
          path("${sample_id}.aligned.bai"), emit: aligned_bam

    script:
    // Calculate memory for Java based on task memory, leaving some headroom
    def java_mem = Math.max(2, (task.memory.toGiga() - 1))
    """
    # Convert uBAM to FASTQ
    gatk --java-options "-Xmx${java_mem}g" SamToFastq \\
        -I ${ubam} \\
        --FASTQ ${sample_id}_R1.fastq \\
        --SECOND_END_FASTQ ${sample_id}_R2.fastq \\
        --VALIDATION_STRINGENCY LENIENT

    # Align with BWA
    bwa mem -t ${params.bwa_threads} ${params.bwa_mem_params} ${ref_fasta} \\
        ${sample_id}_R1.fastq ${sample_id}_R2.fastq > ${sample_id}.sam

    # Convert SAM to BAM and merge with uBAM
    gatk --java-options "-Xmx${java_mem}g" MergeBamAlignment \\
        --REFERENCE_SEQUENCE ${ref_fasta} \\
        --UNMAPPED_BAM ${ubam} \\
        --ALIGNED_BAM ${sample_id}.sam \\
        --OUTPUT ${sample_id}.aligned.bam \\
        --CREATE_INDEX true \\
        --ADD_MATE_CIGAR true \\
        --CLIP_ADAPTERS false \\
        --CLIP_OVERLAPPING_READS true \\
        --INCLUDE_SECONDARY_ALIGNMENTS true \\
        --MAX_INSERTIONS_OR_DELETIONS -1 \\
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
        --ATTRIBUTES_TO_RETAIN XS \\
        --VALIDATION_STRINGENCY LENIENT

    # Clean up intermediate files
    rm ${sample_id}.sam ${sample_id}_R*.fastq
    """
}

process MARK_DUPLICATES {
    tag "${sample_id}"
    label 'alignment'
    publishDir "${params.outDIR}/04_alignment/03_marked_duplicates_bam", mode: 'copy', pattern: "*.{bam,bai}"
    publishDir "${params.outDIR}/05_qc_aligned_reads/alignment_metrics", mode: 'copy', pattern: "*.metrics.txt"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bai"), emit: markdup_bam
    path "${sample_id}.markdup.metrics.txt", emit: metrics

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.markdup.bam \\
        -M ${sample_id}.markdup.metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    """
}

process BASE_RECALIBRATOR {
    tag "${sample_id}"
    label 'alignment'
    publishDir "${params.outDIR}/04_alignment/04_recalibrated_bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple path(ref_fasta), path(ref_fai), path(ref_amb), path(ref_ann),
          path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    tuple val(sample_id), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bai"), emit: recal_bam
    path "${sample_id}.recal.table", emit: recal_table

    script:
    """
    # Generate recalibration table
    gatk --java-options "-Xmx${task.memory.toGiga()}g" BaseRecalibrator \\
        -I ${bam} \\
        -R ${ref_fasta} \\
        --known-sites ${params.REFSNP1} \\
        --known-sites ${params.REFSNP2} \\
        --known-sites ${params.REFSNP3} \\
        -O ${sample_id}.recal.table

    # Apply base quality score recalibration
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyBQSR \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${sample_id}.recal.table \\
        -O ${sample_id}.recal.bam
    """
}

process COLLECT_ALIGNMENT_METRICS {
    tag "${sample_id}"
    label 'alignment_metrics'
    publishDir "${params.outDIR}/05_qc_aligned_reads/alignment_metrics", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple path(ref_fasta), path(ref_fai), path(ref_amb), path(ref_ann),
          path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    path "${sample_id}.alignment_metrics.txt", emit: alignment_metrics
    path "${sample_id}.insert_metrics.txt", emit: insert_metrics
    path "${sample_id}.insert_size_histogram.pdf", emit: insert_histogram
    path "${sample_id}.coverage_metrics.txt", emit: coverage_metrics

    script:
    """
    # Collect alignment metrics
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CollectAlignmentSummaryMetrics \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -O ${sample_id}.alignment_metrics.txt

    # Collect insert size metrics
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CollectInsertSizeMetrics \\
        -I ${bam} \\
        -O ${sample_id}.insert_metrics.txt \\
        -H ${sample_id}.insert_size_histogram.pdf

    # Calculate coverage statistics
    samtools coverage ${bam} > ${sample_id}.coverage_metrics.txt
    """
}
