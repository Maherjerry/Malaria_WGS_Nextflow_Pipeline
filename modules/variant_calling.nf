process INDEX_RESOURCE_VCFS {
    publishDir "${params.known_sites}", mode: 'copy'
    
    input:
    tuple val(id), path(vcf)
    
    output:
    tuple val(id), path(vcf), path("${vcf}.tbi")
    
    script:
    """
    gatk IndexFeatureFile -I ${vcf}
    """
} 

process HAPLOTYPE_CALLER {
    tag "${sample_id}"
    label 'variant_calling'
    publishDir "${params.outDIR}/06_variants/01_gvcf", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple path(ref), path(ref_fai), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)
    
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf
    
    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        -R ${ref} \\
        -I ${bam} \\
        -O ${sample_id}.g.vcf.gz \\
        -ERC GVCF

    gatk IndexFeatureFile \\
        -I ${sample_id}.g.vcf.gz
    """
}

process COMBINE_GVCFS {
    label 'variant_calling_high_mem'
    publishDir "${params.outDIR}/06_variants/02_combined_gvcf", mode: 'copy'
    
    input:
    path gvcfs
    path gvcf_indices
    tuple path(ref), path(ref_fai), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)
    
    output:
    tuple path("cohort.g.vcf.gz"), path("cohort.g.vcf.gz.tbi"), emit: combined_gvcf
    
    script:
    def gvcf_params = gvcfs.collect { "-V $it" }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CombineGVCFs \\
        -R ${ref} \\
        ${gvcf_params} \\
        -O cohort.g.vcf.gz
    """
}

process GENOTYPE_GVCFS {
    label 'variant_calling_high_mem'
    publishDir "${params.outDIR}/06_variants/03_final_vcf", mode: 'copy'
    
    input:
    tuple path(gvcf), path(gvcf_idx)
    tuple path(ref), path(ref_fai), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)
    
    output:
    tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi"), emit: final_vcf
    
    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \\
        -R ${ref} \\
        -V ${gvcf} \\
        -O cohort.vcf.gz
    """
}

process VARIANT_RECALIBRATOR {
    label 'variant_calling_high_mem'
    publishDir "${params.outDIR}/06_variants/04_variant_recalibration", mode: 'copy'
    
    input:
    tuple path(vcf), path(vcf_idx)
    tuple path(ref), path(ref_fai), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)
    
    output:
    tuple path("merged_samples.recal"), path("merged_samples.recal.idx"), emit: recal_table
    path "merged_samples.tranches", emit: tranches
    
    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantRecalibrator \\
        -R ${ref} \\
        -V ${vcf} \\
        --resource:cross1,known=false,training=true,truth=true,prior=15.0 ${params.REFSNP1} \\
        --resource:cross2,known=false,training=true,truth=false,prior=12.0 ${params.REFSNP2} \\
        --resource:cross3,known=false,training=true,truth=false,prior=12.0 ${params.REFSNP3} \\
        -an QD \\
        -an MQ \\
        -an ReadPosRankSum \\
        -an FS \\
        -an SOR \\
        -mode SNP \\
        -O merged_samples.recal \\
        --tranches-file merged_samples.tranches
    """
}

process APPLY_VQSR {
    label 'variant_calling'
    publishDir "${params.outDIR}/06_variants/05_filtered_vcf", mode: 'copy'
    
    input:
    tuple path(vcf), path(vcf_idx)
    tuple path(ref), path(ref_fai), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)
    tuple path(recal), path(recal_idx)
    path tranches
    
    output:
    tuple path("filtered_merged_samples.vcf.gz"), path("filtered_merged_samples.vcf.gz.tbi"), emit: filtered_vcf
    
    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyVQSR \\
        -R ${ref} \\
        -V ${vcf} \\
        -O filtered_merged_samples.vcf.gz \\
        --truth-sensitivity-filter-level 99.0 \\
        --tranches-file ${tranches} \\
        --recal-file ${recal} \\
        -mode SNP
    """
}

process BCFTOOLS_STATS {
    label 'variant_calling'
    publishDir "${params.outDIR}/06_variants/06_vcf_stats", mode: 'copy'

    input:
    tuple path(vcf), path(vcf_idx)
    tuple path(ref), path(ref_fai), path(ref_amb), path(ref_ann), 
          path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_dict)

    output:
    path "bcftools_stats.txt", emit: bcftools_stats

    script:
    """
    bcftools stats \
        --fasta-ref ${ref} \
        ${vcf} > bcftools_stats.txt
    """
}

process VARIANT_MULTIQC {
    label 'variant_calling'
    publishDir "${params.outDIR}/07_qc_variants", mode: 'copy'

    input:
    tuple path(combined_vcf), path(combined_vcf_idx)
    tuple path(filtered_vcf), path(filtered_vcf_idx)

    output:
    path "combined_vcf_multiqc.html"
    path "filtered_vcf_multiqc.html"
    path "combined_vcf_multiqc_data", optional: true
    path "filtered_vcf_multiqc_data", optional: true

    script:
    """
    # Create directories for organizing VCF stats
    mkdir -p combined_vcf filtered_vcf

    # Generate stats for combined VCF
    bcftools stats ${combined_vcf} > combined_vcf/bcftools_stats.txt

    # Generate stats for filtered VCF
    bcftools stats ${filtered_vcf} > filtered_vcf/bcftools_stats.txt

    # Generate MultiQC reports
    multiqc combined_vcf -n "combined_vcf_multiqc.html" -f
    multiqc filtered_vcf -n "filtered_vcf_multiqc.html" -f
    """
}

process VCF_ANALYSIS {
    tag "vcf_analysis"
    label 'variant_calling'
    publishDir "${params.outDIR}/08_VCF_Analysis_Outputs", mode: 'copy'
    
    input:
    tuple path(vcf), path(vcf_idx)
    
    output:
    path "snps_only.vcf.gz", emit: snps_vcf
    path "snps_only.vcf.gz.tbi", emit: snps_vcf_idx
    path "genotype_table.csv", emit: genotypes
    path "allelic_depth.csv", emit: allelic_depth
    path "depth_per_site.csv", emit: site_depth
    path "site_qualities.csv", emit: site_qualities
    path "sample_names.csv", emit: sample_names
    path "allele_frequencies.frq", emit: allele_freq
    path "mean_depth_individual.idepth", emit: indv_depth
    path "mean_depth_site.ldepth.mean", emit: mean_site_depth
    path "missing_individuals.imiss", emit: missing_indv
    path "missing_sites.lmiss", emit: missing_sites
    path "heterozygosity.het", emit: heterozygosity
    
    script:
    """
    echo "Starting VCF summary analysis..."

    ##########################
    # BCFTOOLS-based exports
    ##########################

    echo "Extracting SNPs..."
    bcftools view -v snps -Oz -o snps_only.vcf.gz ${vcf}
    bcftools index snps_only.vcf.gz

    echo "Extracting genotype table..."
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ${vcf} > genotype_table.csv

    echo "Extracting allelic depth..."
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' ${vcf} > allelic_depth.csv

    echo "Extracting read depth..."
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DP]\\n' ${vcf} > depth_per_site.csv

    echo "Extracting site quality..."
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\n' ${vcf} > site_qualities.csv

    echo "Extracting sample names..."
    bcftools query -l ${vcf} > sample_names.csv

    ##########################
    # VCFTOOLS-based exports
    ##########################

    echo "Calculating allele frequency..."
    vcftools --gzvcf ${vcf} --freq2 --out allele_frequencies --max-alleles 2

    echo "Calculating mean depth per individual..."
    vcftools --gzvcf ${vcf} --depth --out mean_depth_individual

    echo "Calculating mean depth per site..."
    vcftools --gzvcf ${vcf} --site-mean-depth --out mean_depth_site

    echo "Calculating missing data per individual..."
    vcftools --gzvcf ${vcf} --missing-indv --out missing_individuals

    echo "Calculating missing data per site..."
    vcftools --gzvcf ${vcf} --missing-site --out missing_sites

    echo "Calculating heterozygosity and inbreeding coefficient..."
    vcftools --gzvcf ${vcf} --het --out heterozygosity
    """
}

process VCF_ANALYSIS_REPORT {
    tag "vcf_report"
    label 'variant_calling'
    publishDir "${params.outDIR}/08_VCF_Analysis_Outputs", mode: 'copy'
    
    input:
    path genotypes
    path allelic_depth
    path site_depth
    path site_qualities
    path sample_names
    path allele_freq
    path indv_depth
    path mean_site_depth
    path missing_indv
    path missing_sites
    path heterozygosity
    
    output:
    path "vcf_analysis_report.html", emit: report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load libraries
    library(rmarkdown)
    
    # Create R Markdown report
    cat <<EOF > vcf_report.Rmd
    ---
    title: "VCF Analysis Report"
    author: "Plasmodium falciparum WGS Analysis Pipeline"
    date: "`r Sys.Date()`"
    output: 
      html_document:
        theme: cosmo
        toc: true
        toc_float: true
    ---
    
    ```{r setup, include=FALSE}
    knitr::opts_chunk\$set(echo = FALSE, warning = FALSE, message = FALSE)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(knitr)
    ```
    
    ## Overview
    
    This report summarizes the variant calls from Plasmodium falciparum whole genome sequencing.
    
    ## Sample Information
    
    ```{r}
    samples <- read.csv("${sample_names}", header=FALSE)
    colnames(samples) <- c("Sample")
    kable(samples, caption = "List of Samples")
    ```
    
    ## Variant Quality Metrics
    
    ```{r}
    qualities <- read.csv("${site_qualities}", sep="\t", header=FALSE)
    colnames(qualities) <- c("CHROM", "POS", "REF", "ALT", "QUAL")
    
    ggplot(qualities, aes(x=QUAL)) + 
      geom_histogram(bins=50, fill="steelblue") +
      theme_minimal() +
      labs(title="Distribution of Variant Quality Scores",
           x="Quality Score", y="Count")
    ```
    
    ## Missing Data
    
    ```{r}
    missing_indv <- read.csv("${missing_indv}", sep="\t")
    
    ggplot(missing_indv, aes(x=INDV, y=F_MISS)) +
      geom_col(fill="coral") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title="Missing Data by Sample", 
           x="Sample", y="Fraction of Missing Genotypes")
    ```
    
    ## Depth of Coverage
    
    ```{r}
    indv_depth <- read.csv("${indv_depth}", sep="\t")
    
    ggplot(indv_depth, aes(x=INDV, y=MEAN_DEPTH)) +
      geom_col(fill="forestgreen") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title="Mean Depth per Sample", 
           x="Sample", y="Mean Depth")
    ```
    
    ## Heterozygosity
    
    ```{r}
    het <- read.csv("${heterozygosity}", sep="\t")
    
    ggplot(het, aes(x=INDV, y=F)) +
      geom_col(fill="purple") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title="Inbreeding Coefficient by Sample", 
           x="Sample", y="Inbreeding Coefficient")
    ```
    
    EOF
    
    # Render the report
    rmarkdown::render("vcf_report.Rmd", output_file = "vcf_analysis_report.html")
    """
}

