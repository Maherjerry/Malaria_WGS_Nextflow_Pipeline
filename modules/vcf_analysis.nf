process VCF_ANALYSIS {
    tag "vcf_analysis"
    label 'variant_calling'
    publishDir "${params.outDIR}/08_VCF_Analysis_Outputs", mode: 'copy'
    
    input:
    tuple path(vcf), path(vcf_idx)
    
    output:
    path "snps_only.vcf.gz", emit: snps_vcf
    path "snps_only.vcf.gz.csi", emit: snps_vcf_idx
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
    // Define the dollar sign variable to properly escape it in the R script
    def dollar = '$'
    
    """
    #!/usr/bin/env Rscript
    
    # Load libraries
    library(rmarkdown)
    
    # Create R Markdown report
    rmd_content <- c(
      "---",
      "title: \\"VCF Analysis Report\\"",
      "author: \\"Plasmodium falciparum WGS Analysis Pipeline\\"",
      "date: \\"`r Sys.Date()`\\"",
      "output: ",
      "  html_document:",
      "    theme: cosmo",
      "    toc: true",
      "    toc_float: true",
      "---",
      "",
      "```{r setup, include=FALSE}",
      "knitr::opts_chunk${dollar}set(echo = FALSE, warning = FALSE, message = FALSE)",
      "library(ggplot2)",
      "library(dplyr)",
      "library(tidyr)",
      "library(knitr)",
      "```",
      "",
      "## Overview",
      "",
      "This report summarizes the variant calls from Plasmodium falciparum whole genome sequencing.",
      "",
      "## Sample Information",
      "",
      "```{r}",
      "samples <- read.csv(\\"${sample_names}\\", header=FALSE)",
      "colnames(samples) <- c(\\"Sample\\")",
      "kable(samples, caption = \\"List of Samples\\")",
      "```",
      "",
      "## Variant Quality Metrics",
      "",
      "```{r}",
      "qualities <- read.csv(\\"${site_qualities}\\", sep=\\"\\\\t\\", header=FALSE)",
      "colnames(qualities) <- c(\\"CHROM\\", \\"POS\\", \\"REF\\", \\"ALT\\", \\"QUAL\\")",
      "",
      "ggplot(qualities, aes(x=QUAL)) + ",
      "  geom_histogram(bins=50, fill=\\"steelblue\\") +",
      "  theme_minimal() +",
      "  labs(title=\\"Distribution of Variant Quality Scores\\",",
      "       x=\\"Quality Score\\", y=\\"Count\\")",
      "```",
      "",
      "## Missing Data",
      "",
      "```{r}",
      "missing_indv <- read.csv(\\"${missing_indv}\\", sep=\\"\\\\t\\")",
      "",
      "ggplot(missing_indv, aes(x=INDV, y=F_MISS)) +",
      "  geom_col(fill=\\"coral\\") +",
      "  theme_minimal() +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +",
      "  labs(title=\\"Missing Data by Sample\\", ",
      "       x=\\"Sample\\", y=\\"Fraction of Missing Genotypes\\")",
      "```",
      "",
      "## Depth of Coverage",
      "",
      "```{r}",
      "indv_depth <- read.csv(\\"${indv_depth}\\", sep=\\"\\\\t\\")",
      "",
      "ggplot(indv_depth, aes(x=INDV, y=MEAN_DEPTH)) +",
      "  geom_col(fill=\\"forestgreen\\") +",
      "  theme_minimal() +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +",
      "  labs(title=\\"Mean Depth per Sample\\", ",
      "       x=\\"Sample\\", y=\\"Mean Depth\\")",
      "```",
      "",
      "## Heterozygosity",
      "",
      "```{r}",
      "het <- read.csv(\\"${heterozygosity}\\", sep=\\"\\\\t\\")",
      "",
      "ggplot(het, aes(x=INDV, y=F)) +",
      "  geom_col(fill=\\"purple\\") +",
      "  theme_minimal() +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +",
      "  labs(title=\\"Inbreeding Coefficient by Sample\\", ",
      "       x=\\"Sample\\", y=\\"Inbreeding Coefficient\\")",
      "```"
    )
    
    # Write R markdown content to file
    writeLines(rmd_content, "vcf_report.Rmd")
    
    # Render the report
    rmarkdown::render("vcf_report.Rmd", output_file = "vcf_analysis_report.html")
    """
} 