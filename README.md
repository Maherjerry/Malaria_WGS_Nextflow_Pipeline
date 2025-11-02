# Plasmodium falciparum WGS Analysis Pipeline 🧬

A Nextflow pipeline for processing Whole Genome Sequencing (WGS) data from the malaria parasite Plasmodium falciparum, implementing GATK best practices from raw FASTQ to variant calling.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10-brightgreen.svg)](https://www.nextflow.io/)
[![GATK](https://img.shields.io/badge/GATK-4.2+-blue.svg)](https://gatk.broadinstitute.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This pipeline performs comprehensive analysis of Plasmodium falciparum Whole Genome Sequencing data, following industry-standard best practices. It's designed to be user-friendly while maintaining high-quality results for both small and large-scale sequencing projects, with optimizations specific to P. falciparum genomics.
![Overview Diagram](https://github.com/Maherjerry/Malaria_WGS_Nextflow_Pipeline/blob/main/images/Untitled%20Diagram3.jpg)

### Key Features
- 🔄 Complete end-to-end processing from FASTQ to filtered VCF
- 📊 Automated quality control at multiple pipeline stages
- 🛠️ Industry-standard tools (GATK4, BWA-MEM, FastQC)
- 📈 Detailed quality metrics and visualization
- ⚡ Parallelized execution for optimal performance
- 🔧 Modular design for easy customization
- 💻 Compatible with local, HPC, and cloud environments
- 🦟 Optimized for P. falciparum genome characteristics

### Workflow Steps
1. **Quality Control** (FastQC, MultiQC) - Comprehensive QC of raw sequence data
2. **Read Trimming** (Trim Galore) - Adapter removal and quality-based trimming
3. **Alignment** (BWA-MEM) - High-performance read mapping to P. falciparum reference genome
4. **Read Processing** (GATK) - Marking duplicates and base quality recalibration
5. **Variant Calling** (GATK HaplotypeCaller) - Sensitive variant detection
6. **Joint Genotyping** - Population-based variant refinement
7. **Variant Filtering** (VQSR) - Statistical variant quality recalibration
8. **VCF Analysis** - Comprehensive analysis and visualization of variant data

## Installation

### Prerequisites

#### Required Software
- [Nextflow](https://www.nextflow.io/) (>=24.10)
- [Java](https://www.oracle.com/java/technologies/downloads/) (>=11)
- [GATK4](https://github.com/broadinstitute/gatk/releases) (>=4.2.0)
- [BWA](https://github.com/lh3/bwa) (>=0.7.17)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (>=0.11.9)
- [MultiQC](https://multiqc.info/) (>=1.12)
- [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (>=0.6.7)
- [Samtools](https://www.htslib.org/) (>=1.13)
- [BCFtools](https://samtools.github.io/bcftools/) (>=1.13)
- [VCFtools](https://vcftools.github.io/) (>=0.1.16)
- [R](https://www.r-project.org/) (>=4.0) with packages: ggplot2, dplyr, tidyr, knitr, rmarkdown
- [Pandoc](https://pandoc.org/) (>=1.12.3) - Required for R Markdown reports

#### Hardware Recommendations
- **RAM**: 8GB minimum, 16GB+ recommended
- **CPU**: 4+ cores recommended
- **Storage**: 10GB+ per sample (P. falciparum has a ~23Mb genome, significantly smaller than human)
- **Network**: Fast connection for reference data download

### Quick Installation

#### Using Conda/Mamba (Recommended)
```bash
# Create environment with all dependencies
conda create -n pf-wgs-pipeline -c bioconda -c conda-forge nextflow=24.10 gatk4 bwa fastqc multiqc trim-galore samtools bcftools vcftools r-ggplot2 r-dplyr r-tidyr r-knitr r-rmarkdown pandoc
conda activate pf-wgs-pipeline
```

#### Manual Installation
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

# Install other dependencies through system package manager or manually
# Instructions vary by system
```

### Getting the Pipeline
```bash
# Clone the repository
git clone https://github.com/username/pf-wgs-pipeline.git
cd pf-wgs-pipeline

# Test installation
nextflow run nf_test.nf --help
```

## Usage

### Input Requirements

#### FASTQ Files
- **Format**: Paired-end, gzipped FASTQ files (`.fastq.gz`)
- **Naming Convention**: `*_R{1,2}.fastq.gz` (e.g., `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`)
- **Quality Encoding**: Phred+33 (Illumina 1.8+)

#### Reference Genome
- **Required Files**: P. falciparum reference genome FASTA (`.fasta` or `.fa`)
- **Recommended**: Pre-indexed reference for BWA (will be created if not provided)
- **Known Sites**: VCF files for known P. falciparum variants (required for BQSR)
- **Default Reference**: PlasmoDB P. falciparum 3D7 reference

#### Directory Structure
```
data/
├── fastq/                       # Input FASTQ files
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
│
└── reference/
    ├── genome/                  # Reference genome
    │   └── Pf_ref.fasta         # P. falciparum reference
    │
    └── known_sites/             # Known variant sites
        ├── 3d7_hb3.combined.final.vcf.gz
        ├── 7g8_gb4.combined.final.vcf.gz
        └── hb3_dd2.combined.final.vcf.gz
```

### Running the Pipeline

#### Basic Run
```bash
nextflow run nf_test.nf --readDIR /path/to/fastq --outDIR /path/to/results
```

#### With Custom Parameters
```bash
nextflow run nf_test.nf \
    --readDIR /path/to/fastq \
    --outDIR /path/to/results \
    --fasta /path/to/pf_reference.fasta \
    --known_sites /path/to/known_sites \
    --cpus 8 \
    --bwa_threads 4 \
    --trim_quality 20
```

#### On HPC (SLURM)
```bash
nextflow run nf_test.nf \
    -profile cluster \
    --readDIR /path/to/fastq \
    --outDIR /path/to/results
```

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--readDIR` | Directory containing input FASTQ files | `data/fastq` |
| `--outDIR` | Output directory for results | `results` |
| `--fasta` | Path to P. falciparum reference genome FASTA | `data/reference/genome/Pf_ref.fasta` |
| `--known_sites` | Directory with known P. falciparum variants | `data/reference/known_sites` |
| `--cpus` | Number of CPUs for the pipeline | `4` |
| `--trim_quality` | Quality threshold for read trimming | `20` |
| `--trim_length` | Minimum read length after trimming | `20` |
| `--bwa_threads` | Threads for BWA alignment | `4` |
| `--variant_recalibrator_sensitivity` | VQSR sensitivity threshold | `99.0` |

For a complete list of parameters, see `nextflow_schema.json` or run with `--help`.

## Pipeline Details

### Module Structure
The pipeline is organized into modules for easy maintenance and customization:

```
modules/
├── alignment.nf      # Alignment and post-processing
├── fastqc.nf         # Quality control
├── multiqc.nf        # Summary reports
├── trimming.nf       # Read trimming and filtering
└── variant_calling.nf # Variant detection and filtering
```

### Process Flow
1. **Input Validation**: Checks input files and reference
2. **Read QC**: Quality assessment of raw FASTQ files
3. **Trimming**: Adapter and quality trimming
4. **Post-Trim QC**: Quality validation after trimming
5. **Alignment**: Reference mapping to P. falciparum genome, followed by various BAM processing steps:
   - Converting to unmapped BAM
   - Aligning with BWA
   - Marking duplicates
   - Base quality recalibration
6. **Alignment QC**: Quality metrics for aligned reads
7. **Variant Calling**: Per-sample variant identification (GVCF mode)
8. **Joint Genotyping**: Multi-sample variant consolidation
9. **Variant Filtering**: Quality-based variant refinement
10. **Final QC**: Quality assessment of variant calls

### Output Structure

```
results/
├── 00_qc_raw_reads/          # Initial FastQC results
├── 01_trimmed_reads/         # Processed FASTQ files
├── 03_qc_trimmed_reads/      # Post-trimming QC
├── 04_alignment/             # Aligned reads (BAM files)
│   ├── 01_unmapped_bam/      # Unmapped reads
│   ├── 02_aligned_bam/       # Aligned BAM files
│   ├── 03_marked_duplicates_bam/ # Duplicate-marked BAMs
│   └── 04_recalibrated_bam/  # Recalibrated BAMs
├── 05_qc_aligned_reads/      # Alignment quality metrics
├── 06_variants/              # Called variants
│   ├── 01_gvcf/              # Per-sample GVCF files
│   ├── 02_combined_gvcf/     # Multi-sample GVCF
│   ├── 03_final_vcf/         # Raw variant calls
│   ├── 04_variant_recalibration/ # VQSR files
│   └── 05_filtered_vcf/      # Final filtered variants
├── 07_qc_variants/           # Variant quality reports
└── pipeline_info/            # Execution reports
```

### Key Output Files
- **Trimmed Reads**: `01_trimmed_reads/*.fq.gz`
- **Aligned Reads**: `04_alignment/04_recalibrated_bam/*.recal.bam`
- **Variant Calls**: `06_variants/05_filtered_vcf/filtered_merged_samples.vcf.gz`
- **QC Reports**: `00_qc_raw_reads/multiqc_report.html`, `05_qc_aligned_reads/multiqc_report.html`, etc.
- **Pipeline Reports**: `pipeline_info/pipeline_report.html`, `pipeline_info/timeline.html`

## Performance and Optimization

### Resource Usage
- **Memory**: Peak usage typically during variant calling (~4GB per sample for P. falciparum)
- **CPU**: Most intensive during alignment and variant calling
- **Disk I/O**: Highest during BAM sorting and merging
- **Time**: ~1-3 hours for P. falciparum WGS at 50x coverage (system-dependent)

### Optimization Tips
1. **Increase Parallelization**: Set appropriate `--cpus` and `--bwa_threads`
2. **Memory Management**: Adjust Java heap size for GATK in config
3. **Disk Performance**: Use SSD for work directory when possible
4. **Resume Failed Runs**: Use `-resume` to continue from last successful step

## Troubleshooting

### Common Issues

#### Insufficient Memory
```
ERROR ~ Error executing process > 'ALIGN_READS'
Caused by: Process requirement exceeds available memory
```

**Solution**: Reduce Java memory in config or increase available system memory

#### Reference Indexing Fails
```
Error: Reference genome index not found
```

**Solution**: Manually create indices or check reference path

#### VCF Analysis fails with missing .tbi index
```
Missing output file(s) `snps_only.vcf.gz.tbi` expected by process `VCF_ANALYSIS`
```

**Solution**: This was fixed by modifying the VCF_ANALYSIS process to either:
- Expect a .csi index file (bcftools default) instead of .tbi
- Add the `-t` flag to bcftools index command to create .tbi index

#### R Markdown report generation fails
```
Error: pandoc version 1.12.3 or higher is required and was not found
```

**Solution**: Install pandoc using your package manager:
```bash
sudo apt-get install pandoc  # Debian/Ubuntu
```
Then verify installation:
```bash
pandoc --version
Rscript -e "rmarkdown::pandoc_available()"  # Should return TRUE
```

## Advanced Usage

### Custom Reference Genomes
```bash
# Using a different P. falciparum reference strain
nextflow run nf_test.nf \
    --fasta /path/to/pf_dd2.fasta \
    --readDIR /path/to/fastq \
    --known_sites /path/to/custom_known_sites
```
### Cloud Execution......PENDING!!!
<!-- ```bash
# AWS Batch example
nextflow run nf_test.nf \
    -profile aws \
    --readDIR s3://mybucket/fastq \
    --outDIR s3://mybucket/results
``` -->

### Resource Configuration
Edit `nextflow.config` to customize resource allocation:
```nextflow
process {
    withLabel: 'alignment_high_mem' {
        cpus = 8
        memory = '16 GB'
    }
}
```


## Interpreting Results

### Quality Assessment
Review the MultiQC reports to evaluate:
- **Read Quality**: Phred scores, GC content, duplication rates
- **Alignment Quality**: Mapping rate, coverage, insert size distribution
- **Variant Quality**: Ti/Tv ratio, depth distribution, quality scores

### Key Metrics for Success
- **Mapping Rate**: >90% for good quality P. falciparum data
- **Mean Coverage**: >30x recommended for P. falciparum
- **Duplication Rate**: <20% for PCR-free libraries
- **GC Content**: ~19-20% for P. falciparum (AT-rich genome)

### Future Downstream Analysis
The filtered VCF file can be used for:
- Drug resistance marker identification
- Population genetics studies
- Transmission studies
- Genome diversity analysis

<!-- ## Development and Contributing

### How to Contribute
1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Commit changes: `git commit -m 'Add amazing feature'`
4. Push to branch: `git push origin feature/amazing-feature`
5. Open a Pull Request -->

<!-- ### Testing
```bash
# Run test dataset
nextflow run nf_test.nf -profile test
```

### Code Structure
- **Main Workflow**: `nf_test.nf`
- **Configuration**: `nextflow.config`
- **Process Modules**: `modules/`
- **Parameter Schema**: `nextflow_schema.json`

## Version History

| Version | Date | Key Changes |
|---------|------|-------------|
| 1.0.0   | 2023-10-01 | Initial release |
| 1.1.0   | 2024-01-15 | Added VQSR, improved memory management |
| 1.2.0   | 2024-04-23 | Updated to Nextflow 24.10, fixed channel issues | -->

## Citations
1. Sanger's parasite_wgs_genotyping pipeline
2. An Optimized GATK4 Pipeline for Plasmodium falciparum Whole Genome Sequencing Variant Calling and Analysis
3. GATK best practices
<!-- If you use this pipeline in your research, please cite:

```
Author, A. (2024). Plasmodium falciparum WGS Analysis Pipeline: A Nextflow implementation of GATK best practices.
GitHub: https://github.com/yourusername/pf-wgs-pipeline
```

### Tools Used
- GATK: Van der Auwera et al. (2013). From FastQ data to high confidence variant calls. Current Protocols in Bioinformatics, 43, 11.10.1-11.10.33
- BWA-MEM: Li H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2
- Nextflow: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319

## License

This project is licensed under the MIT License - see the LICENSE file for details. -->

## Contact and Support

- **Issues**: Please use the GitHub issue tracker
- **Email**: muhabah@mrc.gm
<!-- - **Twitter**: [@YourTwitterHandle](https://twitter.com/YourTwitterHandle) -->

---
Updated: April 23, 2025

