# Sunflower Variant Calling Pipeline

This repository contains a pipeline for variant calling in sunflower genomes using a series of PBS scripts. The pipeline is split into three main stages: trimming, alignment, and duplication removal; variant calling with GATK HaplotypeCaller; and joint genotyping across multiple samples. 

## Pipeline Overview

The pipeline consists of the following modules:

1. **Trim, Align, and Mark Duplicates (01_trim_align_markdup.pbs)**:
   - Trims raw FASTQ data.
   - Aligns the data to the sunflower reference genome (Ha 412 HO v2).
   - Marks duplicate reads and collects alignment and depth metrics.

2. **Variant Calling (02_HaplotypeCaller.pbs)**:
   - Calls variants using GATK's HaplotypeCaller in GVCF mode.
   - Parallelizes variant calling across chromosomes using GNU parallel.

3. **GenomicsDB Import and Joint Genotyping (03_GenomicsDBImport_GenotypeGVCFs.pbs)**:
   - Imports GVCF files into a GenomicsDB for joint genotyping.
   - Produces cohort VCF files, one for each chromosome.

## Requirements

- PBS job scheduler (for job array management).
- GATK (Genome Analysis Toolkit).
- GNU parallel for parallelization.
- FASTQ data, BAM files, and GVCF files as inputs.
- fastp
- bwa-mem2
- mosdepth_d4
- gatk v4.6.0

## Detailed Description of Each Script

### 1. 01_trim_align_markdup.pbs

This script handles the following steps:

- **Required Inputs**:
  - `SAMPLE_LIST`: List of samples to process.
  - `FASTQ_DIR`: Directory containing raw FASTQ files.
  - `OUT_DIR`: Directory for output files.
  
- **Key Points**:
  - A PBS job array is used to process each sample in parallel.
  - The number of threads is set to 16 per job, as this is the maximum number supported by `fastp`.
  - The job array size must be updated manually in the PBS header based on the number of samples.

### 2. 02_HaplotypeCaller.pbs

This script performs variant calling using GATK's HaplotypeCaller in GVCF mode.

- **Required Inputs**:
  - `SAMPLE_LIST`: List of samples.
  - `BAM_DIR`: Directory containing the BAM files from the previous step.
  - `OUT_DIR`: Output directory for the GVCF files.

- **Key Points**:
  - The script parallelizes variant calling across chromosomes using GNU parallel, with 17 threads per job (one for each chromosome).
  - It skips variant calling for the chloroplast, mitochondria, and unplaced scaffolds.
  - Repetitive regions are masked using a BED file from an updated repeat annotation.
  - The heterozygosity parameters (`--heterozygosity` and `--heterozygosity-stdev`) are increased to reflect sunflower population genetic diversity, as per Todesco et al. 2020.

### 3. 03_GenomicsDBImport_GenotypeGVCFs.pbs

This script imports individual GVCF files into a GenomicsDB and performs joint genotyping.

- **Required Inputs**:
  - `SAMPLE_LIST`: List of samples to process.
  - `GVCF_DIR`: Directory containing GVCF files.
  - `COHORT`: Cohort name, used as a prefix for final VCF files.
  - `OUT_DIR`: Output directory for the final VCF files.

- **Key Points**:
  - A PBS job array is used, but each job corresponds to a chromosome, not a sample.
  - This step generates one VCF file per chromosome for the entire cohort.

## Usage

### Script 01: `01_trim_align_markdup.pbs`

1. Prepare a `SAMPLE_LIST` containing the sample names you want to process.
2. Set the `FASTQ_DIR` and `OUT_DIR` variables.
3. Update the job array size in the PBS header for the number of samples.
4. Submit the script to the PBS scheduler.

### Script 02: `02_HaplotypeCaller.pbs`

1. Set the `SAMPLE_LIST`, `BAM_DIR`, and `OUT_DIR` variables.
2. Submit the script to the PBS scheduler. The script will run parallel jobs across chromosomes.

### Script 03: `03_GenomicsDBImport_GenotypeGVCFs.pbs`

1. Set the `SAMPLE_LIST`, `GVCF_DIR`, `COHORT`, and `OUT_DIR` variables.
2. Submit the script to the PBS scheduler for joint genotyping across chromosomes.

## Outputs

- **01_trim_align_markdup.pbs**: BAM files, alignment metrics, and depth metrics.
- **02_HaplotypeCaller.pbs**: GVCF files for each sample and chromosome.
- **03_GenomicsDBImport_GenotypeGVCFs.pbs**: Joint genotyped VCF files for each chromosome in the cohort.

## Notes

- Ensure that the `SAMPLE_LIST` and directory paths are correctly set before running each script.
- This pipeline is designed for processing sunflower genomic data, but it can be adapted for other species with appropriate reference genomes and modifications.
- If you have any issues or suggestions for improvements, please open an issue in this repository.

## Contributing
Contributions are welcome! For major changes, please open an issue first to discuss what you would like to change. Please ensure to update tests and documentation as appropriate.

## License
This work is released into the public domain under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication.
To the extent possible under law, all copyright and related or neighboring rights to this work have been waived. This work is published from the United States.
For more information, please visit: https://creativecommons.org/publicdomain/zero/1.0/
