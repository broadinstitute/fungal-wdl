# GATK3 variant calling workflow
This repository contains GATK3 based variant calling pipelines used by the Fungal Genomics group in the Infectious Disease and Microbiome Program, the Broad Institute of MIT and Harvard.

The pipeline input is set of unaligned (or aligned) BAM files and a reference genome. Reads are converted to FASTQ and aligned with BWA, alignments sorted and marked for duplicates with Picard, SNPs joint-called using GATK HaplotypeCaller, and hard filtration using GATK variant filtration. The output is a multi-sample VCF file.

## Synopsis
* README.md: this file
* docs: documentations and tutorials
* tests: files used for testing the workflow
* workflows: variant calling workflow. Note that `fungal_variant_calling_gatk3.wdl` is the full WDL, while `haplotype_caller_scatter_gatk3.wdl` is a module that the former uses.
* Dockerfile: docker image for the workflow. The docker is also available from dockerhub [here](https://hub.docker.com/r/broadinstitute/fungi-gatk3).

This method was deployed on Terra Method Repository [here](https://portal.firecloud.org/?return=firecloud#methods/broad-fungal-methods/funga-variant-call-gatk3/1). Example Terra workspace can be found [here](https://firecloud.terra.bio/#workspaces/broad-fungal-firecloud/broad-fungal-gatk3). 

## Software versions
```sh
PICARD_VER=1.782
GATK37_VER=3.7-93-ge9d8068
SAMTOOLS_VER=1.3.1
BWA_VER=0.7.12
TABIX_VER=0.2.5_r1005
BGZIP_VER=1.3
```

## About
This repo is developed by Xiao Li from IDMP, the Broad Institute. Use `issues` tag to report any bugs.
