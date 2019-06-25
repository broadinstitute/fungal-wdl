# Analysis pipelines for fungal genomic analysis
This repository contains GATK3 based variant calling pipelines used by the Fungal Genomics group in the Infectious Disease and Microbiome Program, the Broad Institute of MIT and Harvard.

The pipeline input is set of unaligned (or aligned) bam files and a reference genome. Reads are converted to fastq and aligned with BWA, alignments sorted and marked for duplicates with Picard, SNPs joint-called using GATK HaplotypeCaller, and hard filtration using GATK variant filtration. The output is a multi-sample VCF file.

## Synopsis
* README.md: this file
* docs: documentations and tutorials
* tests: files used for testing the workflow
* workflows: variant calling workflow. Note that `fungal_variant_calling_gatk3.wdl` is the full WDL, while `haplotype_caller_scatter_gatk3.wdl` is a module that the former uses.
* Dockerfile: docker image for the workflow
* LICENSE: MIT license

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

This software is under MIT license.
