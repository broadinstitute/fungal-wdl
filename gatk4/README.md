# Analysis pipelines for fungal genomic analysis
This subfolder contains the GATK4-based variant calling pipelines used by the Fungal Genomics group in the Infectious Disease and Microbiome Program, the Broad Institute of MIT and Harvard.

The pipeline input is set of unaligned bam files and a reference genome. Reads are converted to fastq and aligned with BWA, alignments sorted and marked for duplicates and adapters with GATK, SNPs joint-called using GATK HaplotypeCaller, and hard filtration using GATK variant filtration, followed by a custom in house filtering script. The output is a multi-sample VCF file.

## Synopsis
* README.md: this file
* docs: Tutorial, and example input files and variables.
* workflows: variant calling workflow. See Workflows for description of all the possible workflow options, and the best practices.
* Dockerfiles: docker image for the workflow, as well as the python docker image use in the custom filtering step.
* LICENSE: MIT license

## Software versions
```sh
PICARD_VER=2.21.2
GATK37_VER=v4.1.4.1
SAMTOOLS_VER=1.1
BWA_VER=0.7.12
```

## Workflows

There are two possible options to run the pipeline.
1. Run the whole pipeline in one go: Run `fungal_variant_calling-gatk4FULL` on a participant set.
2. Run the pipeline in two steps: First, run the `fungal_variant_calling-gatk4Part1` on individual participants to call variants using HaplotypeCaller. Then, form a participant set and produce a joint called vcf with `fungal_variant_calling-gatk4Part2`. This allows samples to added to joint-calls at a later date without having to run the computationally intensive HaplotypeCaller step on the whole set again.

## About
This pipeline is developed by Aina Martinez Zurita, after the original GATK3 pipeline by Xiao Li from IDMP, the Broad Institute. Use `issues` tag to report any bugs.

This software is under MIT license.
