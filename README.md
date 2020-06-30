# Analysis pipelines for fungal genomic analysis
This repository contains workflows used by the Fungal Genomics group in the Infectious Disease and Microbiome Program (IDMP), the Broad Institute of MIT and Harvard.

Each directory in this repository corresponds to a workflow. Workflows are written in [Workflow Description Language (WDL)](https://github.com/openwdl/wdl). Analysis environment was tracked using [docker](https://www.docker.com). Each workflow could be used on premise with the [Cromwell Engine](https://github.com/broadinstitute/cromwell), or on the cloud with the [Terra platform](https://app.terra.bio).

## Workflow synopsis
* [gatk3](gatk3/README.md): GATK3 fungal variant calling workflow

## About
This repo is developed and maintained by Fungal Genomics Group, IDMP, Broad Institute. Use `issues` tag to report any bugs.

This software is under BSD III license.
