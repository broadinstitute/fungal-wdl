# Variant calling pipeline for fungal haploid genomes
# Developed by Xiao Li (xiaoli@broadinstitute.org) and Aina Martinez Zurita (amartine@broadinstitute.org)
# Fungal Genomics Group, Infectious Disease and Microbiome Program.
# The Broad Institute of MIT and Harvard


workflow GATK4_Germline_Variants {

  # Config Parameters

  String gatk_docker
  String python2_docker
  String gatk_path
  Int preemptible_general
  Int preemptible_HC

  # Reference files
  File ref
  File dict
  File amb
  File ann
  File bwt
  File fai
  File pac
  File sa

  # Input Files
  Array[String] input_samples
  Array[File] input_uBAMs

  scatter(i in range(length(input_samples))) {
    String sample_name = input_samples[i]
    String input_bam = input_uBAMs[i]

    #Document Tool
    call MarkIlluminaAdapters {
        input:
          input_bam=input_bam
          input_sample_name=sample_name
          gatk_path=gatk_path
          docker=docker
          preemptible=preemptible_general
    }

    #Document tool
    call SamToFastqAllignMerge {
      input:

    }

    call MarkDuplicates{

    }

    call ReorderBam {


    }

    call HaplotypeCaller {


    }


  }

  call CombineGVCFs {

  }

  call GenotypeGVCFs {


  }

  call HardFiltration {


  }


  call CostumFiltration {


  }


  output {

  }
}


### Task Definitions ###

task MarkIlluminaAdapters {
  File input_bam
  String input_sample_name

  String gatk_path
  String docker
  Int mem_size_gb
  Int preemptible
  Int disk_size


  command {
    # e - exit when a command fails
    # u - exit when script tries to use undeclared variables
    # x - trace what gets executed
    # o - exit status of the last command that threw a non-zero exit code is returned
    set -euxo pipefail

    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MarkIlluminaAdapters -I=${input_bam} -O=${input_sample_name}_markilluminaadapters.bam -M=${input_sample_name}_markilluminaadapters_metrics.txt

  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk "+ disk_size + " HDD"
  }

  output {
    File bam = "${input_sample_name}_markilluminaadapters.bam"
    File metrics_adapters = "${input_sample_name}_markilluminaadapters_metrics.txt"

  }
}


task SamToFastqAllignMerge {
  File input_bam
  String input_sample_name

  String docker
  String gatk_path

  Int disk_size
  Int mem_size_gb
  Int preemptible

  File ref
  File dict
  File amb
  File ann
  File bwt
  File fai
  File pac
  File sa

  String read_group = "'@RG\\tID:FLOWCELL_${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:LIB_${sample_name}'"

  command {
    # e - exit when a command fails
    # u - exit when script tries to use undeclared variables
    # x - trace what gets executed
    # o - exit status of the last command that threw a non-zero exit code is returned
    set -euxo pipefail

    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" SamToFastq -I=${input_bam} \
    --FASTQ=/dev/stdout --CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | \
    bwa mem -M -t 7 -p /path/Homo_sapiens_assembly19.fasta /dev/stdin

  }

  runtime {


  }

  output {

  }


}
