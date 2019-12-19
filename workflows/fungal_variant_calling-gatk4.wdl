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
  Int general_disk_size
  Int general_mem_size_gb

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
          docker=gatk_docker
          preemptible=preemptible_general
          disk_size=general_disk_size
          mem_size_gb=general_mem_size_gb
    }

    #Document tool
    call SamToFastqAllignMerge {
      input:
        input_bam=MarkIlluminaAdapters.bam
        input_sample_name=sample_name
        input_unmapped_bam=input_bam
        gatk_path=gatk_path
        docker=gatk_docker
        preemptible=preemptible_general

        ref=ref
        dict=dict
        amb=amb
        ann=ann
        bwt=bwt
        fai=fai
        pac=pac
        sa=sa
    }

    #Document Tool
    call MarkDuplicates{
      input:
        input_sorted_bam=SamToFastqAllignMerge.bam
        input_sample_name=sample_name
        disk_size=general_disk_size
        mem_size_gb=general_mem_size_gb
        preemptible=preemptible_general
        docker=gatk_docker
        gatk_path=gatk_path

    }

    call ReorderBam {
      input:
      input_bam=MarkDuplicates.bam
      input_sample_name=sample_name

      ref=ref
      dict=dict

      disk_size=general_disk_size
      mem_size_gb=general_mem_size_gb
      preemptible=preemptible
      docker=gatk_docker
      gatk_path=gatk_path

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

    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MarkIlluminaAdapters -I=${input_bam} \
    -O=${input_sample_name}_markilluminaadapters.bam -M=${input_sample_name}_markilluminaadapters_metrics.txt
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
  File input_unmapped_bam

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

    #Piped command containing SamToFastq | bwa mem | samtools view | MergeBamAlignment
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" SamToFastq -I=${input_bam} \
    --FASTQ=/dev/stdout --CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | \
    bwa mem -M -R ${read_group} -p ${ref} /dev/stdin | samtools view -1 | \
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MergeBamAlignment \
    -ALIGNED=/dev/stdin -UNMAPPED=${input_unmapped_bam} \
    -O=${input_sample_name}.sorted.bam -R=${ref}
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk "+ disk_size + " HDD"
  }

  output {
    File bam=${input_sample_name}.sorted.bam
  }
}


task MarkDuplicates {
  File input_sorted_bam
  String input_sample_name

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  command {
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MarkDuplicates \
    -I=${input_sorted_bam} -O=${input_sample_name}.marked_duplicates.bam \
    -M=${input_sample_name}.marked_duplicates.metrics
  }

  output {
    File bam = "${input_sample_name}.marked_duplicates.bam"
    File metrics_duplicates = "${input_sample_name}.marked_duplicates.metrics"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}


task ReorderBam {
  File input_bam
  String input_sample_name

  File ref
  File dict

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  Int cmd_mem_size_gb = mem_size_gb - 1

  command {
    # reorder bam
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" ReorderSam \
    -I=${input_bam} -O=${input_sample_name}.reordered.bam -SD=${ref}

    # then index
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" BuildBamIndex \
    -I=${input_sample_name}.reordered.bam
  }

  output {
    File bam = "${bam_prefix}.reordered.bam"
    File bai = "${bam_prefix}.reordered.bai"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
