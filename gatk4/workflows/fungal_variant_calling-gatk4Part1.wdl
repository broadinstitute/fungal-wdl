# Variant calling pipeline for fungal haploid genomes: Part I Single BAM
# Developed by Xiao Li (xiaoli@broadinstitute.org) and Aina Martinez Zurita (amartine@broadinstitute.org)
# Fungal Genomics Group, Infectious Disease and Microbiome Program.
# The Broad Institute of MIT and Harvard


workflow GATK4_BAM_VariantCalling_Workflow {
  # Config Parameters
  String gatk_docker
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
  String sample_name
  File input_bam

  #Document Tool
  call MarkIlluminaAdapters {
      input:
        input_bam=input_bam,
        input_sample_name=sample_name,
        gatk_path=gatk_path,
        docker=gatk_docker,
        preemptible=preemptible_general,
        disk_size=general_disk_size,
        mem_size_gb=general_mem_size_gb
  }

  #Document tool
  call SamToFastqAllignMerge {
    input:
      input_bam=MarkIlluminaAdapters.bam,
      input_sample_name=sample_name,
      input_unmapped_bam=input_bam,
      gatk_path=gatk_path,
      docker=gatk_docker,
      preemptible=preemptible_general,

      ref=ref,
      dict=dict,
      amb=amb,
      ann=ann,
      bwt=bwt,
      fai=fai,
      pac=pac,
      sa=sa
  }

  #Document Tool
  call MarkDuplicates{
    input:
      input_sorted_bam=SamToFastqAllignMerge.bam,
      input_sample_name=sample_name,
      disk_size=general_disk_size,
      mem_size_gb=general_mem_size_gb,
      preemptible=preemptible_general,
      docker=gatk_docker,
      gatk_path=gatk_path

  }

  call ReorderBam {
    input:
      input_bam=MarkDuplicates.bam,
      input_sample_name=sample_name,

      ref=ref,
      dict=dict,

      disk_size=general_disk_size,
      mem_size_gb=general_mem_size_gb,
      preemptible=preemptible_general,
      docker=gatk_docker,
      gatk_path=gatk_path
  }

  call HaplotypeCaller {
    input:
      input_bam=ReorderBam.bam,
      input_bam_index=ReorderBam.bai,
      input_sample_name=sample_name,

      ref_dict=dict,
      ref=ref,
      ref_index=fai,

      preemptible=preemptible_HC,
      docker=gatk_docker,
      gatk_path=gatk_path
  }

  output {
    File processed_bam=ReorderBam.bam
    File processed_bai=ReorderBam.bai

    File gvcf=HaplotypeCaller.output_gvcf
    File gvcf_index=HaplotypeCaller.output_gvcf_index
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

    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MarkIlluminaAdapters -I ${input_bam} \
    -O ${input_sample_name}_markilluminaadapters.bam -M ${input_sample_name}_markilluminaadapters_metrics.txt
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

  String read_group = "'@RG\\tID:FLOWCELL_${input_sample_name}\\tSM:${input_sample_name}\\tPL:ILLUMINA\\tLB:LIB_${input_sample_name}'"

  command {
    # e - exit when a command fails
    # u - exit when script tries to use undeclared variables
    # x - trace what gets executed
    # o - exit status of the last command that threw a non-zero exit code is returned
    set -euxo pipefail

    #We need to assign RG to all the unmapped reads
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" AddOrReplaceReadGroups \
    -I ${input_unmapped_bam} -O ${input_sample_name}.readgroups_unmapped.bam -ID FLOWCELL_${input_sample_name} \
    -LB LIB_${input_sample_name} -PL ILLUMINA -SM ${input_sample_name} \
    -PU unit1

    #Piped command containing SamToFastq | bwa mem | samtools view | MergeBamAlignment
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" SamToFastq -I ${input_bam} \
    --FASTQ /dev/stdout  -CLIP_ATTR XT  -CLIP_ACT 2  -INTER true -NON_PF true | \
    bwa mem -M -R ${read_group} -p ${ref} /dev/stdin | samtools view -1 | \
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MergeBamAlignment \
    -ALIGNED /dev/stdin -UNMAPPED ${input_sample_name}.readgroups_unmapped.bam \
    -O ${input_sample_name}.sorted.bam -R ${ref}

  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk "+ disk_size + " HDD"
  }

  output {
    File bam="${input_sample_name}.sorted.bam"
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
    -I ${input_sorted_bam} -O ${input_sample_name}.marked_duplicates.bam \
    -M ${input_sample_name}.marked_duplicates.metrics
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
    -I ${input_bam} -O ${input_sample_name}.reordered.bam -SD ${ref}

    # then index
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" BuildBamIndex \
    -I ${input_sample_name}.reordered.bam
  }

  output {
    File bam = "${input_sample_name}.reordered.bam"
    File bai = "${input_sample_name}.reordered.bai"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}

task HaplotypeCaller {
  File input_bam
  File input_bam_index
  String input_sample_name
  String gvcf_name = "${input_sample_name}.g.vcf.gz"
  String gvcf_index = "${input_sample_name}.g.vcf.gz.tbi"

  File ref_dict
  File ref
  File ref_index

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  command {
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" HaplotypeCaller \
    -R ${ref} -I ${input_bam} -O ${gvcf_name} -ERC GVCF -ploidy 1
  }

  output {
    File output_gvcf = "${gvcf_name}"
    File output_gvcf_index = "${gvcf_index}"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

}
