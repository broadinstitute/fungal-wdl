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
  String run_name

  scatter(i in range(length(input_samples))) {
    String sample_name = input_samples[i]
    String input_bam = input_uBAMs[i]

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
  }

  call CombineGVCFs {
    input:
      vcf_files = HaplotypeCaller.output_gvcf,
      vcf_index_files = HaplotypeCaller.output_gvcf_index,

      ref = ref,
      ref_dict = dict,
      ref_index = fai,

      docker = gatk_docker,
      gatk_path = gatk_path,
      mem_size_gb = general_mem_size_gb,
      disk_size = general_disk_size,
      preemptible = preemptible_general

  }

  call GenotypeGVCFs {
    input:
      vcf_file = CombineGVCFs.out,
      vcf_index_file = CombineGVCFs.out_index,

      ref = ref,
      ref_dict = dict,
      ref_index = fai,

      docker = gatk_docker,
      gatk_path = gatk_path,
      mem_size_gb = general_mem_size_gb,
      disk_size = general_disk_size,
      preemptible = preemptible_general
  }

  call HardFiltration {
    input:
    ref = ref,
    ref_dict = dict,
    ref_index = fai,

    vcf = GenotypeGVCFs.output_vcf_name,
    vcf_index = GenotypeGVCFs.output_vcf_index_name,
    output_filename = "${run_name}.hard_filtered.vcf.gz",

    preemptible = preemptible_general,
    docker = gatk_docker,
    gatk_path = gatk_path

  }


  #call CostumFiltration {
  #
  #
  #}


  output {
    File filtered_vcf=HardFiltration.all_filtered_variants
    File filtered_snps=HardFiltration.snps
    File filtered_indels=HardFiltration.indels

    Array[File] processed_bams=ReorderBam.bam
    Array[File] processed_bais=ReorderBam.bai
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

  String read_group = "'@RG\\tID:FLOWCELL_${input_sample_name}\\tSM:${input_sample_name}\\tPL:ILLUMINA\\tLB:LIB_${input_sample_name}'"

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


task CombineGVCFs {
  File ref
  File ref_dict
  File ref_index
  Array[File] vcf_files
  Array[File] vcf_index_files

  String gvcf_out = "combined_gvcfs.vcf.gz"
  String gvcf_out_index = "combined_gvcfs.vcf.gz.tbi"

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  Int cmd_mem_size_gb = mem_size_gb - 1

  command {
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" CombineGVCFs \
    -R ${ref} -O ${gvcf_out} --variant ${sep=" --variant " vcf_files}
  }

 output {
   File out = gvcf_out
   File out_index = gvcf_out_index
 }

 runtime {
   preemptible: preemptible
   docker:docker
   memory: mem_size_gb + " GB"
   disks: "local-disk " + disk_size + " HDD"
 }
}


task GenotypeGVCFs {
  File ref
  File ref_dict
  File ref_index
  File vcf_file
  File vcf_index_file
  String gvcf_out = "genotyped_gvcfs.vcf.gz"

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  Int cmd_mem_size_gb = mem_size_gb - 1

  command {
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" GenotypeGVCFs \
    -R ${ref} -O ${gvcf_out} -V ${vcf_file}
  }

  output {
    File output_vcf_name = gvcf_out
    File output_vcf_index_name = "${gvcf_out}.tbi"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}


task HardFiltration {
  File ref
  File vcf
  File vcf_index
  File ref_dict
  File ref_index
  String output_filename

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  Int cmd_mem_size_gb = mem_size_gb - 1

  command {
    #Select variants, both SNPs and raw_indels. Mixed sites are discarted.
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" SelectVariants \
    -V ${vcf}  -R ${ref} -select-type SNP -O snps.vcf.gz

    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" SelectVariants \
    -V ${vcf}  -R ${ref} -select-type INDEL -O indels.vcf.gz

    #Filter variants, different paramters for Indels and SNPs. Recommendations by GATK.
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" VariantFiltration \
    -V snps.vcf.gz \
    -R ${ref} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -O snps_filtered.vcf.gz

    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" VariantFiltration \
    -V indels.vcf.gz \
    -R ${ref} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -O indels_filtered.vcf.gz

    #Combine output vcfs into a single output. The output of VariantFiltration
    # is sorted.
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" MergeVcfs \
    -R ${ref} \
    -I=snps_filtered.vcf.gz \
    -I=indels_filtered.vcf.gz \
    -O=${output_filename}

  }

  output {
    File indels="indels_filtered.vcf.gz"
    File snps="snps_filtered.vcf.gz"
    File all_filtered_variants="${output_filename}"

  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"

  }
}

#CostumFiltration task
