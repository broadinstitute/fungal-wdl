# Variant calling pipeline for fungal haploid genomes: Part II Joint Genotyping
# Developed by Xiao Li (xiaoli@broadinstitute.org) and Aina Martinez Zurita (amartine@broadinstitute.org)
# Fungal Genomics Group, Infectious Disease and Microbiome Program.
# The Broad Institute of MIT and Harvard



workflow GATK4_JointGenotyping_Workflow {
  # Config Parameters
  String gatk_docker
  String python2_docker
  String gatk_path
  Int preemptible_general
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
  String run_name
  Array[File] gvcf_files
  Array[File] gvcf_index_files

  call CombineGVCFs {
    input:
      vcf_files = gvcf_files,
      vcf_index_files = gvcf_index_files,

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

  call MakeReadableVCF {
    input:
      vcf_file=HardFiltration.all_filtered_variants,
      output_filename=run_name,
      docker = gatk_docker,
      mem_size_gb = general_mem_size_gb,
      disk_size = general_disk_size,
      preemptible = preemptible_general
  }

  call CostumeVCFFilter{
    input:
      vcf_file=MakeReadableVCF.vcf_output,
      output_filename=run_name,
      docker = python2_docker,
      mem_size_gb = general_mem_size_gb,
      disk_size = general_disk_size,
      preemptible = preemptible_general
  }

  output {
    File filtered_vcf=HardFiltration.all_filtered_variants
    File filtered_snps=HardFiltration.snps
    File filtered_indels=HardFiltration.indels
    File costume_filtered_vcf=CostumeVCFFilter.vcf_filtered
    File costume_filtered_stats=CostumeVCFFilter.vcf_filter_stats
    File multialelic_stats=CostumeVCFFilter.multialelic_stats
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
    -filter "QD < 20.0" --filter-name "QD20" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -O snps_filtered.vcf.gz

    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" VariantFiltration \
    -V indels.vcf.gz \
    -R ${ref} \
    -filter "QD < 20.0" --filter-name "QD20" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -O indels_filtered.vcf.gz

    #Combine output vcfs into a single output. The output of VariantFiltration
    # is sorted.
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" MergeVcfs \
    -R ${ref} \
    -I snps_filtered.vcf.gz \
    -I indels_filtered.vcf.gz \
    -O ${output_filename}

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
task MakeReadableVCF {
  File vcf_file
  String output_filename

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker

  command {
    bcftools view ${vcf_file} > ${output_filename}.nonbinary.vcf
  }

  output {
    File vcf_output= "${output_filename}.nonbinary.vcf"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}

task CostumeVCFFilter {
  File vcf_file
  String output_filename

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker

  command {
    python /broad-fungalgroup/scripts/SNPs/filterGatkGenotypes.py --min_GQ 50 \
    --min_percent_alt_in_AD 0.8 --min_total_DP 10 ${vcf_file} \
     > ${output_filename}.filtered_SNPs_GQ50_AD08_DP10.vcf 2> ${output_filename}.variant_qc_genotype_filter.tsv
  }

  output {
    File vcf_filtered= "${output_filename}.filtered_SNPs_GQ50_AD08_DP10.vcf"
    File vcf_filter_stats= "${output_filename}.variant_qc_genotype_filter.tsv"
    File multialelic_stats= "${vcf_file}.multiallelic.txt"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
