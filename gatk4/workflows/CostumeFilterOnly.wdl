

workflow CostumeFiltration {
  # Config Parameters
  String gatk_docker
  String python2_docker
  Int preemptible_general
  Int general_disk_size
  Int general_mem_size_gb

  String run_name
  File input_vcf

  call MakeReadableVCF {
    input:
      vcf_file=input_vcf,
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
    File costume_filtered_vcf=CostumeVCFFilter.vcf_filtered
    File costume_filtered_stats=CostumeVCFFilter.vcf_filter_stats
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
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
