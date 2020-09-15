# General QC Tasks for BAMs used in variant calling
# Developed by Aina Martinez Zurita (amartine@broadinstitute.org)
# Fungal Genomics Group, Infectious Disease and Microbiome Program.
# The Broad Institute of MIT and Harvard


workflow GATK4_BAM_QC_Tasks {
  # Config Parameters
  String gatk_docker
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
  String sample_name
  File input_bam
  File input_fai


  call CollectMultipleMetrics {
    input:
      input_bam=input_bam,
      input_bam_index=input_fai,
      input_sample_name=sample_name,

      ref_dict=dict,
      ref=ref,
      ref_index=fai,

      preemptible=preemptible_general,
      disk_size=general_disk_size,
      mem_size_gb=general_mem_size_gb,
      docker=gatk_docker,
      gatk_path=gatk_path
  }

  call DepthOfCoverage {
    input:
      input_bam=input_bam,
      input_bam_index=input_fai,
      input_sample_name=sample_name,

      ref_dict=dict,
      ref=ref,
      ref_index=fai,

      preemptible=preemptible_general,
      disk_size=general_disk_size,
      mem_size_gb=general_mem_size_gb,
      docker=gatk_docker,
      gatk_path=gatk_path
  }
}

### Task Definitions ###
task CollectMultipleMetrics {
  File input_bam
  File input_bam_index
  String input_sample_name

  File ref_dict
  File ref
  File ref_index

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  command {
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" CollectMultipleMetrics \
    -R ${ref} -I ${input_bam} -O ${input_sample_name}
  }

  output {
    File alignment_summary_metrics = "${input_sample_name}.alignment_summary_metrics"
    File base_distribution_by_cycle_pdf = "${input_sample_name}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "${input_sample_name}.base_distribution_by_cycle_metrics"
    File insert_size_metrics = "${input_sample_name}.insert_size_metrics"
    File insert_size_histogram = "${input_sample_name}.insert_size_histogram.pdf"
    File quality_by_cycle_metrics = "${input_sample_name}.quality_by_cycle_metrics"
    File quality_by_cycle_pdf = "${input_sample_name}.quality_by_cycle.pdf"
    File quality_distribution_pdf = "${input_sample_name}.quality_distribution.pdf"
    File quality_distribution_metrics = "${input_sample_name}.quality_distribution_metrics"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

}


task DepthOfCoverage {
  File input_bam
  File input_bam_index
  String input_sample_name

  File ref_dict
  File ref
  File ref_index

  Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

  command {
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" DepthOfCoverage \
    -R ${ref} -I ${input_bam} -O ${input_sample_name}
  }

  output {
    File sampleSummary = "${input_sample_name}.sample_summary"
		File sampleStatistics = "${input_sample_name}.sample_statistics"
		File sampleCumulativeCoverageProportions = "${input_sample_name}.sample_cumulative_coverage_proportions"
		File sampleCumulativeCoverageCounts = "${input_sample_name}.sample_cumulative_coverage_counts"
  }

  runtime {
    preemptible: preemptible
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

}
