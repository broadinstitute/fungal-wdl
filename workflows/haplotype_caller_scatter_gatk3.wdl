## Copyright Broad Institute, 2017
##
## This WDL workflow runs HaplotypeCaller in GVCF mode on a single sample according to
## the GATK Best Practices (June 2016) scattered across intervals.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##
## Cromwell version support
## - Successfully tested on v28
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.


workflow HaplotypeCallerGvcf_GATK3 {
  # WORKFLOW DEFINITION
  File input_bam
  File input_bam_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File scattered_calling_intervals_list

  String gatk_path
  String picard_path
  String gatk_docker
  String picard_docker

  Int merge_gvcfs_mem_size_gb
  Int haplotypecaller_mem_size_gb
  Int merge_gvcfs_disk_size
  Int haplotypecaller_disk_size

  Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

  String sample_basename = basename(input_bam, ".bam")

  String gvcf_name = sample_basename + ".g.vcf.gz"
  String gvcf_index = sample_basename + ".g.vcf.gz.tbi"
  String genotyped_gvcf_name = sample_basename + ".genotyped.g.vcf.gz"

  # Call variants in parallel over grouped calling intervals
  scatter (interval_file in scattered_calling_intervals) {
    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        interval_list = interval_file,
        gvcf_name = gvcf_name,
        gvcf_index = gvcf_index,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
	      docker_image = gatk_docker,
	      gatk_path = gatk_path,
        mem_size_gb = haplotypecaller_mem_size_gb,
        disk_size = haplotypecaller_disk_size
    }
  }

  # Merge per-interval GVCFs
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indices = HaplotypeCaller.output_gvcf_index,
      vcf_name = gvcf_name,
      vcf_index = gvcf_index,
      docker_image = picard_docker,
      picard_path = picard_path,
      mem_size_gb = merge_gvcfs_mem_size_gb,
      disk_size = merge_gvcfs_disk_size
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_gvcf = MergeVCFs.output_vcf
    File output_gvcf_index = MergeVCFs.output_vcf_index
  }
}


task HaplotypeCaller {
  # TASK DEFINITIONS
  # HaplotypeCaller per-sample in GVCF mode
  # ${default="1" ploidy}
  File input_bam
  File input_bam_index
  String gvcf_name
  String gvcf_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File interval_list
  Int? interval_padding
  Int? contamination
  Int? max_alt_alleles
  Int? ploidy

  Int disk_size
  Int mem_size_gb
  Int cmd_mem_size_gb = mem_size_gb - 1

  String docker_image
  String gatk_path

  command {
    java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -o ${gvcf_name} \
      -L ${interval_list} \
      -ip ${default=100 interval_padding} \
      -contamination ${default=0 contamination} \
      --read_filter OverclippedRead \
      -ERC GVCF \
      -ploidy 1
  }

  runtime {
    preemptible: 3
    docker: docker_image
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_gvcf = "${gvcf_name}"
    File output_gvcf_index = "${gvcf_index}"
  }
}


task MergeVCFs {
  # Merge GVCFs generated per-interval for the same sample
  Array[File] input_vcfs
  Array[File] input_vcfs_indices
  String vcf_name
  String vcf_index

  Int disk_size
  Int mem_size_gb
  Int cmd_mem_size_gb = mem_size_gb - 1
  Int compression_level = 5

  String docker_image
  String picard_path

  command {
    java -Dsamjdk.compression_level=${compression_level} \
      -Xmx${cmd_mem_size_gb}G -jar \
      ${picard_path} \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${vcf_name}
  }

  runtime {
  	preemptible: 3
    docker: docker_image
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
}

  output {
    File output_vcf = "${vcf_name}"
    File output_vcf_index = "${vcf_index}"
  }
}
