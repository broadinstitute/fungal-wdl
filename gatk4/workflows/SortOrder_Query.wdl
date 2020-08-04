# Query Sort BAMs
# Helper workflow to query sort input unmapped BAMs into variant calling pipeline for fungal genomes.
# Mantained by Aina Martinez Zurita (amartine@broadinstitute.org)
# Fungal Genomics Group, Infectious Disease and Microbiome Program.
# The Broad Institute of MIT and Harvard

workflow sortSamWorkflow {
	call sortSam
}

task SortSam {
  File input_bam
  String input_sample_name
  Int mem_size_gb
  Int disk_size


  command {
    gatk --java-options "-Xmx${mem_size_gb}G" SortSam -I=${input_bam} \
    -O=${input_sample_name}_querysorted.bam -SO="queryname" \
    --CREATE_INDEX=true
  }

  runtime {
    preemptible: 3
    docker: "amartine/fungal_gatk4:latest"
    memory: mem_size_gb + " GB"
    disks: "local-disk "+ disk_size + " HDD"
  }

  output {
    File SortedBam = "${input_sample_name}_querysorted.bam"
    File SortedBam_index = "${input_sample_name}_querysorted.bai"
  }
}
