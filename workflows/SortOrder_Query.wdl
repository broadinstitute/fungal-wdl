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
