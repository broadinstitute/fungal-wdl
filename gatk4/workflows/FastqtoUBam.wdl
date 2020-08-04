workflow FastqToSamTest {
	call FastqToSam
}

task FastqToSam {
	File inputFastq1
	File inputFastq2
	String sampleName

	Int disk_size
  Int mem_size_gb
  Int preemptible
  String docker
  String gatk_path

	command <<<
		${gatk_path} --java-options "-Xmx${mem_size_gb}G" FastqToSam \
		--FASTQ=${inputFastq1} \
		--FASTQ2=${inputFastq2} \
		--OUTPUT=${sampleName}.unaligned_new.bam \
		--READ_GROUP_NAME=1 \
		--SAMPLE_NAME=${sampleName} \
		--LIBRARY_NAME=lib1 \
		--PLATFORM_UNIT=unit1 \
		--PLATFORM=illumina
	>>>

	output {
		File Bam = "${sampleName}.unaligned_new.bam"
	}

	runtime {
		preemptible: preemptible
    docker: docker
		memory: mem_size_gb + " GB"
		disks: "local-disk " + disk_size + " HDD"

	}
}
