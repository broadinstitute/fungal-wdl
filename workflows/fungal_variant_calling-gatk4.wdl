# Variant calling pipeline for fungal haploid genomes
# Developed by Xiao Li (xiaoli@broadinstitute.org) and Aina Martinez Zurita (amartine@broadinstitute.org)
# Fungal Genomics Group, Infectious Disease and Microbiome Program.
# The Broad Institute of MIT and Harvard


workflow GATK4_Germline_Variants {

  # Config Parameters

  String docker
  String gatk_path

  # Reference files


  # Input Files
  Array[String] input_samples
  Array[File] input_uBAMs

  scatter(i in range(length(input_samples))) {
    String sample_name = input_samples[i]
    String input_bam = input_uBAMs[i]

    call MarkIlluminaAdapters {
        input:

    }

    call SamToFastqAllignMerge {
      input:

    }

    call MarkDuplicates{

    }

    call ReorderBam {


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
  File

  String gatk_path
  Int mem_size_gb

  command {
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MarkIlluminaAdapters -I=/gsap/cdcfungal/WGS_pipelines/unaligned_bam_files/TestID-CA03_unaligned_read_pairs.bam -O=TestID-CA03_markilluminaadapters.bam -M=TestID-CA03_markilluminaadapters_metrics.txt

  }

  output {

  }
}
