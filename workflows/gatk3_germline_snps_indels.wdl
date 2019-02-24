# docker file for funpipe
#
#  gatk3_germline_snps_indels.wdl:
#     - standardized pipeline for calling germline SNPs and INDELs.
#     - imports haplotypecaller subworkflow to scatter SNP/INDEL calls
#       during scatter over samples.
#     - developed by Malaria Group, IDMP, Broad Institute.
#


## WORKFLOW DEFINITION
workflow GATK3_Germline_Variants {
    ## config params
    # input data
    String run_name

    File ref                   # path to reference file
    File ref_sa
    File ref_bwt
    File ref_amb
    File ref_ann
    File ref_pac
    File ref_dict
    File ref_index
    File sample_paths_file     # .tsv of (sample_name, sample_sam_path)

    # mem size/ disk size params
    Int small_mem_size_gb
    Int med_mem_size_gb
    Int large_mem_size_gb
    Int extra_large_mem_size_gb

    Int disk_size
    Int med_disk_size
    Int large_disk_size
    Int extra_large_disk_size

    String docker

    String picard_path
    String gatk_path

    # align?
    Boolean do_align

    # hard filtering params: both of these params are required
    String snp_filter_expr
    String indel_filter_expr

    # snpeff
    File organism_gff
    Boolean do_snpeff


    ## task calls
    # run pipeline on each sample, in parallel
    scatter(sample in read_tsv(sample_paths_file)) {
        String sample_name = sample[0]

        if ((length(sample) == 2) && do_align) {
            call SamToFastq {
                input:
                in_bam = sample[1],
                sample_name = sample_name,
                disk_size = large_disk_size,
                mem_size_gb = small_mem_size_gb,
                docker = docker,
                picard_path = picard_path
            }

            call AlignAndSortBAM {
                input:
                sample_name = sample_name,
                fq1 = SamToFastq.fq1,
                fq2 = SamToFastq.fq2,

                ref = ref,
                sa = ref_sa,
                bwt = ref_bwt,
                amb = ref_amb,
                ann = ref_ann,
                pac = ref_pac,
                dict = ref_dict,
                fai = ref_index,

                docker = docker,
                mem_size_gb = small_mem_size_gb,
                disk_size = large_disk_size,
                picard_path = picard_path
            }
        }


        call MarkDuplicates {
            input:
            sample_name = sample_name,
            sorted_bam = select_first([
                AlignAndSortBAM.bam,
                sample[1]]),

            docker = docker,
            picard_path = picard_path,
            mem_size_gb = small_mem_size_gb,
            disk_size = disk_size
        }

        call ReorderBam {
            input:
            bam = MarkDuplicates.bam,

            ref = ref,
            dict = ref_dict,

            docker = docker,
            picard_path = picard_path,
            mem_size_gb = med_mem_size_gb,
            disk_size = disk_size
        }

        call RealignerTargetCreator {
            input:
            sample_name = sample_name,
            in_bam = ReorderBam.out,
            in_bam_index = ReorderBam.out_index,

            ref = ref,
            dict = ref_dict,
            amb = ref_amb,
            sa = ref_sa,
            ann = ref_ann,
            bwt = ref_bwt,
            fai = ref_index,
            pac = ref_pac,

            gatk_path = gatk_path,
            docker = docker,
            mem_size_gb = med_mem_size_gb,
            disk_size = med_disk_size
        }

        call IndelRealigner {
            input:
            in_bam = ReorderBam.out,
            in_bam_index = ReorderBam.out_index,
            sample_name = sample_name,
            intervals = RealignerTargetCreator.intervals,

            ref = ref,
            dict = ref_dict,
            amb = ref_amb,
            sa = ref_sa,
            ann = ref_ann,
            bwt = ref_bwt,
            fai = ref_index,
            pac = ref_pac,

            gatk_path = gatk_path,
            mem_size_gb = med_mem_size_gb,
            docker = docker,
            disk_size = med_disk_size

        }

        call HaplotypeCaller {
            input:
            input_bam = IndelRealigner.bam,
            input_bam_index = IndelRealigner.index,
            sample_name = sample_name,

            gvcf_name = "${sample_name}.g.vcf",
            gvcf_index = "${sample_name}.g.vcf.idx",

            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index,

            mem_size_gb = med_mem_size_gb,
            disk_size = med_disk_size,
            docker = docker,
            gatk_path = gatk_path
        }
    }

    call CombineGVCFs {
        input:
        vcf_files = HaplotypeCaller.output_gvcf,
        vcf_index_files = HaplotypeCaller.output_gvcf_index,

        ref = ref,
        ref_dict = ref_dict,
        ref_index = ref_index,

        docker = docker,
        gatk_path = gatk_path,
        mem_size_gb = med_mem_size_gb,
        disk_size = disk_size
    }

    call GenotypeGVCFs {
        input:
        vcf_file = CombineGVCFs.out,
        vcf_index_file = CombineGVCFs.out_index,

        ref = ref,
        ref_dict = ref_dict,
        ref_index = ref_index,

        docker = docker,
        gatk_path = gatk_path,
        mem_size_gb = med_mem_size_gb,
        disk_size = disk_size
    }

    call HardFiltration {
        input:
        vcf = GenotypeGVCFs.output_vcf_name,
        vcf_index = GenotypeGVCFs.output_vcf_index_name,
        snp_filter_expr = snp_filter_expr,
        indel_filter_expr = indel_filter_expr,

        ref = ref,
        ref_dict = ref_dict,
        ref_index = ref_index,

        output_filename = "${run_name}.hard_filtered.vcf.gz",

        docker = docker,
        gatk_path = gatk_path,
        mem_size_gb = extra_large_mem_size_gb,
        disk_size = extra_large_disk_size
    }

    # to do: genotype filtering step
    if (do_snpeff) {
        call SnpEff {
            input:
            ref = ref,
            organism_gff = organism_gff,
            output_vcf_name = "${run_name}.snpeff.vcf.gz",
            vcf = HardFiltration.out,

            mem_size_gb = extra_large_mem_size_gb,
            disk_size = extra_large_disk_size
        }
    }

    output {
        File gvcf = select_first([SnpEff.out, HardFiltration.out])
    }
}


## TASK DEFINITIONS
task SamToFastq {
    File in_bam
    String sample_name

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    String out_fq1 = "${sample_name}.1.fq"
    String out_fq2 = "${sample_name}.2.fq"
    command {
        java -Xmx${mem_size_gb}G -jar ${picard_path} SamToFastq INPUT=${in_bam} FASTQ=${out_fq1} SECOND_END_FASTQ=${out_fq2} VALIDATION_STRINGENCY=LENIENT
    }
    output {
        String done = "Done"
        File fq1 = out_fq1
        File fq2 = out_fq2
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        in_bam: "The bam file to convert to fastq."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        out_fq1: "The fastq file containing the first read of each pair."
        out_fq2: "The fastq file containing the second read of each pair"
    }
}


task AlignAndSortBAM {
    String sample_name

    File fq1
    File fq2

    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    String read_group = "'@RG\\tID:FLOWCELL_${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:LIB_${sample_name}'"
    command {
        bwa mem -R ${read_group} ${ref} ${fq1} ${fq2} | samtools view -bS -> ${sample_name}.aligned.bam
        java -Xmx${mem_size_gb}G -jar ${picard_path} SortSam I=${sample_name}.aligned.bam O=${sample_name}.sorted.bam SO=coordinate
    }

    output {
        File bam = "${sample_name}.sorted.bam"
    }

    runtime {
        task_name: "AlignBAM"
        preemptible: 5
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk "+ disk_size + " HDD"
    }

    parameter_meta {
        ref: "fasta file of reference genome"
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        fq_array: "An array containing the paths to the first and second fastq files."
        read_group: "The read group string that will be included in the bam header."
    }
}


# mark duplicate reads in bam
task MarkDuplicates {
    File sorted_bam
    String sample_name

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${mem_size_gb}G -jar ${picard_path} MarkDuplicates \
            I=${sorted_bam} \
            O=${sample_name}.marked_duplicates.bam \
            M=${sample_name}.marked_duplicates.metrics
    }

    output {
        File bam = "${sample_name}.marked_duplicates.bam"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# reorder and index a bam
task ReorderBam {
    File ref
    File dict
    File bam
    String bam_prefix = basename(bam, '.bam')

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        # reorder bam
        java -Xmx${cmd_mem_size_gb}G -jar ${picard_path} ReorderSam \
            I=${bam} \
            O=${bam_prefix}.reordered.bam \
            R=${ref}

        # then index
        java -Xmx${cmd_mem_size_gb}G -jar ${picard_path} BuildBamIndex \
            I=${bam_prefix}.reordered.bam
    }

    output {
        File out = "${bam_prefix}.reordered.bam"
        File out_index = "${bam_prefix}.reordered.bai"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# merge vcfs before genotyping
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
    String docker
    String gatk_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T CombineGVCFs \
            -R ${ref} \
            -o ${gvcf_out} \
            --variant ${sep=" --variant " vcf_files}
    }
    output {
       File out = gvcf_out
       File out_index = gvcf_out_index
    }

    runtime {
        preemptible: 4
        docker:docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# genotype gvcfs
task GenotypeGVCFs {
    File ref
    File ref_dict
    File ref_index
    File vcf_file
    File vcf_index_file
    String gvcf_out = "genotyped_gvcfs.vcf.gz"

    Int disk_size
    Int mem_size_gb
    String docker
    String gatk_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T GenotypeGVCFs \
            -R ${ref} \
            -o ${gvcf_out} \
            --variant ${vcf_file} \
    }
    output {
        File output_vcf_name = gvcf_out
        File output_vcf_index_name = "${gvcf_out}.tbi"
    }

    runtime {
        preemptible: 4
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


task HardFiltration {
    # hard-filter a vcf, if vqsr not available
    # http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
    File ref
    File vcf
    File vcf_index
    File ref_dict
    File ref_index
    String output_filename

    String snp_filter_expr
    String indel_filter_expr

    Int disk_size
    Int mem_size_gb
    String docker
    String gatk_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        # select snps
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.g.vcf

        # filter snps
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.g.vcf

        # select indels
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path}\
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType INDEL \
            -o raw_indels.g.vcf

        # filter indels
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.g.vcf

        # combine variants
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path}\
            -T CombineVariants \
            -R ${ref} \
            --variant filtered_snps.g.vcf \
            --variant filtered_indels.g.vcf \
            -o ${output_filename} \
            --genotypemergeoption UNSORTED
    }

    output {
        File out = "${output_filename}"
    }

    runtime {
        preemptible: 3
        docker:docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


task HaplotypeCaller {
  File input_bam
  File input_bam_index

  String gvcf_name
  String gvcf_index

  File ref_dict
  File ref
  File ref_index
  String sample_name

  Int disk_size
  Int mem_size_gb
  Int cmd_mem_size_gb = mem_size_gb - 1

  String docker
  String gatk_path

  String out = "${sample_name}.g.vcf"

  command {
    java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
      -T HaplotypeCaller \
      -R ${ref} \
      -I ${input_bam} \
      -o ${gvcf_name} \
      -ERC "GVCF" \
      -ploidy 1 \
      -variant_index_type LINEAR \
      -variant_index_parameter 128000 \
      --read_filter OverclippedRead
  }

  output {
      #To track additional outputs from your task, please manually add them below
      File output_gvcf = "${gvcf_name}"
      File output_gvcf_index = "${gvcf_index}"
  }

  runtime {
    task_name: "HaplotypeCaller"
    preemptible: 3
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  parameter_meta {
      gatk: "Executable jar for the GenomeAnalysisTK"
      ref: "fasta file of reference genome"
      sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
      sample_dir: "The sample-specific directory inside output_dir for each sample."
      in_bam: "The bam file to call HaplotypeCaller on."
      out: "VCF file produced by haplotype caller."
  }
}


task SnpEff {
    # annotate variants
    # Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
    File vcf
    File ref
    File organism_gff
    String organism_name
    String output_vcf_name

    Int disk_size
    Int mem_size_gb
    String snpeff_path = "/opt/snpEff/snpEff.jar"
    String docker = "maxulysse/snpeff:1.3"
    String organism_db_dir = "/opt/snpEff/data/" + organism_name + "/"
    String snpeff_config_path = "/opt/snpEff/snpEff.config"

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        # init database
        echo "${organism_name}.genome : ${organism_name}" >> ${snpeff_config_path}
        mkdir -p ${organism_db_dir}
        mv ${ref} ${organism_db_dir}/sequences.fa
        mv ${organism_gff} ${organism_db_dir}/genes.gff

        # build db
        java -jar ${snpeff_path} build -gff3 -v ${organism_name}

        # run snpeff
        java -Xmx${cmd_mem_size_gb}G -jar ${snpeff_path} \
            -config ${snpeff_config_path} \
            -formatEff -no-downstream -no-intergenic \
            -no-upstream -no-utr -noStats \
            -treatAllAsProteinCoding false \
            ${organism_name} ${vcf} > ${output_vcf_name}

        # gzip result
        gzip ${output_vcf_name}
    }

    output {
        File out = "${output_vcf_name}.gz"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


task RealignerTargetCreator {
    String gatk_path
    String docker
    Int mem_size_gb
    Int disk_size

    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa

    File in_bam
    File in_bam_index

    String sample_name
    String out = "${sample_name}.interval_list"
    command {
        java -Xmx${mem_size_gb}G -jar ${gatk_path} -T RealignerTargetCreator -R ${ref} -I ${in_bam} -o ${out}
    }
    output {
        File intervals = out
    }
    runtime {
        preemptible: 5
        memory: mem_size_gb
        docker: docker
        disks: "local-disk " + disk_size + " HDD"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for the gatk task"
        out: "The intervals list to be used by IndelRealigner"
    }
}


task IndelRealigner {
    String gatk_path
    String docker
    Int mem_size_gb
    Int disk_size

    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa

    File in_bam
    File in_bam_index
    File intervals
    String sample_name
    String out = "${sample_name}.indels_realigned.bam"


    command {
        samtools index ${in_bam}
        java -Xmx${mem_size_gb}G -jar ${gatk_path} -T IndelRealigner -R ${ref} -I ${in_bam} -targetIntervals ${intervals} -o ${out}
        samtools index ${out}
    }
    output {
        File bam = "${sample_name}.indels_realigned.bam"
        File index = "${sample_name}.indels_realigned.bam.bai"
    }
    runtime {
        preemptible: 5
        memory: mem_size_gb
        docker: docker
        disks: "local-disk " + disk_size + " HDD"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for the gatk task"
        intervals: "The intervals list to be used by IndelRealigner"
        out: "the bam including realigned indels."
    }
}
