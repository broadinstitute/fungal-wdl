#
#  gatk3_germline_snps_indels.wdl:
#     - standardized pipeline for calling germline SNPs and INDELs.
#     - imports haplotypecaller subworkflow to scatter SNP/INDEL calls over intervals,
#       during scatter over samples.
#     - developed by Malaria Group, IDMP, Broad Institute.
#


## WORKFLOW DEFINITION
workflow GATK3_Germline_Variants {
    ## config params
    # input data
    File ref                   # path to reference file
    File ref_sa
    File ref_bwt
    File ref_amb
    File ref_ann
    File ref_pac
    File ref_dict
    File ref_index
    String run_name
    File sample_paths_file     # .tsv of (sample_name, sample_sam_path)
    File interval_files_list   # intervals files to be scattered over during HC
    File interval_list         # all intervals

    # mem size/ disk size params
    Int small_mem_size_gb
    Int med_mem_size_gb
    Int large_mem_size_gb
    Int extra_large_mem_size_gb
    Int disk_size
    Int med_disk_size
    Int large_disk_size
    Int extra_large_disk_size

    # gatk/ picard/ genomes-in-the-cloud
    String gatk_docker         # gatk3 docker (e.g. broadinstitute/gatk3:3.8-0)
    String gatk_path_to_gatk
    String gatk4_docker
    String gatk4_path_to_gatk
    String gitc_docker         # genomes-in-the-cloud (gitc) docker
    String gitc_path_to_picard
    String gitc_path_to_gatk
    String gitc_path_to_bwa
    String gitc_path_to_samtools

    # align?
    Boolean do_align

    # variant quality control param
    # either "vqsr" or "hard_filtering"
    String variant_qc

    # hard filtering params
    # if variant_qc == "hard_filtering"
    # both of these params are required
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
            call AlignBwaMem {
                input:
                ref = ref,
                ref_sa = ref_sa,
                ref_bwt = ref_bwt,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_pac = ref_pac,
                ref_dict = ref_dict,
                ref_index = ref_index,
                sample_name = sample_name,
                bam_or_fastq = sample[1],
                gitc_docker = gitc_docker,
                samtools_path_gitc = gitc_path_to_samtools,
                picard_path_gitc = gitc_path_to_picard,
                bwa_path_gitc = gitc_path_to_bwa,
                mem_size_gb = small_mem_size_gb,
                disk_size = large_disk_size
            }
        }

        if (length(sample) == 3){
            call AlignBwaMemR1R2 {
                input:
                ref = ref,
                ref_sa = ref_sa,
                ref_bwt = ref_bwt,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_pac = ref_pac,
                ref_dict = ref_dict,
                ref_index = ref_index,
                sample_name = sample_name,
                fastq_R1 = sample[1],
                fastq_R2 = sample[2],
                gitc_docker = gitc_docker,
                samtools_path_gitc = gitc_path_to_samtools,
                picard_path_gitc = gitc_path_to_picard,
                bwa_path_gitc = gitc_path_to_bwa,
                mem_size_gb = small_mem_size_gb,
                disk_size = large_disk_size
            }
        }

        call MarkDuplicates {
            input:
            sample_name = sample_name,
            sorted_bam = select_first([AlignBwaMem.aligned_bam, AlignBwaMemR1R2.aligned_bam, sample[1]]),
            picard_docker = gitc_docker,
            picard_path = gitc_path_to_picard,
            mem_size_gb = small_mem_size_gb,
            disk_size = disk_size
        }

        call ReorderBam {
            input:
            ref = ref,
            dict = ref_dict,
            bam = MarkDuplicates.bam,
            picard_docker = gitc_docker,
            picard_path = gitc_path_to_picard,
            mem_size_gb = med_mem_size_gb,
            disk_size = disk_size
        }

        call HaplotypeCaller {
            input:
            input_bam = ReorderBam.out,
            input_bam_index = ReorderBam.out_index,

            gvcf_name = "${sample_name}.g.vcf",
            gvcf_index = "${sample_name}.g.vcf.idx",

            ref_fasta = ref,
            ref_dict = ref_dict,
            ref_fasta_index = ref_index,

            mem_size_gb = med_mem_size_gb,
            disk_size = med_disk_size,

            gatk_docker = gatk_docker,
            gatk_path = gatk_path_to_gatk
        }
    }

    call CombineGVCFs {
        input:
        ref = ref,
        ref_dict = ref_dict,
        ref_index = ref_index,
        interval = interval,
        vcf_files = HaplotypeCaller.output_gvcf,
        vcf_index_files = HaplotypeCaller.output_gvcf_index,
        gatk4_docker = gatk4_docker,
        gatk4_path = gatk4_path_to_gatk,
        mem_size_gb = med_mem_size_gb,
        disk_size = disk_size
    }

    call GenotypeGVCFs {
        input:
        ref = ref,
        ref_dict = ref_dict,
        ref_index = ref_index,
        vcf_file = CombineGVCFs.out,
        vcf_index_file = CombineGVCFs.out_index,
        gatk4_docker = gatk4_docker,
        gatk4_path = gatk4_path_to_gatk,
        mem_size_gb = med_mem_size_gb,
        disk_size = disk_size
    }

    if (variant_qc == "hard_filtering") {
        call HardFiltration {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index,
            vcf = GenotypeGVCFs.output_vcf_name,
            vcf_index = GenotypeGVCFs.output_vcf_index_name,
            snp_filter_expr = snp_filter_expr,
            indel_filter_expr = indel_filter_expr,
            output_filename = "${run_name}.hard_filtered.g.vcf.gz",
            gatk_docker = gatk_docker,
            gatk_path = gatk_path_to_gatk,
            mem_size_gb = extra_large_mem_size_gb,
            disk_size = extra_large_disk_size
        }
    }

    # to do: genotype filtering step
    if (do_snpeff) {
        call SnpEff {
            input:
            ref = ref,
            organism_gff = organism_gff,
            output_vcf_name = "${run_name}.snpeff.g.vcf",
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
task AlignBwaMem {
    File ref
    File ref_sa
    File ref_bwt
    File ref_amb
    File ref_ann
    File ref_pac
    File ref_dict
    File ref_index

    File bam_or_fastq
    String sample_name

    Int disk_size
    Int mem_size_gb
    String gitc_docker
    String bwa_path_gitc
    String picard_path_gitc
    String samtools_path_gitc

    String d = "$"

    command <<<
        if [[ ${bam_or_fastq} = *".bam" ]] ; then
            # bam to fastq
            ${samtools_path_gitc} fastq ${bam_or_fastq} > ${sample_name}.fastq
            # bwa mem
            ${bwa_path_gitc} mem -M -R "$(${samtools_path_gitc} view -H ${bam_or_fastq} | grep @RG | sed 's/\t/\\t/g')" -p ${ref} ${sample_name}.fastq > ${sample_name}.sam
            # sam2bam, sort and index
            ${samtools_path_gitc} view ${sample_name}.sam -bS > ${sample_name}.bam
        else
            # bwa mem
            ${bwa_path_gitc} mem -M -p ${ref} ${bam_or_fastq} > ${sample_name}.sam
            # add read group
            if [[ ${bam_or_fastq} = *".gz" ]] ; then
                IFS=':' read -r -a read_group_data <<< "$(zcat ${bam_or_fastq} | head -n 1)"
            else
                IFS=':' read -r -a read_group_data <<< "$(cat ${bam_or_fastq} | head -n 1)"
            fi
            IFS=" " read_group_data=(${d}{read_group_data[@]})
            java -jar ${picard_path_gitc} AddOrReplaceReadGroups \
                I=${sample_name}.sam \
                O=${sample_name}.bam \
                RGID=${d}{read_group_data[0]:1}.${d}{read_group_data[1]} \
                RGPU=${d}{read_group_data[2]}.${d}{read_group_data[3]} \
                RGLB=${sample_name} \
                RGSM=${sample_name} \
                RGPL=ILLUMINA
        fi
        ${samtools_path_gitc} sort -o ${sample_name}.sorted.bam ${sample_name}.bam
        ${samtools_path_gitc} index ${sample_name}.sorted.bam
    >>>

    output {
        File aligned_bam = "${sample_name}.sorted.bam"
        File aligned_bam_index = "${sample_name}.sorted.bam.bai"
    }

    runtime {
        preemptible: 5
        docker: gitc_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task AlignBwaMemR1R2 {
    File ref
    File ref_sa
    File ref_bwt
    File ref_amb
    File ref_ann
    File ref_pac
    File ref_dict
    File ref_index

    File fastq_R1
    File fastq_R2
    String sample_name

    Int disk_size
    Int mem_size_gb
    String gitc_docker
    String bwa_path_gitc
    String samtools_path_gitc
    String picard_path_gitc

    String d = "$"

    command <<<
        # bwa mem
        ${bwa_path_gitc} mem -M ${ref} ${fastq_R1} ${fastq_R2} > ${sample_name}.sam
        # add read groups, convert to bam
        fastq_file=${fastq_R1}
        if [ "${d}{fastq_file##*.}" == "gz" ] ; then
            IFS=':' read -r -a read_group_data <<< "$(zcat ${fastq_R1} | head -n 1)"
        else
            IFS=':' read -r -a read_group_data <<< "$(cat ${fastq_R1} | head -n 1)"
        fi
        IFS=" " read_group_data=(${d}{read_group_data[@]})
        java -jar ${picard_path_gitc} AddOrReplaceReadGroups \
            I=${sample_name}.sam \
            O=${sample_name}.bam \
            RGID=${d}{read_group_data[0]:1}.${d}{read_group_data[1]} \
            RGPU=${d}{read_group_data[2]}.${d}{read_group_data[3]} \
            RGLB=${sample_name} \
            RGSM=${sample_name} \
            RGPL=ILLUMINA
        # sort and index
        ${samtools_path_gitc} sort -o ${sample_name}.sorted.bam ${sample_name}.bam
        ${samtools_path_gitc} index ${sample_name}.sorted.bam
    >>>

    output {
        File aligned_bam = "${sample_name}.sorted.bam"
        File aligned_bam_index = "${sample_name}.sorted.bam.bai"
    }

    runtime {
        preemptible: 5
        docker: gitc_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

# mark duplicate reads in bam
task MarkDuplicates {
    File sorted_bam
    String sample_name

    Int disk_size
    Int mem_size_gb
    String picard_docker
    String picard_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${cmd_mem_size_gb}G -jar ${picard_path} MarkDuplicates \
            I=${sorted_bam} \
            O=${sample_name}.marked_duplicates.bam \
            M=${sample_name}.marked_duplicates.metrics
    }

    output {
        File bam = "${sample_name}.marked_duplicates.bam"
    }

    runtime {
        preemptible: 3
        docker: picard_docker
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
    String picard_docker
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
        docker: picard_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

# base quality score recalibration
task BQSR {
    File ref
    File ref_dict
    File ref_index
    File bam
    File bam_index
    String sample_name
    String output_table_name
    String output_bam_name
    Array[File] known_sites
    Array[File] known_sites_indices
    Array[String] intervals_to_exclude

    Int disk_size
    Int mem_size_gb
    String gitc_docker
    String gatk_path
    String picard_path

    Int cmd_mem_size_gb = mem_size_gb - 1
    Array[String] prefixed_intervals_to_exclude = prefix("--excludeIntervals ", intervals_to_exclude)

    command {
        # build BQSR table
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} -T BaseRecalibrator -nt 1 \
            -R ${ref} -I ${bam} \
            -knownSites ${sep=" -knownSites " known_sites} \
            -o ${output_table_name} \
            ${sep=" " prefixed_intervals_to_exclude}

        # install GATK AnalyzeCovariates R dependencies
R --vanilla << CODE
install.packages("gplots", repos="http://cran.us.r-project.org")
install.packages("gsalib", repos="http://cran.us.r-project.org")
install.packages("reshape", repos="http://cran.us.r-project.org")
CODE

        # AnalyzeCovariates
        java -jar ${gatk_path} \
            -T AnalyzeCovariates \
            -R ${ref} \
            --BQSR ${output_table_name} \
            -plots ${sample_name}.bqsr.pdf

        # clean reads, using bqsr if applicable
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T PrintReads \
            -nt 1 \
            -R ${ref} \
            -I ${bam} \
            --BQSR ${output_table_name} \
            -o ${output_bam_name}

        # build index
        java -Xmx${cmd_mem_size_gb}G -jar ${picard_path} \
            BuildBamIndex \
            I=${output_bam_name} \
            O=${output_bam_name}.bai
    }

    output {
        File out = "${output_bam_name}"
        File out_index = "${output_bam_name}.bai"
        File table = "${output_table_name}"
    }

    runtime {
        preemptible: 3
        docker: gitc_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

# merge vcfs before genotyping
task CombineGVCFs {
    File ref
    File ref_dict
    File ref_index
    String interval
    Array[File] vcf_files
    Array[File] vcf_index_files

    String gcvf_out = "combined_gvcfs.vcf.gz"
    String gvcf_out_index = "combined_gvcfs.vcf.gz.tbi"

    Int disk_size
    Int mem_size_gb
    String gatk4_docker
    String gatk4_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        ${gatk4_path} --java-options "-Xmx${cmd_mem_size_gb}G" \
            CombineGVCFs \
            -R ${ref} \
            -L ${interval} \
            -O ${gcvf_out} \
            --variant ${sep=" --variant " vcf_files}
    }
    output {
       File out = gcvf_out
       File out_index = gvcf_out_index
    }

    runtime {
        preemptible: 4
        docker: gatk4_docker
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
    String gcvf_out = "genotyped_gvcfs.vcf.gz"

    Int disk_size
    Int mem_size_gb
    String gatk4_docker
    String gatk4_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        ${gatk4_path} --java-options "-Xmx${cmd_mem_size_gb}G" \
            GenotypeGVCFs \
            -R ${ref} \
            -O ${gcvf_out} \
            --variant ${vcf_file} \
    }
    output {
       String output_vcf_name = "genotyped_gvcfs.vcf.gz"
       String output_vcf_index_name = "genotyped_gvcfs.vcf.gz.tbi"
    }

    runtime {
        preemptible: 4
        docker: gatk4_docker
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
    String gatk_docker
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
        docker: gatk_docker
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
  File ref_fasta
  File ref_fasta_index
  String sample_name

  Int disk_size
  Int mem_size_gb
  Int cmd_mem_size_gb = mem_size_gb - 1

  String gatk_docker
  String gatk_path

  #File index
  String ? intervals
  File ? bqsr_file
  Int ? ploidy
  String ? erc
  String ? extra_hc_params
  String out = "${sample_name}.g.vcf"

  command {
    java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -o ${gvcf_name} \
      ${"--intervals " + intervals} \
      ${"-BQSR " + bqsr_file} \
      -ERC ${default="GVCF" erc} \
      -ploidy ${default="1" ploidy} \
      --interval_padding 100 \
      -variant_index_type LINEAR \
      -variant_index_parameter 128000 \
      ${default="\n" extra_hc_params} \
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
    docker: gatk_docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  parameter_meta {
      gatk: "Executable jar for the GenomeAnalysisTK"
      ref: "fasta file of reference genome"
      sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
      sample_dir: "The sample-specific directory inside output_dir for each sample."
      in_bam: "The bam file to call HaplotypeCaller on."
      intervals: "An array of intervals to restrict processing to."
      bqsr_file: "The full path to the BQSR file."
      erc: "Mode for emitting reference confidence scores."
      extra_hc_params: "A parameter that allows users to pass any additional paramters to the task."
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
    String snpeff_docker = "maxulysse/snpeff:1.3"
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
        docker: snpeff_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
