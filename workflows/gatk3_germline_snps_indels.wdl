#
#  gatk3_germline_snps_indels.wdl:
#     - standardized pipeline for calling germline SNPs and INDELs.
#     - imports haplotypecaller subworkflow to scatter SNP/INDEL calls over intervals,
#       during scatter over samples.
#     - developed by Malaria Group, IDMP, Broad Institute.
#
import 'https://api.firecloud.org/ga4gh/v1/tools/broad-fungal-methods:haplotypecaller_gatk3-scatter_intervals/versions/1/plain-WDL/descriptor' as HaplotypeCaller
# import '../gatk3_haplotype_caller/haplotypecaller_gatk3-scatter_intervals.wdl' as HaplotypeCaller

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

    # for base quality score recalibration
    Boolean do_bqsr
    Array[File] known_sites
    Array[File] known_sites_indices
    Array[String] bqsr_intervals_to_exclude

    # variant quality control param
    # either "vqsr" or "hard_filtering"
    String variant_qc

    # vqsr params
    # if variant_qc == "vqsr"
    # all of these params are required
    Float ts_filter_snp
    Float ts_filter_indel
    Int snp_max_gaussians
    Int indel_max_gaussians
    Int vqsr_mapping_qual_cap
    Array[String] snp_resources
    Array[String] indel_resources
    Array[String] snp_annotations
    Array[String] indel_annotations

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

        if (do_bqsr) {
            # base quality score recalibration
            call BQSR {
                input:
                ref = ref,
                ref_dict = ref_dict,
                ref_index = ref_index,
                bam = ReorderBam.out,
                bam_index = ReorderBam.out_index,
                sample_name = sample_name,
                known_sites = known_sites,
                known_sites_indices = known_sites_indices,
                intervals_to_exclude = bqsr_intervals_to_exclude,
                output_table_name = sample_name + ".bqsr.table",
                output_bam_name = sample_name + ".bsqr.bam",
                gitc_docker = gitc_docker,
                gatk_path = gitc_path_to_gatk,
                picard_path = gitc_path_to_picard,
                mem_size_gb = large_mem_size_gb,
                disk_size = disk_size
            }
        }
        # to do: targetcreator
        # to do: indel realignment

        # call subworkflow
        call HaplotypeCaller.HaplotypeCallerGvcf_GATK3 {
            input:
            ref_fasta = ref,
            ref_dict = ref_dict,
            ref_fasta_index = ref_index,
            input_bam = select_first([BQSR.out, ReorderBam.out]),
            input_bam_index = select_first([BQSR.out_index, ReorderBam.out_index]),
            scattered_calling_intervals_list = interval_files_list,
            merge_gvcfs_mem_size_gb = large_mem_size_gb,
            haplotypecaller_mem_size_gb = med_mem_size_gb,
            merge_gvcfs_disk_size = med_disk_size,
            haplotypecaller_disk_size = med_disk_size,
            gatk_docker = gatk_docker,
            gatk_path = gatk_path_to_gatk,
            picard_docker = gitc_docker,
            picard_path = gitc_path_to_picard
        }
    } # end scatter

    scatter (interval in read_lines(interval_list)) {
        call CombineGVCFs {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index,
            interval = interval,
            vcf_files = HaplotypeCallerGvcf_GATK3.output_gvcf,
            vcf_index_files = HaplotypeCallerGvcf_GATK3.output_gvcf_index,
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
            interval = interval,
            vcf_file = CombineGVCFs.out,
            vcf_index_file = CombineGVCFs.out_index,
            gatk4_docker = gatk4_docker,
            gatk4_path = gatk4_path_to_gatk,
            mem_size_gb = med_mem_size_gb,
            disk_size = disk_size
        }
    }

    call GatherVCFs {
        input:
        vcf_files = GenotypeGVCFs.out,
        gitc_docker = gitc_docker,
        picard_path_gitc = gitc_path_to_picard,
        tabix_path_gitc = "/usr/gitc/tabix",
        mem_size_gb = extra_large_mem_size_gb,
        disk_size = large_disk_size
    }

    # variant quality control
    if (variant_qc == "vqsr") {
        # variant quality score recalibration
        # snp vqsr
        call VQSR as SnpVQSR {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index,
            gvcf = GatherVCFs.out,
            gvcf_index = GatherVCFs.out_index,
            output_filename = "${run_name}.snp_vqsr.g.vcf.gz",
            output_index_filename = "${run_name}.snp_vqsr.g.vcf.gz.tbi",

            mode = "SNP",
            resources = snp_resources,
            resource_files = known_sites,
            resource_file_indices = known_sites_indices,
            annotations = snp_annotations,
            ts_filter = ts_filter_snp,
            max_gaussians = snp_max_gaussians,
            mapping_qual_cap = vqsr_mapping_qual_cap,

            gatk4_docker = gatk4_docker,
            gatk4_path = gatk4_path_to_gatk,
            mem_size_gb = med_mem_size_gb,
            disk_size = disk_size
        }
        # indel vqsr
        call VQSR as IndelVQSR {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index,
            gvcf = SnpVQSR.out,
            gvcf_index = SnpVQSR.out_index,
            output_filename = "${run_name}.indel_vqsr.g.vcf.gz",
            output_index_filename = "${run_name}.indel_vqsr.g.vcf.gz.tbi",

            mode = "INDEL",
            resources = indel_resources,
            resource_files = known_sites,
            resource_file_indices = known_sites_indices,
            annotations = indel_annotations,
            ts_filter = ts_filter_indel,
            max_gaussians = indel_max_gaussians,
            mapping_qual_cap = vqsr_mapping_qual_cap,

            gatk4_docker = gatk4_docker,
            gatk4_path = gatk4_path_to_gatk,
            mem_size_gb = med_mem_size_gb,
            disk_size = disk_size
        }
    }
    if (variant_qc == "hard_filtering") {
        call HardFiltration {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index,
            vcf = GatherVCFs.out,
            vcf_index = GatherVCFs.out_index,
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
    # add variant annotations using SnpEff
    if (do_snpeff) {
        call SnpEff {
            input:
            ref = ref,
            organism_gff = organism_gff,
            output_vcf_name = "${run_name}.snpeff.g.vcf",
            vcf = select_first([IndelVQSR.out, HardFiltration.out]),
            mem_size_gb = extra_large_mem_size_gb,
            disk_size = extra_large_disk_size
        }
    }

    output {
        File gvcf = select_first([SnpEff.out, IndelVQSR.out, HardFiltration.out])
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
    String interval
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
            -L ${interval} \
            -O ${gcvf_out} \
            --variant ${vcf_file} \
            --only-output-calls-starting-in-intervals \
            --use-new-qual-calculator
    }
    output {
       File out = gcvf_out
    }

    runtime {
    	preemptible: 4
        docker: gatk4_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task GatherVCFs {
    Array[File] vcf_files

    Int disk_size
    Int mem_size_gb
    String gitc_docker
    String picard_path_gitc
    String tabix_path_gitc

    String output_vcf_name = "gathered.vcf.gz"
    String output_vcf_index_name = "gathered.vcf.gz.tbi"
    Int cmd_mem_size_gb = mem_size_gb - 1
    #File vcf_file_list = write_lines(vcf_files)

    command {
        java "-Xmx${cmd_mem_size_gb}G" -jar ${picard_path_gitc} \
            GatherVcfs \
            I=${sep=" I=" vcf_files} \
            O=${output_vcf_name}

        ${tabix_path_gitc} -p vcf ${output_vcf_name}
    }
    output {
        File out = output_vcf_name
        File out_index = output_vcf_index_name
    }

    runtime {
        preemptible: 4
        docker: gitc_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

# variant quality score recalibration
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
task VQSR {
    File ref
    File ref_dict
    File ref_index
    File gvcf
    File gvcf_index
    String output_filename
    String output_index_filename

    String mode
    Float ts_filter
    Array[String] resources
    Array[String] annotations
    Array[File] resource_files
    Array[File] resource_file_indices

    Int max_gaussians
    Int mapping_qual_cap

    Int disk_size
    Int mem_size_gb
    String gatk4_docker
    String gatk4_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    String vqsr_file = "${mode}.recal"
    String rscript_file = "${mode}.plots.R"
    String tranches_file = "${mode}.tranches"

    command {
        # build vqsr file
        ${gatk4_path} --java-options "-Xmx${cmd_mem_size_gb}G" \
            VariantRecalibrator \
            -R ${ref} \
            -V ${gvcf} \
            -mode ${mode} \
            -O ${vqsr_file} \
            --tranches-file ${tranches_file} \
            --rscript-file ${rscript_file} \
            --resource ${sep=" --resource " resources} \
            -an ${sep=" -an " annotations} \
            --max-gaussians ${max_gaussians} \
            -mq-cap ${mapping_qual_cap}

        # apply vqsr
        ${gatk4_path} --java-options "-Xmx${cmd_mem_size_gb}G" \
            ApplyVQSR \
            -R ${ref} \
            -V ${gvcf} \
            --truth-sensitivity-filter-level ${ts_filter} \
            --tranches-file ${tranches_file} \
            --recal-file ${vqsr_file} \
            -mode ${mode} \
            -O ${output_filename}
    }

    output {
        File vqsr = vqsr_file
        File rscript = rscript_file
        File tranches = tranches_file
        File out = output_filename
        File out_index = output_index_filename
    }

    runtime {
    	preemptible: 3
        docker: gatk4_docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

# hard-filter a vcf, if vqsr not available
# http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
task HardFiltration {
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

# annotate variants
# Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
task SnpEff {
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
