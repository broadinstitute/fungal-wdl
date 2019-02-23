#!/usr/bin/env sh

# export directory where cromwell and wdltool locate to JARS
set -e

WDL=workflows/gatk3_germline_snps_indels.wdl

printf "##########\nBegin test\n##########\n"

java -jar ${JARS}/wdltool-0.14.jar validate $WDL
printf "##########\nSyntax test done. \n##########\n"

java -jar ${JARS}/cromwell-30.1.jar run -i tests/mini_test.json $WDL
printf "Mini test done. \n##########"

java -jar ${JARS}/cromwell-30.1.jar run -i tests/test.json $WDL
printf "\n##########\nFull test done.\n##########\n"

# scp -r cromwell-executions/GATK3*/*/call-HaplotypeCaller/shard-0/execution/*g.vcf xiaoli@login:/gsap/cdcfungal/fc_pipeline
