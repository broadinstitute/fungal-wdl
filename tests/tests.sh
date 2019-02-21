#!/usr/bin/env sh

# export directory where cromwell and wdltool locate to JARS

WDL=workflows/gatk3_germline_snps_indels.wdl

java -jar ${JARS}/wdltool-0.14.jar validate $WDL
java -jar ${JARS}/cromwell-30.1.jar run -i tests/test.json $WDL 
