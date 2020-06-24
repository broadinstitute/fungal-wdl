# Fungal Firecloud (Terra) tutorial. SNP calling GATK3 workflow.

### Getting Started in Terra documentation could be found [here](https://support.terra.bio/hc/en-us/categories/360001728852-Getting-Started-in-Terra).

## Quick start
1. Register a Terra account [here](https://app.terra.bio/#workspaces). Simply sign in via your gmail account.

<a href="https://imgur.com/CADZ6dT"><img src="https://imgur.com/CADZ6dT.png" title="source: imgur.com" /></a>

2. In My Workspace tab, search for `broad-fungal-gatk3` workspace. Click the workspace name and Clone it.

<a href="https://imgur.com/RfzuC02"><img src="https://imgur.com/RfzuC02.png" title="source: imgur.com" /></a>

3. Rename the workspace for your own analysis project (e.g. `my_gatk3_workspace`). Use your specific billing project (e.g. `broad-fungal-firecloud`) and authorization domain. 

<a href="https://imgur.com/CdSAnuR"><img src="https://imgur.com/CdSAnuR.png" title="source: imgur.com" /></a>

There should be a free $300 credit along with your newly registered account. Use the billing project with that to test it.

The workspace contains reference files for Candida auris strain B8441 (Clade I), including the fasta file and its indice and gff files (as workspace attributes), two example BAMs, and a method configuration `fungal-variant-call-gatk3` for illustration purposes. You can customize them for your own analysis. Read more [below](#CUS).

4. Set up reference and sample files: 

go to `data` configuration, click the `Workspace data` to set up REFERENCE files. Define keys using exactly the same name as in the example and link to the corresponding file in google cloud (see below how to import REFERENCE and sample files into google cloud).

<a href="https://imgur.com/p29VGBx"><img src="https://imgur.com/p29VGBx.png" title="source: imgur.com" /></a>

To add samples (read samples as participant set) make sure that `data` tables are empty. Then, in `data` click `+` to add a table. Then drag or select the `*.tsv` file with the sample name (i.e. participant names) and google cloud path to the google cloud bucket (`gs://broad_fungal_data_sets/my_dataset`).

<a href="https://imgur.com/0iR0op5"><img src="https://imgur.com/0iR0op5.png" title="source: imgur.com" /></a>

The `*.tsv` file should look like this. Make sure to include `entity:participant_id	bam_file`:

```
entity:participant_id	bam_file
B13273	gs://broad_fungal_data_sets/cand_hae/B13273.sorted.bam
B13444	gs://broad_fungal_data_sets/cand_hae/3001030171.sorted.bam
B13704	gs://broad_fungal_data_sets/cand_hae/3001033241.sorted.bam
B13909	gs://broad_fungal_data_sets/cand_hae/3000029013.sorted.bam
B15318	gs://broad_fungal_data_sets/cand_hae/3001501541.sorted.bam
```

For GATK4. The `*.tsv` file could be something like:
```
membership:participant_set_id	participant
funga-variant-call-gatk4_2020-01-28	ACTA3157-D
funga-variant-call-gatk4_2020-01-28	ACTA3092-D
``` 


5. Set up the analysis:

go to `Tools` and click the method configuration `fungal-variant-call-gatk3`.

<a href="https://imgur.com/40RH2EM"><img src="https://imgur.com/40RH2EM.png" title="source: imgur.com" /></a>

select the dataset (samples | participant_set).

<a href="https://imgur.com/s6C3A7y"><img src="https://imgur.com/s6C3A7y.png" title="source: imgur.com" /></a>

`variables` in `Tools` should not be changed if `Keys` and `participant_id	bam_file` were specified as described above. The only `variable` that could be changed is `run_name` where you can define the name of your analysis in the `Attribute` field.

click `RUN_ANALYSIS`

<a href="https://imgur.com/gOnNzdO"><img src="https://imgur.com/gOnNzdO.png" title="source: imgur.com" /></a>



6. Monitoring job progress:

Click the `View` text in `JOB HISTORY`, and you can see the details of your job.

<a href="https://imgur.com/PYwgK1p"><img src="https://imgur.com/PYwgK1p.png" title="source: imgur.com" /></a>

<a href="https://imgur.com/h2hGlW3"><img src="https://imgur.com/h2hGlW3.png" title="source: imgur.com" /></a>

Click individual steps in `LIST VIEW` and `TIMING DIAGRAM` to see more details of the analysis run and progress.

The last step of the analysis is `HardFiltration`.

The output is `[your_run_name].hard_filtered.vcf.gz`. The google cloud path for this file can be find in `OUTPUTs`.

Copy that file into the server using gsutil.

Go to the command line and copy this combined VCF:
```
use .google-cloud-sdk-20160101
gsutil cp gs://fc-secure-f76aaa59-a5c8-42e5-8ec4-db44ad254b33/8cf8a219-f855-4e93-b0c9-89795dc6d37d/GATK3_Germline_Variants/5bec624a-afb8-4f69-9682-58b1cd790cea/call-HardFiltration/Cand_hae.hard_filtered.vcf.gz /my_path/
```

Enjoy!!!




To use your own reference and dataset, you will need to create a google bucket and upload files there, and update the JSON files accordingly. Read more [below](#REF).

## <a name="CUS">Customize your own workspace</a>
### <a name="REF">Reference files</a>

Create index file for your reference (${REF}.fasta)
```
# Step 1:
bwa index ${REF}
samtools faidx ${REF}
java -jar /cil/shed/apps/external/picard/current/bin/CreateSequenceDictionary.jar R= ${REF}.fasta O= ${REF}.dict
```

If you are using GATK4. The input bams should be unaligned bam. Use Picard tools to create a unaligned bams:
```
#picard RevertSam I=x.bam O=x.sorted.bam SO=queryname
java -jar /cil/shed/apps/external/picard/current/bin/picard.jar RevertSam I=B1-{SAMPLE}.sorted.bam O={SAMPLE}.unaligned.bam SO=queryname
``` 

### <a name="DATA">Import your data</a>

Go to Google Cloud Platform
https://console.cloud.google.com/storage/browser?project=gcid-cromwell

Upload BAM files:
Create a new folder for your dataset (`my_dataset`) in `broad_fungal_data_sets`
<a href="https://imgur.com/XsDXltb"><img src="https://imgur.com/XsDXltb.png" title="source: imgur.com" /></a>

<a href="https://imgur.com/8vFBd1l"><img src="https://imgur.com/8vFBd1l.png" title="source: imgur.com" /></a>

Go to the command line and upload data (bam files) to this new folder:
```
use .google-cloud-sdk-20160101
gsutil -m cp *bam gs://broad_fungal_data_sets/my_dataset # -m is multithreading
```

Upload REFERENCE files:
Create a new folder for your reference (`my_reference`) in `broad_fungal_references`
<a href="https://imgur.com/UXJcuVl"><img src="https://imgur.com/UXJcuVl.png" title="source: imgur.com" /></a>

Go to the command line and upload all files created in Step1 to this new folder:
```
use .google-cloud-sdk-20160101
gsutil -m cp ${REF}.* gs://broad_fungal_references/my_reference
```
