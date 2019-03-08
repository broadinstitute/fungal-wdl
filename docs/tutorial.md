# Fungal Firecloud tutorial

## Quick start
1. Register a Firecloud account [here](https://portal.firecloud.org). Simply sign in via your gmail account.

<a href="https://imgur.com/s2nbP1D"><img src="https://i.imgur.com/s2nbP1D.png" title="source: imgur.com" /></a>

2. In My Workspace tab, search for `broad-fungal-gatk3` workspace. Click the workspace name to get in the workspace page.

<a href="https://imgur.com/TOwmxV2"><img src="https://i.imgur.com/TOwmxV2.png" title="source: imgur.com" /></a>

3. Clone the workspace and rename it for your own analysis project.

<a href="https://imgur.com/POp7SvG"><img src="https://i.imgur.com/POp7SvG.png" title="source: imgur.com" /></a>

There should be a free $300 credit along with your newly registered account. Use the billing project with that to test it.

<a href="https://imgur.com/xFmJbZI"><img src="https://i.imgur.com/xFmJbZI.png" title="source: imgur.com" /></a>

The workspace contains reference files, including the fasta file and its indice and gff files (as workspace attributes), two example BAMs, and a method configuration `fungal-variant-call-gatk3` for illustration purposes. You can customize them for your own analysis. Read more [below](#CUS).

4. Launch an analysis: go to method configuration, click the `fungal-variant-call-gatk3`.

<a href="https://imgur.com/MpPFahu"><img src="https://i.imgur.com/MpPFahu.png" title="source: imgur.com" /></a>

Select participant set used for the analysis `TestSet` here as an example. Click launch to submit the job.

<a href="https://imgur.com/ivKFcj9"><img src="https://i.imgur.com/ivKFcj9.png" title="source: imgur.com" /></a>

5. Monitoring job progress:

Click the `View` text, and you can see the details of your job.

<a href="https://imgur.com/QoxCEPo"><img src="https://i.imgur.com/QoxCEPo.png" title="source: imgur.com" /></a>

Click `show` tag to see more details, such as time graph and individual job status.

<a href="https://imgur.com/YbzdUD9"><img src="https://i.imgur.com/YbzdUD9.png" title="source: imgur.com" /></a>
To use your own reference and dataset, you will need to create a google bucket and upload files there, and update the JSON files accordingly. Read more [below](#REF).

## <a name="CUS">Customize your own workspace</a>
### <a name="REF">Reference files</a>


### <a name="DATA">Import your data</a>

Using data model.

### <a name="CONF">Job configuration</a>


Via Firecloud APIs.
