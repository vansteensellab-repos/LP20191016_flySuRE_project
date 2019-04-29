# 190425

## get data from Heidelberg

Get data and store in `~/projects/LP190425_flySuRE/data/LP190425_iPCR_67`:

```
while read -r f ; do echo $f; wget $f; done < fastq_filelist_ftp_LP190425.txt
```

## Fastqc

```
procHTS_LP170705.sh -c 20
```

# 190426

I wrote an email to Matteo, Mattia to establish clear overview of project details.

Also talked to Marcel about the project. He pointed out an overview of the iPCR
samples for the 10 fly strains in his labguru journal:
<https://nki.labguru.com/knowledge/projects/14/milestones/24>.

Also talked to Joris. He explained that the "drop in complexity" mentioned in
an email by Matteo probably referred to the sequencing machine they used being
incapable of dealing with the low sequence complexity in the read structure due
to the plasmid backbone in the reads.

# 190429

## Fastqc / multiqc

Running `fastqc` took some investment in setting up a new conda env, `NGS_QC`,
as fastqc was generating errors (missing font) in the base env.

After this fastqc ran just fine:
```
procHTS_LP170705.sh -c 20
```
Next, `multiqc` also ran without probels:
```
multiqc ../../data/LP190425_iPCR_67/
```
The resulting report (analyses/LP190425_check_iPCR_fastq/multiqc_report.html)
suggests there might be some problems with a few samples:
* Some samples have relatively high adapter content (Illumina universal adapter)
* For all samples adapter content appears to increase significantly after position 102 or so!
* Some sample have somewhat higher duplication levels.

The indicated samples are all the GRP samples!

One possibility might be that the GRP samples have somewhat shorter fragments.
This explains that they tend to read into adapter more frequently.  
To confirm I need to parse the reads to check the gDNA length distribution of the samples.


