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
