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
