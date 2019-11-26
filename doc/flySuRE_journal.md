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

# 190606

## Data looks very bad


### FASTQ files
I looked at fastq files and eye-balled the expected fixed read structure. I can
see that the expected sequences are hardly present, and if present mostly not
at expected position. I revisited the FASTQC reports and the overal read
quality is lower than normal (max 35, dropping to 18 at the regions where fixed
backbone is expected. Normally we see ~38, slowly dropping to 24 over the
readlength).

It seems to me that the reported problems wrt the sequencing machine and the
lack of complexity in the sequences (due to fixed backbone sequence) was very
detrimental.

### Alignments
I aligned some sequences using BLAST and they perfectly align with Dme, so no
problem there. 

Also aligned 1e6 reads using bowtie against dm5 and the results look good:
```
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/raw/LP190425_iPCR_67/raw_data
bowtie2 -s 10000 -u 1000000 -5 10 -x /DATA/data/bowtie2/dmel_r5_genome_2L_YHet -U 67_SuREseq_DM09_GDL_I02_2.fastq.gz -p 30 | samtools view - -b -o /tmp/tt.bam                                                                                                                                                                                                  
1000000 reads; of these:
  1000000 (100.00%) were unpaired; of these:
    254188 (25.42%) aligned 0 times
    645797 (64.58%) aligned exactly 1 time
    100015 (10.00%) aligned >1 times
74.58% overall alignment rate
```

### Coverage

I can estimate the coverage of the data over the genome. I use the 1e6 alignments of the reverse reads only:

```
samtools depth -a /tmp/tt.srt.bam | awk ' {depth[$3]++} END{for(d in depth){print d, depth[d];tot+=d*depth[d]} print tot/NR}'

.
.
.
511 2
513 1
514 4
517 3
519 3
520 2
521 1
0.801608
```
Roughly 50% of the genome is not at all covered. The average coverage is 0.8, and the overall distribution looks good to me. 

Actual coverage:
* actual fragment length is ~ twice the readlength in calculation above, so cov=1.6
* total number of reads is 25e6, with a unique alignment rate of 65%: 16e6 reads in total, so cov=25

# 190726

## New data after Phi-X spike in

I received a mail from Matteo (dd 20190723 and subsequent days) announcing that the sequencing is finished. There were some issues downloading the data as the remote ftp server was updated or something.\\
I downloaded the data into a new data dir, checked file checksums, checked fastqc html reports, etc:
```
cd /home/ludo/projects/LP190425_flySuRE/data/raw/
mkdir LP190724_iPCR_complexity_test
cd LP190724_iPCR_complexity_test
wget  -r ftp://ftp-exchange.embl.de/pub/exchange/perino/iPCR/
mv ./ftp-exchange.embl.de/pub/exchange/perino/iPCR/ .
rm -rf ftp-exchange.embl.de
```

## Estimate barcode complexity of new data

The reads are 66bp long. The only things I can do is:

* check the fixed backbone sequence
* check the diversity (complexity) of the barcode sequence

### Fixed backbone sequence

First I will trim the reads on the backbone sequence 'CCTAGCTAACTATAACGGTCCTAAGGTAGCGAACCAGTGAT'

```
cd /home/ludo/projects/LP190425_flySuRE/data/intermediate
mkdir LP190726_iPCR_190724_trimmed_fastq
cd LP190726_iPCR_190724_trimmed_fastq
for f in ../../raw/LP190724_iPCR_complexity_test/*_sequence.txt.gz; do
  fnew=$( basename $f )
  echo $fnew
  cutadapt --discard-trimmed -j 1 -O5 -g CCTAGCTAACTATAACGGTCCTAAGGTAGCGAACCAGTGAT -o ${fnew%.txt.gz}.fq $f --info-file ${fnew%.txt.gz}.info &> ${fnew%.txt.gz}.stats
done
```

* number of untrimmed per sample
  Approximately 2% of the reads do not contain a (good enough) backbone sequence
* ~95% of all reads contain a recognizabble barcode of exactly 20bp
`for f in *info; do echo $f; cat $f | awk -F"\t" ' $3==20{c++} END{print 100*(1-(c/NR))} '; done | awk 'NR%2==1{nm=$1; next}{print nm"\t"$1}' > LP20190726_percentage_untrimmed.txt`

### Barcode complexity

```
echo -e "file\t#reads\t#BC_20\tperc_uniq_BC\tBC_count_distrib" > LP20190726_barcode_complexity.txt
for f in *info; do
  awk -F"\t" ' 
  BEGIN{OFS="\t"} 
  $2==-1{next} 
  { if($5~/[ACGT]{20}/){c++; a[$5]++} } 
  END{
    for(k in a){b[a[k]]++}
    st="1:"b[1];
    for(k=2;k<=length(b);k++){st=st","k":"b[k]} 
    print FILENAME,NR,c,100*(length(a)/c),st
  } ' $f
done >> LP20190726_barcode_complexity.txt
```

Approximately 2-3% of barcodes are duplicated. This seems very good.

### Conclusion

The data looks very good, with good complexity of barcodes.

# 20190805

## Report to Matteo

Matteo sent me an email and asked about the results. I sent him an email saying
all is perfect and he should continue with the sequencing.

# 20190906

Matteo sent me an email saying new data is available, both iPCR data (mail dd
Sep 3, 2019, 11:21 AM) and cDNA data (mail dd Sep 5, 2019, 2:47 PM)!

## iPCR data
```
mkdir /DATA/usr/ludo/projects/LP190425_flySuRE/data/raw/LP20190906_iPCR_resequencing_I
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/raw/LP20190906_iPCR_resequencing_I
dos2unix iPCR_ftp
cat  iPCR_ftp.txt | grep iPCR_mix | while read -r f ; do echo $f; wget $f; done
```

## cDNA data
```
mkdir /DATA/usr/ludo/projects/LP190425_flySuRE/data/raw/LP20190906_cDNA_resequencing_I
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/raw/LP20190906_cDNA_resequencing_I
dos2unix SuRE_ftp.txt
cat  SuRE_ftp.txt | while read -r f ; do echo $f; wget $f; done
```

# 20190909

At first the data was not world-readable, which Matteo fixed. I downloaded the data over the weekend.

## MultiQC

I ran multiQC on the fastQC files which Matteo already generated, and the data
looks good. The only thing that strikes me is that  there appears to be quite
high level of duplicate data, but I'm not sure if it is higher than normal???

```
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/raw/LP20190906_iPCR_resequencing_I
# sent output to the analysis directory so that it will be put under git control
multiqc -o /DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20190909_resequencing-I_1st_samples/ . 

cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/raw/LP20190906_cDNA_resequencing_I
multiqc -o /DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20190909_resequencing-I_1st_samples/cDNA_resequencing_I_multiQC/ .
```

## Retrieve VCF data from EMBL

Mattia sent me a link to a dropbox file. I downloaded the vcf and the tbi file
to spinel and moved it to
`lovelace:/DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190909_VCF_fly_Mattia`.

## Genome refrence sequence

Mattia used a specific ref seq for his variant analysis; `dm6.AE017196.fasta`.
I sent him an email requesting a copy.

# 20190912

## Retrieve ref seq

I received an email from Mattia (dd Wed, Sep 11, 11:06 AM). He explained that
he used a combination of dm6 and a Wolbachia genome sequence for the variant
annotation. I don't know why he did so.

He included a link for the download of the combined fasta:

```
mkdir /DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190912_DM6_Wolbachia_refSeq_Mattia
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190912_DM6_Wolbachia_refSeq_Mattia
wget ftp://ftp-exchange.embl.de/pub/exchange/forneris/Ludo/dm6.AE017196.fasta
```

Aha; i see from the fasta file that Wolbachia is a:
*endosymbiont of Drosophila melanogaster*\\
Now it makes more sense.

# 20190913

## Bowtie2 indices
```
mkdir /DATA/usr/ludo/projects/LP190425_flySuRE/data/intermediate/LP20190913_Dmel_Wolb_bowtie2Indices
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/intermediate/LP20190913_Dmel_Wolb_bowtie2Indices
bowtie2-build ../../external/LP20190912_DM6_Wolbachia_refSeq_Mattia/dm6.AE017196.fasta DM6_Wolbachia
```

# 20190916

## Generate parental reference sequences for each of the cell lines

I need to follow the methodology I used for the INDEL pipeline on the SuRE42-54
data in
`projects/LP140430_SureSeq_JvArensbergen/analyses/LP171023_createAltRefs`.

This work is done in `projects/LP190425_flySuRE/analyses/LP20190916_generateParentalRefSeq`.

# 20190925

## Generate parental seq's, con't

I fuond that the VCF files are not phased. But after discussing with Joris and
with Mattia it is clear that the flies are (near) homozygous. So, the idea will
be to generate a single ref seq per fly strain. Before doing this I want to
filter the VCF files such that all non-homozygous loci are discarded.

## Filter VCF files

There is a single VCF file, containing 13 columns with GT specs. The format of
the columns is specified in the column 'GT', and is: `GT:AD:DP:GQ:PL`. The data
is not phased so the separator is `/`. I need to split the sample columns using
the pattern (from https://github.com/lpagie/SuRE-K562-pipeline/commit/65a6bd1799d1e17d0a9a705d9af2cedace72351a):
```
  gt=gensub(/([^:]+):.*/,"\\1","g",$GT)
  split(gt, alleles, /[|\\]/)
```

To discard the heterozygous SNPs I will loop through all, check the sample
column for all 13 samples, check for heterozygosity, and discard the
heterozygous ones.

```
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190909_VCF_fly_Mattia

zcat Haplotype_joint_call_SuRe_lines_lenient_filtering.vcf.gz | \
gawk -F"\t" '
# copy header lines to output
/^#/ { 
  # print; 
  next 
}
{
  homoz=1 # homozygosity flag for current SNP
  for (i=10; i<=NF; i++) { # loop over sample columns
    gt=gensub(/([^:]+):.*/,"\\1","g",$i) # discard all fields except first
    split(gt, alleles, /[|\/]/) # split GT into alleles
    if (alleles[1] != alleles[2]) {
      print
      next
    }
  }
  # print
  # exit
} '

```

### Some stats on VCF and filtering

The number of heterozygous SNPs is ~1e6:

```
zcat Haplotype_joint_call_SuRe_lines_lenient_filtering.vcf.gz | gawk -F"\t" '
# copy header lines to output
/^#/ {  
  # print;
  next 
}
{
  homoz=1 # homozygosity flag for current SNP
  for (i=10; i<=NF; i++) { # loop over sample columns
    gt=gensub(/([^:]+):.*/,"\\1","g",$i) # discard all fields except first
    split(gt, alleles, /[|\/]/) # split GT into alleles
    if (alleles[1] != alleles[2]) {
      print
      next
    }
  } 
  # print
  # exit
} ' | wc -l
# 905559
```

The total number of SNPs is less than 3e6:
```
zcat Haplotype_joint_call_SuRe_lines_lenient_filtering.vcf.gz | wc -l
# 2707597
```

So, it is a bit much to discard 33% of all SNPs

# 20190926

## Con't stats of VCF and discarding SNPs

The number and percentage of SNPs that will be discarded per sample, with the percentage of SNPs left after combining all (assuming independent overlapping heterozygous SNPs):

```
zcat Haplotype_joint_call_SuRe_lines_lenient_filtering.vcf.gz | gawk -F"\t" '
# copy header lines to output
/^#/ {  
  # print;
  next 
}         
{         
  homoz=1 # homozygosity flag for current SNP
  for (i=10; i<=NF; i++) { # loop over sample columns
    gt=gensub(/([^:]+):.*/,"\\1","g",$i) # discard all fields except first
    split(gt, alleles, /[|\/]/) # split GT into alleles
    if (alleles[1] != alleles[2]) {
      print i
      #next
    }
  } 
  # print
  # exit
} ' | sort | uniq -c | awk ' BEGIN{left=100} { p=$1/2707597; left=left*(1-p); print $1"\t"$2"\t"p*100 } END{print left}'
34724   10      1.28247
18078   11      0.667677
17719   12      0.654418
220849  13      8.15664
202699  14      7.48631
15969   15      0.589785
55019   16      2.03202
51531   17      1.9032
68844   18      2.54262
158633  19      5.85881
35979   20      1.32882
305927  21      11.2988
15263   22      0.56371
63.1417
```

The predicted 63.14% left after discarding all SNPs heterozygous in any of the samples is very close to the observed `(1-905559/2707597)*100 = 66.55%`

### CONCLUSION

I need to discard SNPs per sample. When analyzing significance of SNPs the
number of underlying samples will differ per SNP. I guess this is not a real
problem, as the loss of power will affect the statistical significance exactly
the same as if there are fewer gDNA fragments overlapping the SNP in the
library.

## Allele frequency stats

Per SNP; how many samples and how many reference and alternative alleles are available:

```
zcat Haplotype_joint_call_SuRe_lines_lenient_filtering.vcf.gz | gawk -F"\t" '
BEGIN {
  ORS="\t"
}
# copy header lines to output
/^#/ {
  # print;
  next
}
{
  for(i=0;i<100;i++){cnt[i]=0}
  for (i=10; i<=NF; i++) { # loop over sample columns
    gt=gensub(/([^:]+):.*/,"\\1","g",$i) # discard all fields except first
    split(gt, alleles, /[|\/]/) # split GT into alleles
    if (alleles[1]!="." && alleles[1] == alleles[2]) {
      cnt[alleles[1]]++
    }
  }
  split($5,alt,",")
  nalleles=length(alt)+1
  ORS="\t"
  print nalleles
  for(k=0;k < nalleles;k++) {print cnt[k]}
  ORS=""
  print "\n"
} ' | sort | uniq -c | sort -gr | gawk ' {c+=$1; print c"\t"(c/2705697)*100"\t"$0} '
742436  27.4397  742436 2       12      1       
989174  36.5589  246738 2       11      2       
1157339 42.7742  168165 2       12      0       
1306742 48.2959  149403 2       10      3       
1413194 52.2303  106452 2       11      1       
1518180 56.1105  104986 2       9       4       
1598757 59.0885   80577 2       8       5       
1675120 61.9108   76363 2       10      2       
1742302 64.3938   67182 2       7       6       
1806495 66.7663   64193 2       9       3       
1861637 68.8043   55142 2       6       7       
1915227 70.785    53590 2       8       4       
1963657 72.5749   48430 2       5       8       
2009312 74.2623   45655 2       7       5       
2053020 75.8777   43708 2       4       9       
2093087 77.3585   40067 2       3       10      
2131993 78.7964   38906 2       6       6       
2170436 80.2173   38443 2       2       11      
2208714 81.632    38278 2       1       12      
2241676 82.8502   32962 2       5       7       
2268720 83.8497   27044 2       4       8       
2290423 84.6519   21703 2       3       9       
2312078 85.4522   21655 2       11      0       
2331918 86.1855   19840 2       10      1       
.
.
.
.
```

By far most SNPs have ref/alt allele distribution stronly skewed towards ref being overrepresented in fly strains.

# 20191004

## Generate per sample VCF files

Given the distribution of the heterozygous SNPs I need to discard the
heterozygous SNPs per sample. Not, as I planned inititally, for the entire set
of fly-strains. So, first I'll split the VCF per strain and then remove
heterozygous SNPs.

### Split VCF file per sample

**Note:** Later (dd 20191008) I found that the VCFs contain a blank line before
the #CHROM line. I fixed that below but for future reference I should print the
header (all ^## lines) without an additional newline, probably by changing ORS
to "")

```
awk -F"\t" '
BEGIN {
  OFS="\t"
}
/^##/ {
  header=header""$0"\n"
  next
}
/^#CHROM/ {
  for (i=1; i<=NF; i++) {
    if ( index("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT", $i)==0 ) {
      sampleCol[$i]=i
    }
  }
  for (k in sampleCol) {
    print header > "Haplotype_joint_call_SuRe_lines_lenient_filtering_" k ".vcf"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" k >> "Haplotype_joint_call_SuRe_lines_lenient_filtering_" k ".vcf"
  }
  next
}
{
  for (k in sampleCol) {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$sampleCol[k] >> "Haplotype_joint_call_SuRe_lines_lenient_filtering_" k ".vcf"
  }
}
' <( zcat Haplotype_joint_call_SuRe_lines_lenient_filtering.vcf.gz )

for f in *vcf; do
  cat $f | awk -F"\t" '
    BEGIN { OFS="\t" }
    /#/ { print; next }
    {
      sample=gensub(/([^:]+):.*/,"\\1","g",$10)
      split(sample, allele, /[|\/]/)
      if (allele[1] == allele[2]) {print}
    }
  ' > ${f%.vcf}_homoz.vcf
done

parallel gzip ::: *vcf

find -name \*gz -not -name "*homo*" -not -name Haplotype_joint_call_SuRe_lines_lenient_filtering.vcf.gz -delete
```

# 20191008

## Generate ref sequences per sample

Follow recipe from https://olgabotvinnik.com/blog/how-to-create-a-custom-genome-fasta-from-a-vcf-file/
```
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190912_DM6_Wolbachia_refSeq_Mattia
cat dm6.AE017196.fasta | awk '/Wolbachia/{OUT="wolbachia.fa";print >OUT; next} /^>chr/{OUT=substr($0,2) ".fa"}; OUT{print >OUT}'
for f in chr2L.fa chr2R.fa chr3L.fa chr3R.fa chr4.fa chrM.fa chrX*.fa chrY*.fa chrUn*.fa; do cat $f >> dm6.AE017196_karyoSort.fasta; done
rm -f chr*fa
samtools faidx dm6.AE017196_karyoSort.fasta
picard CreateSequenceDictionary R=dm6.AE017196_karyoSort.fasta O=dm6.AE017196_karyoSort.dict

# remove empty lines from VCF files:
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190909_VCF_fly_Mattia
for f in *filtering_*vcf; do mv $f tmp; awk NF tmp > $f; rm -f tmp; done

# create refseq fasta files

**BELOW CODE IS WRONG!!** As it turns out GATK creates genome sequences taking
the alternative allele at every SNP position, regardless haplotype. BCFtools
consensus can do what I want it to do. See dd20191015 for continued processing.

```
mkdir /DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/fasta

conda activate gatk
REFSEQ=/DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190912_DM6_Wolbachia_refSeq_Mattia/dm6.AE017196_karyoSort.fasta
declare -a SAMPLES=("B04" "DGRP-304" "DGRP-324" "DGRP-360" "DGRP-362" "DGRP-57" "DGRP-714" "I02" "I33" "N02" "T01" "ZH23" "vgn")

for sample in "${SAMPLES[@]}"; do
  echo "$sample"
  VAR=/DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190909_VCF_fly_Mattia/Haplotype_joint_call_SuRe_lines_lenient_filtering_${sample}_homoz.vcf
  OUT=/DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/fasta/${sample}_LP20191008.fa
  gatk3 -T FastaAlternateReferenceMaker -R ${REFSEQ} --variant ${VAR} -o ${OUT}
done
```

# 20191015

## Create refseq fastta, using bcftools consensus

```
# clean up previous data which are wrong
rm -rf /DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/fasta/*

conda activate vcftools
REFSEQ=/DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190912_DM6_Wolbachia_refSeq_Mattia/dm6.AE017196.fasta
declare -a SAMPLES=("B04" "DGRP-304" "DGRP-324" "DGRP-360" "DGRP-362" "DGRP-57" "DGRP-714" "I02" "I33" "N02" "T01" "ZH23" "vgn")

for sample in "${SAMPLES[@]}"; do
  echo "$sample"
  VAR=/DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190909_VCF_fly_Mattia/Haplotype_joint_call_SuRe_lines_lenient_filtering_${sample}_homoz.vcf
  OUT=/DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/fasta/${sample}_LP20191015.fa

  # compress vcf using bgzip/tabix
  bgzip ${VAR}
  tabix ${VAR}.gz

  CMD="bcftools consensus -H 1 -s ${sample}  -f ${REFSEQ} ${VAR}.gz -o ${OUT} 2> ${OUT%.fa}.err"
  echo ${CMD}
  eval ${CMD}
done

```


## Create bowtie2 indices

```
# set some paths to data files and runfiles
DATE=$(date +"%Y%m%d")
INDIR="/DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/fasta/"
OUTDIR="/DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/bowtie2-indices"
CHRS=( $(seq 1 22) 'X' )
BOWTIE=bowtie2-build
echo "BOWTIE VERSION:"
${BOWTIE} --version
echo -e "\n\n"

declare -a SAMPLES=("B04" "DGRP-304" "DGRP-324" "DGRP-360" "DGRP-362" "DGRP-57" "DGRP-714" "I02" "I33" "N02" "T01" "ZH23" "vgn")

NCORES=1
# loop pver samples and create all ref sequences
for sample in "${SAMPLES[@]}"; do
    indir="${INDIR}/${sample}"
    infile="${INDIR}/${sample}_LP20191015.fa"
    CMD=( "${BOWTIE}" --threads ${NCORES} "${infile}" "${OUTDIR}/${sample}/${sample}" ) 
    echo "CMD = ${CMD[@]}"
    mkdir -p "${OUTDIR}/${sample}"
    eval "${CMD[@]}"
done
```

# 20191028

## Generate chain and reverse chain files

```
cd /home/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_VCF_per_sample/

declare -a SAMPLES=("B04" "DGRP-304" "DGRP-324" "DGRP-360" "DGRP-362" "DGRP-57" "DGRP-714" "I02" "I33" "N02" "T01" "ZH23" "vgn")
REFSEQ=/DATA/usr/ludo/projects/LP190425_flySuRE/data/external/LP20190912_DM6_Wolbachia_refSeq_Mattia/dm6.AE017196.fasta
mkdir -p /DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_VCF_per_sample/chain/
for sample in "${SAMPLES[@]}"; do
  echo "$sample"
  VAR=/DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_VCF_per_sample/Haplotype_joint_call_SuRe_lines_lenient_filtering_${sample}_homoz.vcf.gz
  OUT=/DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_VCF_per_sample/chain/${sample}_LP20191028.chain

  CMD="bcftools consensus -H 1 -s ${sample} -f ${REFSEQ} ${VAR} -c ${OUT} 2> ${OUT%.chain}.err > /dev/null"
  conda activate vcftools
  echo "CMD = ${CMD[@]}"
  eval "${CMD[@]}"

  CMD="chainSwap <( grep -v \"^0 0 0\" ${OUT} ) ${OUT%.chain}_reverse.chain"
  conda activate base
  echo "CMD = ${CMD[@]}"
  eval "${CMD[@]}"
done
```

# 20191029 

## Overview of adapting the pipeline for current fly data

I adapted the pipeline from commit `84097a803e405a233a4b9ed25d930354581216e8` onward, in branch `fly-embl`.

The adaptations revolved around two aspects:

* adapting for input data, eg ref seq, VCF files, etc
* adapting for the homozygosity of the fly data

The adaptations start with commit `0aef3d3` and continue (at least) to commit `8bf8731`.

## Overview of workflow for processing data

The processing is done from directory
`../analyses/LP20191017_processing_cDNA_iPCR`. The entire project directory
(`../../`) is under git control, in the github repos
[LP20191016_flySuRE_project](https://github.com/vansteensellab-repos/LP20191016_flySuRE_project).

Wrt the processing I commit the command wrapper (eg `cmd_I02_LP20191028.sh`)
and the resulting log file(s) (eg
`config-Dm08_B04_LP20191028_run-LP20191028_0957.log`) to the repos.

The output is stored per fly strain in `../../data/intermediate/`, eg. `LP20191017_Dm12_T01_pipelineOutput`.

## Overview of current results.

The testing of the pipeline was done while processing the data of Dm12_T01. Now
that it is finished I checked the results, basically checking the counts for
iPCR and cDNA in the counts table
(`../../data/intermediate/LP20191017_Dm12_T01_pipelineOutput/count_tables/11_sorted/chr2R.bedpe.gz`).

I find that there are nearly no cDNA counts!!

The cDNA trimmed-table file suggests there are plenty of cDNAs, and the count
distribution looks good. Also the count distribution of the iPCR counts looks
good. It appears as if the cDNA counts in the final counts table is only due to
random matches. This could indicate that the samples are mixed...

### Checking whether (cDNA) samples are mixed

I want to parse all cDNA data, extract frequent barcode sequences from each,
and check overlap of each with the iPCR barcodes I have so far (T01, B04,
DGRP-304).

### Prep config and cmd wrappers for all strains

Here is an overview of strains and sample numbers:

| sample nr (based on github wiki) | strain short name | strain name in outputdir | strain long name in iPCR fastq| cDNA fastq fname | processing started | processing finished |
| --------- | ----------------- | --------------------- | ---------------- | --- |
| 03  | DGRP-324 | Dm03_DGRP-324  | DM03DGRP324  | SuRE_III_Dm3  |   |   |
| 04  | DGRP-360 | Dm04_DGRP-360  | DM04DGRP360  | SuRE_III_Dm4  |   |   |
| 05  | DGRP-362 | Dm05_DGRP-362  | DM05DGRP362  | SuRE_III_Dm5  |   |   |
| 06  | DGRP-714 | Dm06_DGRP-714  | DM06DGRP714  | SuRE_III_Dm6  |   |   |
| 08  | B04      | Dm08_B04       | DM08GDLB04   | SuRE_III_Dm8  | X |   |
| 09  | I02      | Dm09_I02       | DM09GDLI02   | SuRE_III_Dm9  | X |   |
| 10  | I33      | Dm10_I33       | DM10GDLI33   | SuRE_III_Dm10 |   |   |
| 11  | N02      | Dm11_N02       | DM11GDLN02   | SuRE_III_Dm11 |   |   |
| 12  | T01      | Dm12_T01       | DM12GDLT01   | SuRE_III_Dm12 | X |   |
| 13  | ZH23     | Dm13_ZH23      | DM13GDLZH23  | SuRE_III_Dm13 |   |   |
| **No libraries:** | 
| --  | DGRP-304 |   |   |   |
| --  | DGRP-57 |   |   |   |
| --  | vgn |   |   |   |

I created all remaining config files and command wrappers, among others with below code qbut also with lots of manual editing.




```
declare -a SAMPLES=("Dm03_DGRP-324" "Dm04_DGRP-360" "Dm05_DGRP-362" "Dm06_DGRP-714" "Dm10_I33" "Dm11_N02" "Dm13_ZH23")

pushd /DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191017_processing_cDNA_iPCR
for sample in "${SAMPLES[@]}"; do
  # copy template config file
  cp config-Dm08_B04_LP20191028.yml config-${sample}_LP20191029.yml
  # create new command wrapper from template, including replacing the correct config file
  sed 's/config-Dm08_B04_LP20191028.yml/config-'"${sample}"'_LP20191029.yml/g' cmd_B04_LP20191028.sh | sed 's/TARGET="bedpe"/TARGET="trim_cDNA"/g' - > cmd_${sample#Dm*_}_LP20191029.sh
done
popd
```

```
declare -a SAMPLES=("Dm03_DGRP-324" "Dm04_DGRP-360" "Dm05_DGRP-362" "Dm06_DGRP-714" "Dm10_I33" "Dm11_N02" "Dm13_ZH23")
for sample in "${SAMPLES[@]}"; do
  smpl=${sample#Dm*_}
  sed -i 's#T01: "chain/T01_LP20191028_reverse.chain"#'"${smpl}"': "chain/'"${smpl}"'_LP20191028_reverse.chain"#' config-${sample}_LP20191029.yml
done
```

# 20191031

## Sample *are* swapped

I finished processing all data. I checked the cDNA-iPCR sample correspondence in R (`../analyses/LP20191017_processing_cDNA_iPCR/check_cDNA-iPCR_SampleCorrelation_20191031.R', run in `../data/intermediate`) and found that only in 3 cases the barcodes of cDNA and iPCR data have a high overlap. For in total 6 strains there is a obvious correspondence between cDNA and iPCR samples (median overlap is 5e-5, high overlap >0.5, this includes the three correct samples). The remaining 4 strains have max overlap scores of:
```
0.02080964
0.06932100
0.09863326
0.33418182
```
Not so good.

Strains with lower overlap tend to have worse sequencing depth of the cDNA.

I sent an email to Matteo/Mattia with the bad news.

# 20191104

## Need for cross correlating iPCR data sets

Matteo replied he doesn't know what happened. He can't say whether either the
cDNA or the iPCR samples are swapped. So, I need to check the iPCR samples as
well.

For the iPCR data sets I want to extract reads for a single chromosome from the
bam files with previously aligned reads. Then re-align these reads to all other
9 genome sequences. I assume that the alignment against the correct genome
sequence will have the maximal alignment scores.

## Extract reads from chr3L, for T01

```
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/intermediate/
mkdir -p LP20191104_iPCR_cross_correlation
cd LP20191104_iPCR_cross_correlation

for f in ../LP2019*_pipelineOutput/iPCR/*_B1_T01/05_split_bam/*chr3L.bam; do
# for f in "../LP20191017_Dm12_T01_pipelineOutput/iPCR/T01_B1_T01/05_split_bam/*chr3L.bam"; do
  tt=$( basename ${f} )
  f1=${tt%.bam}_forw.fq
  f2=${tt%.bam}_rev.fq
  CMD="samtools bam2fq -1 ${f1} -2 ${f2} -F 2304 -@ 4 ${f}"
  echo "CMD = ${CMD[@]}"
  eval "${CMD[@]}"
  for b in /DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/bowtie2-indices/*; do
    index=${b}/$( basename ${b} )

    CMD="(bowtie2  --sam-no-qname-trunc -p 35 -x ${index} -1 ${f1} -2 ${f2} -X 500 | samtools view -S - | awk -F\"\t\" '{as=gensub(\"AS:i:(.*)\$\",\"\\\\1\",\"g\",\$12); print \$1\"\t\"as}' > ${tt%.bam}_$( basename $index ).tsv) 2>>${tt%.bam}.stats"
    echo "CMD = ${CMD[@]}"
    eval "${CMD[@]}"
  done
done
```

I made a mistake: The B04 sample used T01 as label. I renamed all files and dirs in `data/intermediate/LP20191028_Dm08_B04_pipelineOutput/iPCR/T01_B1_T01/05_split_bam/` from T01 to B04!!!

```
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/intermediate/
# mkdir -p LP20191104_iPCR_cross_correlation
cd LP20191104_iPCR_cross_correlation

for f in ../LP2019*_[BT]0[14]_pipelineOutput/iPCR/*_B1_T01/05_split_bam/*chr3L.bam; do
# for f in "../LP20191017_Dm12_T01_pipelineOutput/iPCR/T01_B1_T01/05_split_bam/*chr3L.bam"; do
  tt=$( basename ${f} )
  f1=${tt%.bam}_forw.fq
  f2=${tt%.bam}_rev.fq
  CMD="samtools bam2fq -1 ${f1} -2 ${f2} -F 2304 -@ 4 ${f}"
  echo "CMD = ${CMD[@]}"
  eval "${CMD[@]}"
  for b in /DATA/usr/ludo/projects/LP190425_flySuRE/data/processed/LP20191008_refseq_per_sample/bowtie2-indices/*; do
    index=${b}/$( basename ${b} )

    CMD="(bowtie2  --sam-no-qname-trunc -p 35 -x ${index} -1 ${f1} -2 ${f2} -X 500 | samtools view -S - | awk -F\"\t\" '{as=gensub(\"AS:i:(.*)\$\",\"\\\\1\",\"g\",\$12); print \$1\"\t\"as}' > ${tt%.bam}_$( basename $index ).tsv) 2>>${tt%.bam}.stats"
    echo "CMD = ${CMD[@]}"
    eval "${CMD[@]}"
  done
done
```

# 20191106

## Checking cross contamination of iPCR samples

Analysis of correspondence of iPCR samples shows that B04 and ZH23 have lower
correspondence to the appearent correct data files than other samples. One
possibility is that samples got cross contaminated. I want to check whether I
see a larger overlap of iPCR barcodes of the B04/ZH23 samples with other data
sets, as compared to the other samples. For this I will collect barcodes of the
iPCR data sets anddetermine overlap with other data sets, for all samples.

```
cd data/intermediate/LP20191106_iPCR_cross_contamination
for f in ../*_pipelineOutput/iPCR/*T01/01*/*for*info*; do 
  (strain=$(bn=$( basename $f ); 
   echo ${bn%_B1_T01_forw.info.gz}); 
   zcat $f | \
     awk -F"\t" '$5~/[ACGT]+/{print $5}' | \
     sort | \
     uniq -c | \
     sort -rg > ${strain}_BC.tbl& ) ; 
done

```

```
for t1 in *tbl; do 
  for t2 in *tbl; do 
    [ $t1 == $t2 ] && continue; 
    awk ' length($2)!=20{next} 
          NR==FNR{a[$2]=1; next} 
          {if($2 in a){a[$2]+=2}else{a[$2]=2}} 
          END{for(k in a){b[a[k]]++} 
              for(k in b){print k, b[k]}}' $t1 $t2 > ${t1%_BC.tbl}_${t2%_BC.tbl}.tbl; 
  done; 
done
```

## Fixing swaps

There are two types of swaps, and a possible contamination:
- cDNA samples 3,4,5,6,8 are swapped to 8,3,4,5,6
- iPCR samples 12, 13 are swapped to 13,12
- iPCR samples Dm08_B04 and Dm13_ZH23 (filename Dm12_T01) may have been cross-contaminated


| ------------ | ------------- | ------------- |
| sample       | iPCR filename | cDNA filename |
| ------------ | ------------- | ------------- |
| Dm03_DGRP324 | Dm03          | Dm04          |
| Dm04_DGRP360 | Dm04          | Dm05          |
| Dm05_DGRP362 | Dm05          | Dm06          |
| Dm06_DGRP714 | Dm06          | Dm08          |
| Dm08_B04     | Dm08          | Dm03          |
| Dm09_I02     | Dm09          | Dm09          |
| Dm10_I33     | Dm10          | Dm10          |
| Dm11_N02     | Dm11          | Dm11          |
| Dm12_T01     | Dm13          | Dm12          |
| Dm13_ZH23    | Dm12          | Dm13          |
| ------------ | ------------- | ------------- |

