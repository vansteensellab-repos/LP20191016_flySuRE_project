---
title: "Overview of pipeline results"
output: 
  html_document: 
    keep_md: yes
    toc: yes
---

## Init

```r
library(data.table)
library(parallel)
```

Note: In this document I'm using data.table partly as exercise. With hindsight for this kind of small work a data.frame is probably easier. 



## Introduction

This document will show some statistics on the iPCR and cDNA data and the results of processing these data with the SuRE-INDEL pipeline.

The results in this document are taking into account the following sample - iPCR datafile - cDNA datafile correspondence into account:

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

### Association strains and filenames


```r
assoc <- data.table(
  strain=c("Dm03_DGRP324", "Dm04_DGRP360", "Dm05_DGRP362", "Dm06_DGRP714", "Dm08_B04", "Dm09_I02", "Dm10_I33", "Dm11_N02", "Dm12_T01", "Dm13_ZH23"),
  ipcr=c("Dm03", "Dm04", "Dm05", "Dm06", "Dm08", "Dm09", "Dm10", "Dm11", "Dm13", "Dm12"),
  cdna=c("Dm04", "Dm05", "Dm06", "Dm08", "Dm03", "Dm09", "Dm10", "Dm11", "Dm12", "Dm13"),
  stringsAsFactors = FALSE)
assoc <- cbind(assoc, short=sub("(Dm..).*","\\1",assoc[["strain"]]))
setkey(assoc, short)
assoc
```

```
##           strain ipcr cdna short
##  1: Dm03_DGRP324 Dm03 Dm04  Dm03
##  2: Dm04_DGRP360 Dm04 Dm05  Dm04
##  3: Dm05_DGRP362 Dm05 Dm06  Dm05
##  4: Dm06_DGRP714 Dm06 Dm08  Dm06
##  5:     Dm08_B04 Dm08 Dm03  Dm08
##  6:     Dm09_I02 Dm09 Dm09  Dm09
##  7:     Dm10_I33 Dm10 Dm10  Dm10
##  8:     Dm11_N02 Dm11 Dm11  Dm11
##  9:     Dm12_T01 Dm13 Dm12  Dm12
## 10:    Dm13_ZH23 Dm12 Dm13  Dm13
```

## Raw read counts of iPCR and cDNA fastq files

Statistics of the read counts are collected from the multiQC output (ie overview of the fastQC reports).

### iPCR

```r
# read ipcr data
ipcrStats <- fread('/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20190909_resequencing-I_1st_samples/iPCR_resequencing_I_multiQC/multiqc_data/multiqc_general_stats.txt')
# short name based on filename
ipcrStats$shortname <- paste0(rep(each=2, sprintf("Dm%0.2d",c(3,4,5,6,8,9,10,11,12,13))), c("","_2"))
setkey(ipcrStats, shortname)
# rownames(ipcrStats) <- paste0(ipcrStats$shortname, c("","_2"))
# readcount scaled to million
ipcrStats[,6] = signif(ipcrStats[,6]/1e6, digits = 3)
# rename readcount column
names(ipcrStats)[6]="totalReads (M)"
# show correspondence of given names and shortnames
owidth <- options(width=500)
ipcrStats[,c(1,7)]
```

```
##                                                                   Sample shortname
##  1: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM03DGRP324_1_sequence      Dm03
##  2: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM03DGRP324_2_sequence    Dm03_2
##  3: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM04DGRP360_1_sequence      Dm04
##  4: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM04DGRP360_2_sequence    Dm04_2
##  5: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM05DGRP362_1_sequence      Dm05
##  6: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM05DGRP362_2_sequence    Dm05_2
##  7: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM06DGRP714_1_sequence      Dm06
##  8: H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM06DGRP714_2_sequence    Dm06_2
##  9:  H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM08GDLB04_1_sequence      Dm08
## 10:  H3C3HBGXC_iPCR_mix1_19s003538-1-1_Perino_lane1DM08GDLB04_2_sequence    Dm08_2
## 11:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM09GDLI02_1_sequence      Dm09
## 12:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM09GDLI02_2_sequence    Dm09_2
## 13:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM10GDLI33_1_sequence      Dm10
## 14:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM10GDLI33_2_sequence    Dm10_2
## 15:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM11GDLN02_1_sequence      Dm11
## 16:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM11GDLN02_2_sequence    Dm11_2
## 17:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM12GDLT01_1_sequence      Dm12
## 18:  H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM12GDLT01_2_sequence    Dm12_2
## 19: H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM13GDLZH23_1_sequence      Dm13
## 20: H5JNHBGXC_iPCR_mix2_19s003539-1-1_Perino_lane1DM13GDLZH23_2_sequence    Dm13_2
```

```r
# options(owidth)
# print short name/readcount columns
ipcrStats[seq(1,20, by=2),c(7,6)]
```

```
##     shortname totalReads (M)
##  1:      Dm03           71.1
##  2:      Dm04           52.3
##  3:      Dm05           66.9
##  4:      Dm06           74.6
##  5:      Dm08           70.3
##  6:      Dm09           87.3
##  7:      Dm10           75.1
##  8:      Dm11           70.5
##  9:      Dm12           78.1
## 10:      Dm13           79.3
```
### cDNA

```r
# read cdna data
cdnaStats <- fread('/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20190909_resequencing-I_1st_samples/cDNA_resequencing_I_multiQC/multiqc_data/multiqc_general_stats.txt')
# create short sample name
cdnaStats$shortname <- sprintf("Dm%0.2d",c(10,11,12,13,3,4,5,6,8,9))
setkey(cdnaStats, shortname)
# show correspondence of given names and shortnames
cdnaStats[,c(1,7)]
```

```
##                     Sample shortname
##  1:  SuRE_III_Dm3_sequence      Dm03
##  2:  SuRE_III_Dm4_sequence      Dm04
##  3:  SuRE_III_Dm5_sequence      Dm05
##  4:  SuRE_III_Dm6_sequence      Dm06
##  5:  SuRE_III_Dm8_sequence      Dm08
##  6:  SuRE_III_Dm9_sequence      Dm09
##  7: SuRE_III_Dm10_sequence      Dm10
##  8: SuRE_III_Dm11_sequence      Dm11
##  9: SuRE_III_Dm12_sequence      Dm12
## 10: SuRE_III_Dm13_sequence      Dm13
```

```r
# readcounts scaled to million
cdnaStats[,6] <- signif(cdnaStats[,6]/1e6, digits=3)
# rename readcount column
names(cdnaStats)[6] <- "totalReads (M)"
# print sort name/readcount columns
cdnaStats[,c(7,6)]
```

```
##     shortname totalReads (M)
##  1:      Dm03          34.60
##  2:      Dm04          25.90
##  3:      Dm05          20.60
##  4:      Dm06          26.90
##  5:      Dm08          33.10
##  6:      Dm09          34.40
##  7:      Dm10          33.50
##  8:      Dm11          33.00
##  9:      Dm12           5.12
## 10:      Dm13          21.80
```

## Combine iPCR/cDNA counts, using correct data sets per strain


```r
stats <- data.table(libComplFlat=0, libCompl=0,cdnaComplFlat=0,cdnaCompl=0,iPCRraw=0, cDNAraw=0)[-1,]
for (i in seq.int(nrow(assoc))) {
  stats = rbind(stats,list(0,0,0,0,0,0))
  set(stats, i=i, j="iPCRraw", value = ipcrStats[assoc[i,.(ipcr)][[1]],.(`totalReads (M)`)][[1]])
  set(stats, i=i, j="cDNAraw", value = cdnaStats[assoc[i,.(cdna)][[1]],.(`totalReads (M)`)][[1]])
}
stats <- cbind(stats, short=assoc$short,long=assoc$strain)
setkey(stats, short)
stats
```

```
##     libComplFlat libCompl cdnaComplFlat cdnaCompl iPCRraw cDNAraw short         long
##  1:            0        0             0         0    71.1   25.90  Dm03 Dm03_DGRP324
##  2:            0        0             0         0    52.3   20.60  Dm04 Dm04_DGRP360
##  3:            0        0             0         0    66.9   26.90  Dm05 Dm05_DGRP362
##  4:            0        0             0         0    74.6   33.10  Dm06 Dm06_DGRP714
##  5:            0        0             0         0    70.3   34.60  Dm08     Dm08_B04
##  6:            0        0             0         0    87.3   34.40  Dm09     Dm09_I02
##  7:            0        0             0         0    75.1   33.50  Dm10     Dm10_I33
##  8:            0        0             0         0    70.5   33.00  Dm11     Dm11_N02
##  9:            0        0             0         0    79.3    5.12  Dm12     Dm12_T01
## 10:            0        0             0         0    78.1   21.80  Dm13    Dm13_ZH23
```

## Add iPCR and cDNA processing statistics


```r
fnames <- grep(value=TRUE, pattern="LP20191128.*11_sorted",dir(path='/DATA/usr/ludo/projects/LP190425_flySuRE/data/intermediate/',recursive=T, include.dirs=T,pattern=".*gz",full.names=T))
fname.sample <- sub("^.*LP20191128_(Dm..).*_pipelineOutp.*$","\\1",fnames)
samples <- unique(fname.sample)

data4viz <- list()
for (s in samples) {
  fs <- fnames[fname.sample %in% s]
  dt <- rbindlist(lapply(fs, fread, select=c(2,3,4,5,14,15), sep="\t"))
  sname <- setdiff(colnames(dt), c("BC","chrom","start","end","strand","count"))
  stats[s,c('libCompl', 'cdnaCompl', 'libComplFlat', 'cdnaComplFlat') := 
          dt[,.('libCompl'=sum(count), 'cdnaCompl'=sum(get(sname)), 
                'libComplFlat'=sum(count>0), 'cdnaComplFlat'=sum(0+(get(sname)>0)))]]

  data4viz[[s]] <- list(iPCR=dt[,.N,by="count"][order(count)],
                        cDNA=dt[,.N,by=eval(sname)][order(get(sname))],
                        ratio=dt[get(sname)>0,.(ratio=get(sname)/count)][,.N, by=ratio][order(ratio)],
                        chrom=dt[,.N,by=chrom][order(chrom)],
                        fragLen=dt[!(is.na(start) | is.na(end)),.(len=end-start+1)][,.N, by=len][order(len)],
                        sampledCompl=list(sizes=seq(1e6, dt[,sum(count)], length.out = 10),
                                             compl=sapply(seq(1e6, dt[,sum(count)], length.out = 10),
                                                          function(s){length(unique(sample(rep(dt[,.I], dt$count), size=s)))})))
}
```

```
## Registered S3 method overwritten by 'R.oo':
##   method        from       
##   throw.default R.methodsS3
```

## Visualization

### iPCR counts and SuRE fragments


```r
#for (type in c("iPCR","cDNA","ratio","fragLen")) {
col <- rep(1:5,2)
names(col) <- names(data4viz)
pch <- rep(c(4,19), each=5)
names(pch) <- names(data4viz)
type = "iPCR"
xlim <- range(unlist(lapply(data4viz, function(e) {e[[type]]$count})))
ylim <- range(unlist(lapply(data4viz, function(e) {e[[type]]$N})))
# ylim <- c(1, 7e7)
plot(NA,ty='p', log='yx', pch=19, cex=.5, xlim=xlim, ylim=ylim, 
     ylab='frequency',xlab='iPCR count', main=sprintf("Distribution iPCR counts per SuRE fragment"))
for (s in names(data4viz)) {
  points(data4viz[[s]][[type]], pch=pch[s], cex=.5, col=col[s])
}
# abline(h=ipcrStats$`totalReads (M)`[c(T,F)]*1e6, lty=3, col=1:5)
legend(x='topr', inset=0.02, bty='n',legend=names(data4viz), col=col, pch=pch,ncol = 2)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plot.iPCR.distr-1.png" width="1000" />


```r
col <- grey(c(0,.45,.95))
barplot(rbind(ipcrStats$`totalReads (M)`[c(T,F)]*1e6,
              sapply(data4viz, function(e) e[["iPCR"]][,sum(N*count)]),
              sapply(data4viz, function(e) e[["iPCR"]][,sum(N)])), 
        beside = T, main="iPCR: raw read count, sum of iPCR counts, and unique fragment count", ylab="count", ylim=c(0,1e8),las=2,col=col)
legend(x='topl', inset=0.025, bty='n', legend=c("raw read count","aligned read count","unique fragment count"), fill=col)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plot.iPCR.depth-1.png" width="1000" />

**CONCLUSION:**  
All samples have a fragment complexity 15-20e6 fragments. This is according to expectations. Also see below.

### cDNA counts


```r
col <- rep(1:5,2)
names(col) <- names(data4viz)
pch <- rep(c(4,19), each=5)
names(pch) <- names(data4viz)
type = "cDNA"
xlim <- range(unlist(lapply(data4viz, function(e) {e[[type]][[1]]})))
ylim <- range(unlist(lapply(data4viz, function(e) {e[[type]]$N})))
plot(NA,ty='p', log='y', pch=19, cex=.5, xlim=xlim, ylim=ylim, 
     ylab='frequency',xlab='cDNA count', main=sprintf("Distribution cDNA counts per SuRE fragment"))
for (s in names(data4viz)) {
  points(data4viz[[s]][[type]], pch=pch[s], cex=.5, col=col[s])
}
legend(x='topr', inset=0.02, bty='n',legend=names(data4viz), col=col, pch=pch,ncol = 2)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plot.cDNA.distr-1.png" width="1200" />

```r
par(mfrow=c(2,5), mar=c(2,2,1,1))
for (s in names(data4viz)) {
  plot(data4viz[[s]][[type]], pch=pch[s], cex=.5, col=col[s],ty='p', log='y', ylim=ylim, 
     ylab='',xlab='', main=sprintf("%s",s))
}
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plot.cDNA.distr-2.png" width="1200" />


**CONCLUSION:**  
Several cDNA samples appear over-amplified; *Dm03, Dm04, Dm05, Dm13*, as they show increasing count frequencies for the lowest counts. Samples *Dm06 - Dm11* are more or less monotenously decreasing with increasing counts.  
*Dm12* appears extremely undersequenced, but also strongly over-amplified.

Below is the number of barcodes (BC's) present in the cDNA data, and the number of these BC's actually seen in the iPCR data. The third bar depicts the proportion of cDNA BC's seen in iPCR (axis on right side).


```bash
cd /DATA/usr/ludo/projects/LP190425_flySuRE/data/intermediate
for d in LP20191128*ut; do 
awk -v sample="$d" ' 
  FNR==NR{a[$2]=1; next} 
  {if($1 in a){c++}} 
  END{print sample"\t"c"\t"length(a)}' \
  <( zcat $d/cDNA/*_B1_T1/*_B1_T1_trimmed_table.txt.gz) \
  <( for f in $d/count_tables*/09_ipcr_cdna_merged/ch*.bedpe.gz ; do zcat $f | awk 'NR>1'; done ); done > /tmp/tt
```


```r
opar <- par(mar=c(5,4,4,4)+.1)
col <- grey(c(0.25,0.45,0.8))
cDNAseen <- read.table("/tmp/tt",header=F, col.names = c('sample','seen','total'))
rownames(cDNAseen) <- sub("LP20191128_(Dm..)_.*_pipelineOutput","\\1",cDNAseen$sample)
barplot(t(as.matrix(cbind(cDNAseen[,3:2],prop=(cDNAseen[[2]]/cDNAseen[[3]])*max(cDNAseen$total)))),beside=TRUE, ylab="BC count",
        density=c(-1,-1,10),col=col, main="cDNA barcodes (BC) in sample and present in iPCR", ylim=c(0,max(cDNAseen$total)*1.35), las=2)
axis(4,at=seq(0,max(cDNAseen$total),length.out = 6),labels=seq(0,1,length.out = 6))
mtext("proportion seen",side=4,line=3)
par(opar)
legend('topl',bty='n',inset=0.025,legend=c("number unique BC in cDNA","number BC seen in iPCR","proportion of BC seen in iPCR"),
       fill=col, density=c(-1,-1,30),xpd=NA)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/cdnacounts-1.png" width="100%" />

In accordance with the undersampling/over-amplification of samples *Dm03, Dm04, Dm05, Dm12, Dm13* these samples (except *Dm13*) also exhibit  low proportion of cDNA BC's which are seen in the iPCR. This is probably due to the overrepresentation of a small number of BC's in these samples.

### Length of SuRE fragments


```r
col <- rep(1:5,2)
names(col) <- names(data4viz)
pch <- rep(c(4,19), each=5)
names(pch) <- names(data4viz)
type = "fragLen"
xlim <- c(1,max(unlist(lapply(data4viz, function(e) {e[[type]][[1]]}))))
ylim <- c(1,max(unlist(lapply(data4viz, function(e) {e[[type]]$N}))))
main <- "Distribution fragment length"
xlab=sprintf("fragment length (bp)")
ylab <- "frequency"
plot(NA,ty='p', log='', xlim=xlim, ylim=ylim, ylab=ylab,xlab=xlab, main=main)
for (s in names(data4viz)) {
  points(data4viz[[s]][[type]], pch=pch[s], cex=.5, col=col[s])
}
legend(x='topr', inset=0.02, bty='n',legend=names(data4viz), col=col, pch=pch,ncol = 2)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plot.fraglen.distr-1.png" width="1200" />

**CONCLUSION:**  
The distribution of the length of the SuRE fragments looks good. It is unclear what causes the dip at length ~ 140bp.

### Distribution of SuRE scores (cDNA_count/iPCR_ count, per SuRE fragment)


```r
col <- rep(1:5,2)
names(col) <- names(data4viz)
pch <- rep(c(4,19), each=5)
names(pch) <- names(data4viz)
type = "ratio"
xlim <- range(unlist(lapply(data4viz, function(e) {e[[type]][[1]]})))
ylim <- c(0,100) #range(unlist(lapply(data4viz, function(e) {e[[type]]$N})))
main <- sprintf("distribution cDNA/iPCR ratio")
xlab <- "ratio(cDNA/iPCR), for cDNA>0"
ylab <- "cumulative percentage"
plot(NA,ty='p', log='x', xlim=xlim, ylim=ylim, ylab=ylab,xlab=xlab, main=main)
for (s in names(data4viz)) {
  with(data4viz[[s]][[type]], points(ratio, 100*cumsum(N)/sum(N), pch=pch[s], cex=.5, col=col[s]))
}
legend(x='topl', inset=0.02, bty='n',legend=names(data4viz), col=col, pch=pch,ncol = 2)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plot.SuRE.score-1.png" width="1200" />



```r
chromcount <- sapply(data4viz, function(e){e[["chrom"]][["N"]]})
rownames(chromcount) <- data4viz[[1]][["chrom"]][["chrom"]]
chromcount <- t(t(chromcount)/colSums(chromcount))
barplot(t(chromcount), beside = T, col=1:5, angle=rep(c(45,-45),times=5), density=10,space = c(.0,1.2))
legend(x='topr',inset=0.025, bty='n', legend=names(data4viz), fill =1:5, angle=rep(c(45,-45),times=5), density=20,cex=1.2)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plot.chrom-1.png" width="1200" />


**CONCLUSION:**  
The over-amplified cDNA samples have relatively high cDNA/iPCR ratio's. The discrete nature of the ratio scores of the other samples is due to the low counts

### iPCR and cDNA complexity


```r
col <- rep(1:5, 2)
pch <- rep(c(4,19), each=5)
plot(stats$iPCRraw, stats$libCompl, col=col, pch=19,log='y', xlab="Raw iPCR read count (in fastq, x1e6)", ylab="SuRE fragment count (x1e6, see legend)", main="Library complexity and iPCR counts vs raw read count",ylim=range(stats$libCompl,stats$libComplFlat))
legend(x="topl", inset = 0.025, bty='n', legend=c("fragment total count","library complexity\n(ie, detected fragment)"), pch=c(19,1))
text(stats$iPCRraw, stats$libCompl,labels=stats$short,pos=c(2,3,1,4),xpd=NA,col=col)
points(stats$iPCRraw, stats$libComplFlat, col=col, pch=1)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plots-1.png" width="100%" />

```r
plot(stats$libCompl, stats$libComplFlat, pch=19, xlab="sum of fragments", ylab="library complexity (unique fragments)", main="iPCR library complexity; fragment counts vs unique fragments (lib complexity)", col=col)
#abline(lm(libComplFlat~libCompl,data = as.data.frame(stats)), lty=3)
text(stats$libCompl, stats$libComplFlat,labels=stats$short,pos=4:1,xpd=NA,col=col)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plots-2.png" width="100%" />

```r
cor(stats$libCompl, stats$libComplFlat)
```

```
## [1] 0.1143634
```

```r
plot(stats$cdnaCompl, stats$cdnaComplFlat, pch=19, xlab="sum of BC counts", ylab="cDNA complexity (unique BC's)", main="cDNA complexity; unique BCs vs total BC count",col=col, ylim=c(-5e5,0.7e7))
#abline(lm(cdnaComplFlat~cdnaCompl,data = as.data.frame(stats)), lty=3)
text(stats$cdnaCompl, stats$cdnaComplFlat,labels=stats$short,pos=1:4,xpd=NA,col=col)
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/plots-3.png" width="100%" />

```r
cor(stats$cdnaCompl, stats$cdnaComplFlat)
```

```
## [1] 0.8444158
```

**Conclusion:**  

* Relation of iPCR data *raw read count* and *total frament count* shows a good correlation.
* *Total  fragment count* and *unique fragment count* is poorly correlated
* *Dm04* and *Dm05* in particular have very low library complexity
* *cDNA BC complexity* correlates well with *cDNA BC counts *
* The under sequenced and over-sampled cDNA samples (*Dm03, Dm04, Dm05, Dm12, Dm13*) have much lower BC complexity than the other samples.


## Estimated saturation of iPCR library complexity

Below estimates use a method developed by Tao Chen in our group. The method is not documented but appears to be based on fitting a negative binomial distribution on the observed complexity in the iPCR data, and in samplings of 6this data in decreasing sample sizes. The iPCR sample complexity is estimated using the parameterization of the fitted distribution.


```r
estimComplex <- function(reads, complexity, plot=FALSE, sample="sample") {
  idx <- 4:10
  #cow = reads[idx]/min(reads[idx]);
  reads = reads[idx]/min(reads[idx]);
  #bull = complexity[idx];
  complexity = complexity[idx];
  uni = vector();
  
  for (i in 1:(length(reads)-1)){
    for (j in (i+1):length(reads)){
      t1=reads[i];
      t2=reads[j];
      
      approx <- function(x,t1=1,t2=2,b=0) {(1-x^t1)/(1-x^t2) - b}
      xirr <- function() {
        tryCatch(
          {irr <- uniroot(f=approx,c(0,3),  b=complexity[i]/complexity[j], t1=reads[i], t2=reads[j])[["root"]];
          return(irr)}, 
          error=function(err){return(NA)}
        )
        return(irr)
      }
      root <- xirr()
      uni <- c(uni,root);
    }
  }

  ttl=mean(complexity/(1-mean(uni)^reads))
  if(plot){
    plot(seq(0,2*max(reads*ttl), length.out=20), ttl*(1-mean(uni)^seq(0,2*max(reads), length.out=20)),pch=19,type='o',
         main =paste(sample,': saturation library complexity\nk =',round(log(mean(uni)),2),'; complexity =',round(ttl)),ylim=c(0,ttl),
         xlab="total fragment count", ylab="complexity (unique fragments)")
    points(c(0,reads*ttl),c(0,complexity),pch=19,col='red');
    legend("bottomr", inset=0.025, bty="n",legend=c("modeled complexity","measured complexity","estimated SuRE library complexity"),col=c("black","red","grey"), 
           pch=c(19,19,NA),lty=c(NA,NA,3))
    abline(h=ttl,lty=2,col='grey')
  }
  return(round(ttl))
}
for(s in names(data4viz)) {
  estimComplex(data4viz[[s]][[6]][[1]],data4viz[[s]][[6]][[2]],plot=TRUE,sample=stats[s]$long)
}
```

<img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-1.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-2.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-3.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-4.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-5.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-6.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-7.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-8.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-9.png" width="100%" /><img src="/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191108_review_results_processing_LP20191106/review_results_processing/review_results_processing_files/figure-html/estimCompl-10.png" width="100%" />

**CONCLUSION:**  
Most libraries are sequenced to near saturation (according to above estimations). Library *Dm04* is most undersequenced, and least saturated. An additional 
# R session


```r
sessionInfo()
```

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.6 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/libblas/libblas.so.3.6.0
## LAPACK: /usr/lib/lapack/liblapack.so.3.6.0
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] data.table_1.12.2
## 
## loaded via a namespace (and not attached):
##  [1] compiler_3.6.2   magrittr_1.5     tools_3.6.2      htmltools_0.3.6  yaml_2.2.0       Rcpp_1.0.1       codetools_0.2-16 stringi_1.4.3    rmarkdown_1.12   knitr_1.22       stringr_1.4.0    xfun_0.6         digest_0.6.18    evaluate_0.13
```
