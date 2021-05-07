library(data.table)
getwd()
bedpe_fnames=dir("07_bedpe_BC/",pattern="gz", full.names=T)
col.names.select <- c("BC", "chrom","start","end","strand","MAPQ1","MAPQ2","AS1","AS2","XS1","XS2","SNP_ABS_POS","SNP_REL_POS","SNP_ID","SNP_SEQ","SNP_VAR","SNP_PARENT","SNP_TYPE","SNP_SUBTYPE")
bedpe <- rbindlist(lapply(bedpe_fnames, fread, nrows=-100, nThread=4, key=c("BC","chrom","start","end"), col.names=col.names.select, select=col.names.select))
setkey(bedpe, BC, chrom, start, end)
col.names.select <- c("BC", "chrom","start","end","strand","MAPQ1","MAPQ2","AS1","AS2","XS1","XS2","SNP_ABS_POS","SNP_REL_POS","SNP_ID","SNP_SEQ","SNP_VAR","SNP_PARENT","SNP_TYPE","SNP_SUBTYPE")
Sys.time(); bedpe.1128 <- bedpe[,':='(MAPQ=MAPQ1+MAPQ2, AS=AS1+AS2)][,':='(n=.N),by=.(BC,start,end)][,.SD[n==max(n), .SD[MAPQ==max(MAPQ), .SD[AS==max(AS), .SD[1], by=AS], by=MAPQ], by=n],by=BC]; Sys.time()
names(bedpe)
bedpe <- rbindlist(lapply(bedpe_fnames, fread, nrows=-100, nThread=4, key=c("BC","chrom","start","end"), col.names=col.names.select, select=col.names.select))
bedpe.org <- copy(bedpe)
Sys.time(); bedpe <- copy(bedpe.org); setkey(bedpe, BC, chrom, start, end); qq <- bedpe[,':='(MAPQ=MAPQ1+MAPQ2, AS=AS1+AS2, rnd=sample.int(.N))][,':='(n=.N,mapqmaxidx=(MAPQ==max(MAPQ))),by=.(BC,start,end)][,':='(nmaxidx=(n==max(n))),by=.(BC)][,':='(asmaxidx=(AS==max(AS))),by=.(BC,start,end,MAPQ)][,.SD[asmaxidx&mapqmaxidx&nmaxidx][rnd==max(rnd)], by=BC]; qq[1:3]; Sys.time() 
str(bedpe)
q)
nrow(q)
nrow(qq)
sessionInfo()
ls()
bedpe.10000 <- copy(bedpe.org)
bedpe.10000 <- bedpe[1:10000]
save(file="/tmp/bedpe.10000.rda",bedpe.10000)
names(bedpe.10000)
names(bedpe.org)
bedpe.10000 <- copy(bedpe.org[1:10000])
save(file="/tmp/bedpe.10000.rda",bedpe.10000)
sindri2 <- function(ttt) {
ttt[ttt[, .N, by = .(BC, chrom,start,end)][, .SD[N == max(N)], by = BC],
    on = .(BC, chrom,start,end)
    ][, .SD[MAPQ1+MAPQ2 == max(MAPQ1+MAPQ2)], by = BC
    ][, .SD[AS1+AS2 == max(AS1+AS2)], by = BC
    ][, .SD[ifelse(N==1,1,sample(.N,1))], by = BC
    ][,!"N"]   
}              
Sys.time(); bedpe <- copy(bedpe.org); sindri2(bedpe); Sys.time()
sindri2
getwd()
setwd("/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191128_processing_cDNA_iPCR_redownload")
savehistory("Rhistory_merge-bedpe_LP20191208.R")
