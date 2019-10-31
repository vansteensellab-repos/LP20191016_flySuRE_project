library(data.table)
setwd('..')
dir(recursive=TRUE, pattern=".bedpe.gz")
dir(recursive=TRUE, pattern=".bedpe.gz",path="11_sorted")
dir(recursive=TRUE, pattern=".bedpe.gz",path="LP20191*_pipelineOutput/count_tables/11_sorted/")
dir(recursive=TRUE, pattern=".bedpe.gz",path="LP20191.*_pipelineOutput/count_tables/11_sorted/")
dir(recursive=TRUE, pattern=".bedpe.gz",path="./LP20191.*_pipelineOutput/count_tables/11_sorted/")
dir(recursive=TRUE, pattern=".bedpe.gz",path="./")
dir(recursive=TRUE, pattern=".*/11_sorted/*.bedpe.gz",path="./")
dir(recursive=TRUE, pattern="*/11_sorted/*.bedpe.gz")
dir(recursive=TRUE, pattern="*/count_tables/11_sorted/*.bedpe.gz")
dir(recursive=TRUE, pattern=".bedpe.gz")
grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"))
grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T)
dirname(grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T))
grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T)
sub("LP20191029_Dm04_","",grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T))
sub("LP2019102._Dm.._","",grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T))
sub("LP201910.._Dm.._","",grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T))
sub("LP201910.._Dm.._(.*)_pipelineOutput.*","",grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T))
sub("LP201910.._Dm.._(.*)_pipelineOutput.*","\\1",grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T))
unique(sub("LP201910.._Dm.._(.*)_pipelineOutput.*","\\1",grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T)))
for(f in grep("11_sort",dir(recursive=TRUE, pattern=".bedpe.gz"), value=T)){ data[[f]]=fread(f)}
data=list()
grep("cDNA",dir(recursive=TRUE, pattern="_B1_T1_trimmed_table.txt.gz"), value=T)
for(f in grep("cDNA",dir(recursive=TRUE, pattern="_B1_T1_trimmed_table.txt.gz"), value=T)){data[[f]]=fread(f)}
length(data)
lapply(data,head)
iPCR=list()
for(f in grep("iPCR",dir(recursive=TRUE, pattern="forw.info.gz"), value=T)){iPCR[[f]]=fread(f)}
iPCR[[1]]
length(iPCR)
f
fread(f,nrow=10)
fread(f,nrow=10,select=c(5),sep="\t",skip=9)
for(f in grep("iPCR",dir(recursive=TRUE, pattern="forw.info.gz"), value=T)){iPCR[[f]]=fread(f,select=c(5),sep="\t",skip=9)}
for(f in grep("iPCR",dir(recursive=TRUE, pattern="forw.info.gz"), value=T)){iPCR[[f]]=fread(f,select=c(5),sep="\t",skip=9,fill=TRUE)}
iPCR[[1]]
iPCR[[1]][,.N]
iPCR[[1]][,.N, by=V5]
iPCR[[1]][1:10,.N, by=V5]
iPCR[[1]][,.N, by=V5]
iPCR[[1]][1:10000,count=.N, by=V5]
iPCR[[1]][1:10000,.(count=.N), by=V5]
iPCR.tbl = lapply(iPCR, function(e) e[,.(count=.N),by=V5])
iPCR[[1]]
iPCR.tbl~[[1]]
iPCR.tbl[[1]]
names(iPCR)
names(cDNA)
names(data)
names(iPCR.tbl)
data[[1]]
all = data.table(BC=union(iPCR.tbl[[1]][,V5], cDNA[[1]][,V2]))
all = data.table(BC=union(iPCR.tbl[[1]][,V5], data[[1]][,V2]))
all
all = data.table(BC=uniq(c(iPCR.tbl[[1]][,V5], data[[1]][,V2])))
all = data.table(BC=unique(c(iPCR.tbl[[1]][,V5], data[[1]][,V2])))
all
nrow(iPCR.tbl[[1]]) + nrow(data[[1]])
head(data[[1]][[V2]])
head(data[[1]][["V2"]])
mean(data[[1]][["V2"]] %in% iPCR.tbl[[1]][["BC"]])
head(data[[1]][["V2"]])
head(iPCR.tbl[[1]][["BC"]])
head(iPCR.tbl[[1]][["V5"]])
mean(data[[1]][["V2"]] %in% iPCR.tbl[[1]][["V5"]])
lapply(data,function(e) mean(e[["V2"]] %in% iPCR.tbl[[1]][["V5"]])
)
lapply(data,function(e) mean(e[["V2"]] %in% iPCR.tbl[[2]][["V5"]]))
?mapply
res=lapply(data,function(cdna) lapply(iPCR.tbl, function(ipcr) mean(cdna[["V2"]] %in% ipcr[["V5"]]))
res=lapply(data,function(cdna) lapply(iPCR.tbl, function(ipcr) mean(cdna[["V2"]] %in% ipcr[["V5"]])))
res
as.data.frame(res)
str(as.data.frame(res))
dim(as.data.frame(res))
sapply(res,length)
sapply(res,unlist)
as.data.frame(sapply(res,unlist))
tt=as.data.frame(sapply(res,unlist))
dimnames(tt)
lapply(dimnames(tt), sub, ".*(Dm..).*","\\1")
lapply(dimnames(tt), function(x) sub(".*(Dm..).*","\\1",x))
dimnames(tt)=lapply(dimnames(tt), function(x) sub(".*(Dm..).*","\\1",x))
tt
signif(tt,digits=2)
tt
signif(tt,digits=2)

signif(tt,digits=1)
plot(1)
heatmap(tt)
heatmap(as.matrix(tt))
?heatmap
heatmap(as.matrix(tt),Rowv=NA,Colv=NA)
heatmap(as.matrix(tt),Rowv=NA,Colv=NA,xlab="bla")
dimnames(as.data.frame(sapply(res,unlist)))
colnames(as.data.frame(sapply(res,unlist)))
heatmap(as.matrix(tt),Rowv=NA,Colv=NA,xlab="cDNA", ylab="iPCR")
heatmap(as.matrix(tt),Rowv=NA,Colv=NA,xlab="cDNA", ylab="iPCR",scale=NA)
heatmap(as.matrix(tt),Rowv=NA,Colv=NA,xlab="cDNA", ylab="iPCR",scale=NULL)
heatmap(as.matrix(tt),Rowv=NA,Colv=NA,xlab="cDNA", ylab="iPCR",scale="none")
heatmap(log(as.matrix(tt)),Rowv=NA,Colv=NA,xlab="cDNA", ylab="iPCR",scale="none")
heatmap(as.matrix(tt),Rowv=NA,Colv=NA,xlab="cDNA", ylab="iPCR",scale="none")
signif(tt,digits=1)
savehistory("check_cDNA-iPCR_SampleCorrelation_20191031.R")
