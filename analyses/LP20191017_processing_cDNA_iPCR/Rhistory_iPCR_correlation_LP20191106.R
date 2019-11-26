fread("T01.correl"")
fread("T01.correl")
library(data.tab)
library(data.table)
fread("T01.correl")
fread("*.correl")
dir(pattern="correl")
fread(dir(pattern="correl"))
data=list()
for (f in dir (pattern="correl")){data[[f]]=fread(f)}
data
for (f in dir (pattern="correl")){data[[f]]=fread(f,col.names=c("sample","count"),key='sample')}
data
names(data)
sub(".correl","",names(data))
names(data) = sub(".correl","",names(data))
data[[1]][names(data),]
foo=function(dt){m=mean(dt$counts); s=stdev(dt$counts); }
foo(data[[1]])
data[[1]]$counts
data[[1]]
foo=function(dt){m=mean(dt$count); s=stdev(dt$count); }
foo(data[[1]])
foo=function(dt){m=mean(dt$count); s=sd(dt$count); }
foo(data[[1]])
foo=function(dt){m=mean(dt[[1]]$count); s=sd(dt[[1]]$count); dt[[1]][names(dt),counts]}
foo(data[[1]])
foo=function(dt){m=mean(dt[[1]]$count); s=sd(dt[[1]]$count); dt[[1]][names(dt),count]}
foo(data[[1]])
data[[1]]$count
foo(data[1])
foo=function(dt){m=mean(dt[[1]]$count); s=sd(dt[[1]]$count); (dt[[1]][names(dt),count]-m)/2}
foo(data[1])
foo=function(dt){m=mean(dt[[1]]$count); s=sd(dt[[1]]$count); (dt[[1]][names(dt),count]-m)/s}
foo(data[1])
data[names(data)[1]]
for(n in names(data)){foo(data[n])}
for(n in names(data)){print(foo(data[n]))}
sapply(names(data),function(n) foo(data[n])}
sapply(names(data),function(n) foo(data[n]))
wideScreen()
sapply(names(data),function(n) foo(data[n]))
foo=function(dt){m=mean(dt[[1]]$count); s=sd(dt[[1]]$count); (max(dt[[1]][,count])-m)/s}
sapply(names(data),function(n) foo(data[n]))
foo=function(dt){m=mean(dt[[1]]$count); s=sd(dt[[1]]$count); c('self'=dt[[1]][names(dt),count]-m)/s,'max'=(max(dt[[1]][,count])-m)/s,'count'=dt[[1]][names(dt),count])}
foo=function(dt){m=mean(dt[[1]]$count); s=sd(dt[[1]]$count); c('self'=(dt[[1]][names(dt),count]-m)/s,'max'=(max(dt[[1]][,count])-m)/s,'count'=dt[[1]][names(dt),count])}
sapply(names(data),function(n) foo(data[n]))
t(sapply(names(data),function(n) foo(data[n])))
signig(t(sapply(names(data),function(n) foo(data[n]))), digit=1)
signif(t(sapply(names(data),function(n) foo(data[n]))), digit=1)
signif(t(sapply(names(data),function(n) foo(data[n]))), digit=2)
foo
names(data)
dir(pattern='correl')
data=list()
for (f in dir (pattern="correl")){data[[f]]=fread(f,col.names=c("sample","count"),key='sample')}
names (data)
for (f in dir (pattern="correl")){data[[f]]=fread(f,col.names=c("sample","count"),key='sample')}
names(data)
names(data) = sub(".correl","",names(data))
names(data)
signif(t(sapply(names(data),function(n) foo(data[n]))), digit=2)
data[['B04']]
data[['T01']]
data[['ZH23']]
data[[c('ZH23','T01','B04')]]
data[c('ZH23','T01','B04')]
wideScreen()
data[c('ZH23','T01','B04')]
do.call(what='cbind',data[c('ZH23','T01','B04')])
data[1]
data[1]$count
data[[1]]$count
mean(data[[1]]$count)
median(data[[1]]$count)
sd(data[[1]]$count)
sd(data[[1]]$count-mean(data[[1]]$count))
sd(data[[1]]$count[-3])
sapply(data,function(e) sd(e$count))
sort(sapply(data,function(e) sd(e$count)))
sapply(data,function(e) mean(e$count))
sapply(data,function(e) c(mean(e$count),sd(e$count)))
sapply(data,function(e) c(mean(e$count),sd(e$count),mean(e$count)/sd(e$count)))
data[['B04']]
sapply(data,function(e) c(mean(e$count),sd(e$count),mean(e$count)/sd(e$count)))
sd(data[[1]]$count[-3])
sd(data[[1]]$count)
names(data)
sd(data[[10]]$count)
sd(data[[10]]$count[-10])
sd(data[[9]]$count)
sd(data[[9]]$count[-9])
sd(data[[8]]$count)
sd(data[[8]]$count[-8])
(data[[8]]$count)
sd(data[[8]]$count[-9])
sd(data[[8]]$count[-12])
getwd()
savehistory("/DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20191017_processing_cDNA_iPCR/Rhistory_iPCR_correlation_LP20191106.R")
