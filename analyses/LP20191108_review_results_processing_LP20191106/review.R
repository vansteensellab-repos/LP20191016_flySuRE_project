library(data.table)
fnames <- grep(value=TRUE, pattern="LP20191106.*11_sorted",dir(path='/DATA/usr/ludo/projects/LP190425_flySuRE/data/intermediate/',recursive=T, include.dirs=T,pattern=".*gz",full.names=T))
fname.sample <- sub("^.*LP20191106_(.*)_pipelineOutp.*$","\\1",fnames)
samples <- unique(fname.sample)
data <- data.frame(libComplFlat=0, libCompl=0,cdnaComplFlat=0,cdnaCompl=0,iPCRraw=0, cDNAraw=0, row.names=s)[-1,]
for (s in samples) {
  data[s,] <- data.frame(libComplFlat=0, libCompl=0,cdnaComplFlat=0,cdnaCompl=0,iPCRraw=0, cDNAraw=0)
  fs <- fnames[fname.sample %in% s]
  for (f in fs) {
    dt <- fread(f, select=c(6,14), sep="\t")
    data[s,c('libComplFlat', 'cdnaComplFlat', 'libCompl', 'cdnaCompl')] <- 
      data[s,c('libComplFlat', 'cdnaComplFlat', 'libCompl', 'cdnaCompl')] + c(colSums(dt), colSums(dt>0))
  }
}
data[,'iPCRraw'] <- c(71.1, 52.3, 66.9, 74.6, 70.3, 87.3, 75.1, 70.5, 79.3, 78.1)
data[,'cDNAraw'] <- c(25.9, 20.6, 26.9, 33.1, 34.6, 34.4, 33.5, 33.0, 5.1, 21.8)

