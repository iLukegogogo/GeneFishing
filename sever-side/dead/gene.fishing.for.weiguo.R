source('code/gene.fishing.package.R')
PEER.raw.data              <- read.delim("./original_data/expression/426sham0covK40addMeanPEERresidualsTranspose.txt", stringsAsFactors=FALSE)
rownames(PEER.raw.data)    <- PEER.raw.data$Gene
PEER.raw.data$Gene         <- NULL
PEER.raw.data.mt           <- t(as.matrix(PEER.raw.data))
PEER.expr.matrix           <- PEER.raw.data.mt


weiguo.gene        <- read.delim("./raw.data/weiguo.gene.txt", header=FALSE, dec=",", stringsAsFactors=FALSE)$V2
CHORI.fish.freq.df <- do.gene.fishing(bait.gene = weiguo.gene,expr.matrix = PEER.expr.matrix,alpha=5,fishing.time=1000)


load('RData/GEUVADIS.peer.normalized.rpkm.RData')
GEUVADIS.fish.freq.df <- do.gene.fishing(bait.gene = weiguo.gene,expr.matrix=peer.normalized.rpkm.20,alpha=5,fishing.time=1000)


save(file='RData/gene.fishing.for.weiguo.RData',list=c('CHORI.fish.freq.df','GEUVADIS.fish.freq.df'))

