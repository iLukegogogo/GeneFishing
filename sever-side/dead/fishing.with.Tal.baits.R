source('code/gene.fishing.package.R')
require(dplyr)
require(foreach)
require(doParallel)
registerDoParallel(30)
load('RData/util_data.RData')
load('RData/GEUVADIS.peer.normalized.rpkm.RData')
load('RData/tal.sample.baits.RData')

PEER.raw.data              <- read.delim("original_data/expression/426sham0covK40addMeanPEERresidualsTranspose.txt", stringsAsFactors=FALSE)
rownames(PEER.raw.data)    <- PEER.raw.data$Gene
PEER.raw.data$Gene         <- NULL
PEER.raw.data.mt           <- t(as.matrix(PEER.raw.data))
PEER.expr.matrix           <- PEER.raw.data.mt[rownames(PEER.expr.matrix),]

CHORI.LCL.sample.baits.result.list <- foreach(tal.bait = tal.sample.baits.list) %do% {
    rs <- do.gene.fishing.v3(bait.gene = tal.bait,expr.matrix = PEER.expr.matrix,alpha=5,fishing.round=1000)
    rs <- rs[rs$fish.id != 'NOTHING',]
    rs
}
names(CHORI.LCL.sample.baits.result.list) <- names(tal.sample.baits.list)


GEUVADIS.LCL.sample.baits.result.list <- foreach(tal.bait = tal.sample.baits.list) %do% {
  rs <- do.gene.fishing.v3(bait.gene = tal.bait,expr.matrix = peer.normalized.rpkm.25,alpha=5,fishing.round=1000)
  rs <- rs[rs$fish.id != 'NOTHING',]
  rs
  r
}
names(GEUVADIS.LCL.sample.baits.result.list) <- names(tal.sample.baits.list)




load('raw.data/GTex/normalized.GTex.with.ncRNA.median.rpkm.larger.than.0.1.RData')
GTex.liver.sample.baits.result.list <- foreach(tal.bait = tal.sample.baits.list) %do% {
  rs <- do.gene.fishing.v3(bait.gene = tal.bait,expr.matrix = normalized.GTex.data[['Liver']],alpha=5,fishing.round=1000)
  rs <- rs[rs$fish.id != 'NOTHING',]
  rs
  r
}
names(GTex.liver.sample.baits.result.list) <- names(tal.sample.baits.list)



save(
     file = 'RData/fishing.with.Tal.baits.RData',
     list = c('CHORI.LCL.sample.baits.result.list','GEUVADIS.LCL.sample.baits.result.list','GTex.liver.sample.baits.result.list')
     )