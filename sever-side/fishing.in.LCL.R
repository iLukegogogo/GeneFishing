source('code/gene.fishing.package.R')
load('RData/util_data.RData')

require(dplyr)
require(foreach)
require(doParallel)
registerDoParallel(30)

c1                         <- read.table("raw.data/c1.csv", quote="\"", stringsAsFactors=FALSE)$V1

PEER.raw.data              <- read.delim("original_data/expression/426sham0covK40addMeanPEERresidualsTranspose.txt", stringsAsFactors=FALSE)
rownames(PEER.raw.data)    <- PEER.raw.data$Gene
PEER.raw.data$Gene         <- NULL
PEER.raw.data.mt           <- t(as.matrix(PEER.raw.data))
PEER.expr.matrix           <- PEER.raw.data.mt[rownames(PEER.expr.matrix),]
CHORI.beth.expr.matrix     <- PEER.expr.matrix


load('RData/CHORI.peer.normalized.rpkm.RData')
CHORI.ke.expr.matrix <- peer.normalized.rpkm.25


load('RData/GEUVADIS.peer.normalized.rpkm.RData')
GEUVADIS.expr.matrix <- peer.normalized.rpkm.25


overlapped.gene                         <- intersect(colnames(CHORI.beth.expr.matrix),colnames(GEUVADIS.expr.matrix)) %>% as.character
CHORI.LCL.intersected.fish.freq.df      <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = CHORI.beth.expr.matrix[,overlapped.gene],alpha=5,fishing.round=1000)
GEUVADIS.LCL.intersected.fish.freq.df   <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = GEUVADIS.expr.matrix  [,overlapped.gene],alpha=5,fishing.round=1000)
CHORI.LCL.fish.freq.df                  <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = CHORI.beth.expr.matrix,                  alpha=5,fishing.round=1000)
GEUVADIS.LCL.fish.freq.df               <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = GEUVADIS.expr.matrix  ,                  alpha=5,fishing.round=1000)


save(file = 'RData/fishing.in.LCL.CHORI.beth.vs.GEUVADIS.RData',list = c('CHORI.LCL.fish.freq.df','GEUVADIS.LCL.fish.freq.df','CHORI.LCL.intersected.fish.freq.df','GEUVADIS.LCL.intersected.fish.freq.df'))



overlapped.gene                         <- intersect(colnames(CHORI.ke.expr.matrix),colnames(GEUVADIS.expr.matrix)) %>% as.character
CHORI.LCL.intersected.fish.freq.df      <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = CHORI.ke.expr.matrix[,overlapped.gene],alpha=5,fishing.round=1000)
GEUVADIS.LCL.intersected.fish.freq.df   <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = GEUVADIS.expr.matrix[,overlapped.gene],alpha=5,fishing.round=1000)
CHORI.LCL.fish.freq.df                  <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = CHORI.ke.expr.matrix,                  alpha=5,fishing.round=1000)
GEUVADIS.LCL.fish.freq.df               <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = GEUVADIS.expr.matrix,                  alpha=5,fishing.round=1000)


save(file = 'RData/fishing.in.LCL.CHORI.ke.vs.GEUVADIS.RData',list = c('CHORI.LCL.fish.freq.df','GEUVADIS.LCL.fish.freq.df','CHORI.LCL.intersected.fish.freq.df','GEUVADIS.LCL.intersected.fish.freq.df'))








