source('code/gene.fishing.package.R')
require(dplyr)
require(foreach)
require(doParallel)
registerDoParallel(30)


load('RData/util_data.RData')
PEER.raw.data              <- read.delim("original_data/expression/426sham0covK40addMeanPEERresidualsTranspose.txt", stringsAsFactors=FALSE)
rownames(PEER.raw.data)    <- PEER.raw.data$Gene
PEER.raw.data$Gene         <- NULL
PEER.raw.data.mt           <- t(as.matrix(PEER.raw.data))
PEER.expr.matrix           <- PEER.raw.data.mt[rownames(PEER.expr.matrix),]
c1                         <- read.table("raw.data/c1.csv", quote="\"", stringsAsFactors=FALSE)$V1


load('RData/GEUVADIS.peer.normalized.rpkm.RData')
overlapped.gene            <- intersect(colnames(PEER.expr.matrix),colnames(peer.normalized.rpkm.25)) %>% as.character

CHORI.jordan <- do.gene.fishing.v3(bait.gene = c1,
                                   expr.matrix = PEER.expr.matrix[,overlapped.gene],       
                                   alpha=5,
                                   fishing.round=1000,
                                   method.name='jordan')

CHORI.random.walk<- do.gene.fishing.v3(bait.gene = c1,
                                   expr.matrix = PEER.expr.matrix[,overlapped.gene],       
                                   alpha=5,
                                   fishing.round=1000,
                                   method.name='random.walk')

CHORI.berkeley.405 <- do.gene.fishing.v3(bait.gene = c1,
                                   expr.matrix = PEER.expr.matrix[,overlapped.gene],       
                                   alpha=5,
                                   fishing.round=1000,
                                   method.name='berkeley.405')



GEUVADIS.jordan <- do.gene.fishing.v3(bait.gene = c1,
                                      expr.matrix = peer.normalized.rpkm.25[,overlapped.gene],       
                                      alpha=5,
                                      fishing.round=1000,
                                      method.name='jordan')

GEUVADIS.random.walk <- do.gene.fishing.v3(bait.gene = c1,
                                                     expr.matrix = peer.normalized.rpkm.25[,overlapped.gene],       
                                                     alpha=5,
                                                     fishing.round=1000,
                                                     method.name='random.walk')

GEUVADIS.berkeley.405 <- do.gene.fishing.v3(bait.gene = c1,
                                            expr.matrix = peer.normalized.rpkm.25[,overlapped.gene],       
                                            alpha=5,
                                            fishing.round=1000,
                                            method.name='berkeley.405')



load('raw.data/GTex/normalized.GTex.with.ncRNA.median.rpkm.larger.than.0.1.RData')
liver.expr <- normalized.GTex.data[['Liver']]
GTex.liver.jordan <- do.gene.fishing.v3(bait.gene = c1,
                                        expr.matrix = liver.expr,       
                                        alpha=5,
                                        fishing.round=1000,
                                        method.name='jordan')

GTex.liver.random.walk <- do.gene.fishing.v3(bait.gene = c1,
                                                       expr.matrix = liver.expr,       
                                                       alpha=5,
                                                       fishing.round=1000,
                                                       method.name='random.walk')

GTex.liver.berkeley.405 <- do.gene.fishing.v3(bait.gene = c1,
                                              expr.matrix = liver.expr,       
                                              alpha=5,
                                              fishing.round=1000,
                                              method.name='berkeley.405')


save(file='RData/compare.pipeline.RData',
           list=c('CHORI.jordan',      'CHORI.random.walk',      'CHORI.berkeley.405',
                  'GEUVADIS.jordan',   'GEUVADIS.random.walk',   'GEUVADIS.berkeley.405',
                  'GTex.liver.jordan', 'GTex.liver.random.walk', 'GTex.liver.berkeley.405'
                 )
           )






