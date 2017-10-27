source('code/gene.fishing.package.R')
require(dplyr)
require(foreach)
require(doParallel)
registerDoParallel(40)

c1                         <- read.table("raw.data/c1.csv", quote="\"", stringsAsFactors=FALSE)$V1


load('raw.data/GTex/normalized.GTex.with.ncRNA.median.rpkm.larger.than.0.1.RData')
GTex.liver.expr <- normalized.GTex.data[['Liver']]
non.bait.genes  <- setdiff(colnames(GTex.liver.expr),c1)
GTex.liver.p.hat.matrix <- foreach(i = 1:20,.combine='rbind') %do% {
    rd.matrix              <- GTex.liver.expr[,non.bait.genes]
    rd.matrix              <- apply(rd.matrix,MARGIN=2,sample)
    permutated.expr.matrix <- cbind(GTex.liver.expr[,c1],rd.matrix)
    GTex.liver.random.walk.permutation  <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = permutated.expr.matrix,alpha=5,fishing.round=5000,method.name='random.walk')
    GTex.liver.berkeley.405.permutation <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = permutated.expr.matrix,alpha=5,fishing.round=5000,method.name='berkeley.405')
    GTex.liver.random.walk.permutation  <- GTex.liver.random.walk.permutation [GTex.liver.random.walk.permutation$fish.id  != 'NOTHING',]
    GTex.liver.berkeley.405.permutation <- GTex.liver.berkeley.405.permutation[GTex.liver.berkeley.405.permutation$fish.id != 'NOTHING',]
    p.hat.random.walk                   <- sum(GTex.liver.random.walk.permutation$fish.freq)/ncol(rd.matrix) # average of capture freq is the esitmators of p
    p.hat.berkeley.405                  <- sum(GTex.liver.berkeley.405.permutation$fish.freq)/ncol(rd.matrix)
    c(p.hat.random.walk,p.hat.berkeley.405)
    
}






#load('RData/CHORI.peer.normalized.rpkm.RData')
#PEER.expr.matrix <- peer.normalized.rpkm.25

load('RData/util_data.RData')
PEER.raw.data              <- read.delim("original_data/expression/426sham0covK40addMeanPEERresidualsTranspose.txt", stringsAsFactors=FALSE)
rownames(PEER.raw.data)    <- PEER.raw.data$Gene
PEER.raw.data$Gene         <- NULL
PEER.raw.data.mt           <- t(as.matrix(PEER.raw.data))
PEER.expr.matrix           <- PEER.raw.data.mt[rownames(PEER.expr.matrix),]

non.bait.genes              <- setdiff(colnames(PEER.expr.matrix),c1)
CHORI.p.hat.matrix <- foreach(i = 1:20,.combine='rbind') %do% {
    rd.matrix              <- PEER.expr.matrix[,non.bait.genes]
    rd.matrix              <- apply(rd.matrix,MARGIN=2,sample)
    permutated.expr.matrix <- cbind(PEER.expr.matrix[,c1],rd.matrix)
    CHORI.random.walk.permutation  <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = permutated.expr.matrix,alpha=5,fishing.round=5000,method.name='random.walk')
    CHORI.berkeley.405.permutation <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = permutated.expr.matrix,alpha=5,fishing.round=5000,method.name='berkeley.405')
    CHORI.random.walk.permutation  <- CHORI.random.walk.permutation [CHORI.random.walk.permutation$fish.id != 'NOTHING',]
    CHORI.berkeley.405.permutation <- CHORI.berkeley.405.permutation[CHORI.berkeley.405.permutation$fish.id != 'NOTHING',]
    p.hat.random.walk              <- sum(CHORI.random.walk.permutation$fish.freq)/ncol(rd.matrix) # average of capture freq is the esitmators of p
    p.hat.berkeley.405             <- sum(CHORI.berkeley.405.permutation$fish.freq)/ncol(rd.matrix)
    c(p.hat.random.walk,p.hat.berkeley.405)
    
}
#save(file = 'RData/permutation.RData',list = c('CHORI.p.hat.matrix'))

save(file = 'RData/permutation.RData',list = c('GTex.liver.p.hat.matrix','CHORI.p.hat.matrix'))











# fish.list <- foreach( tissue  = c('Liver') )%do%{
#     tissue.expr.matrix <- 
#     if(nrow(tissue.expr.matrix) <100){
#         fish <- c('NOTHING')  
#         fish
#     }else{
#         fish.freq.df
#     }
# }
# names(fish.list) <- c('Liver')
# save(file='RData/fishing.in.GTex.with.ncRNA.median.rpkm.larger.than.0.1.permutation.RData',list=c('fish.list'))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# load('RData/util_data.RData')
# PEER.raw.data              <- read.delim("original_data/expression/426sham0covK40addMeanPEERresidualsTranspose.txt", stringsAsFactors=FALSE)
# rownames(PEER.raw.data)    <- PEER.raw.data$Gene
# PEER.raw.data$Gene         <- NULL
# PEER.raw.data.mt           <- t(as.matrix(PEER.raw.data))
# PEER.expr.matrix           <- PEER.raw.data.mt[rownames(PEER.expr.matrix),]
# 
# 
# c1                         <- read.table("raw.data/c1.csv", quote="\"", stringsAsFactors=FALSE)$V1
# CHORI.LCL.fish.freq.df.permutation     <- do.gene.fishing.with.permutation(bait.gene = c1,expr.matrix = PEER.expr.matrix,alpha=5,       fishing.round=1000)
# 








