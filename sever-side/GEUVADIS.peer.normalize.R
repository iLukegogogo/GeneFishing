require(dplyr)
require(peer)

rpkm.data              <- read.delim("raw.data/GEUVADIS/GD660.GeneQuantRPKM.txt", stringsAsFactors=FALSE)
rpkm.data$TargetID     <- NULL
rpkm.data$Chr          <- NULL
rpkm.data$Coord        <- NULL
rpkm.data$Gene_Symbol  <- gsub(x=rpkm.data$Gene_Symbol, pattern="\\.\\d+", replacement="",perl=TRUE)
rownames(rpkm.data)    <- rpkm.data$Gene_Symbol
rpkm.data$Gene_Symbol  <- NULL
rpkm.matrix            <- as.matrix(rpkm.data) %>% t
tmp                    <- apply(rpkm.matrix,MARGIN=2,function(x) median(x))
rpkm.matrix            <- rpkm.matrix[,tmp>=0.1]


# peer.obj <-  PEER()
# PEER_setAdd_mean(peer.obj,TRUE)
# PEER_setNk(peer.obj,10)
# PEER_setPhenoMean(peer.obj,log10(rpkm.matrix+0.01))
# PEER_update(peer.obj)
# peer.normalized.rpkm.10            <- PEER_getResiduals(peer.obj)
# rownames(peer.normalized.rpkm.10)  <- rownames(rpkm.matrix)
# colnames(peer.normalized.rpkm.10)  <- colnames(rpkm.matrix)
# 
# 
# 
# 
# 
# 
# peer.obj <-  PEER()
# PEER_setAdd_mean(peer.obj,TRUE)
# PEER_setNk(peer.obj,15)
# PEER_setPhenoMean(peer.obj,log10(rpkm.matrix+1))
# PEER_update(peer.obj)
# peer.normalized.rpkm.15            <- PEER_getResiduals(peer.obj)
# rownames(peer.normalized.rpkm.15)  <- rownames(rpkm.matrix)
# colnames(peer.normalized.rpkm.15)  <- colnames(rpkm.matrix)
# 
# 
# 
# 
# peer.obj <-  PEER()
# PEER_setAdd_mean(peer.obj,TRUE)
# PEER_setNk(peer.obj,20)
# PEER_setPhenoMean(peer.obj,log(rpkm.matrix+1))
# PEER_update(peer.obj)
# peer.normalized.rpkm.20            <- PEER_getResiduals(peer.obj)
# rownames(peer.normalized.rpkm.20)  <- rownames(rpkm.matrix)
# colnames(peer.normalized.rpkm.20)  <- colnames(rpkm.matrix)




peer.obj <-  PEER()
PEER_setAdd_mean(peer.obj,TRUE)
PEER_setNk(peer.obj,25)
PEER_setPhenoMean(peer.obj,log10(rpkm.matrix+0.01))
PEER_update(peer.obj)
peer.normalized.rpkm.25            <- PEER_getResiduals(peer.obj)
rownames(peer.normalized.rpkm.25)  <- rownames(rpkm.matrix)
colnames(peer.normalized.rpkm.25)  <- colnames(rpkm.matrix)




# peer.obj <-  PEER()
# PEER_setAdd_mean(peer.obj,TRUE)
# PEER_setNk(peer.obj,30)
# PEER_setPhenoMean(peer.obj,log(rpkm.matrix+1))
# PEER_update(peer.obj)
# peer.normalized.rpkm.30            <- PEER_getResiduals(peer.obj)
# rownames(peer.normalized.rpkm.30)  <- rownames(rpkm.matrix)
# colnames(peer.normalized.rpkm.30)  <- colnames(rpkm.matrix)
# 
# 
# peer.obj <-  PEER()
# PEER_setAdd_mean(peer.obj,TRUE)
# PEER_setNk(peer.obj,50)
# PEER_setPhenoMean(peer.obj,log(rpkm.matrix+1))
# PEER_update(peer.obj)
# peer.normalized.rpkm.50            <- PEER_getResiduals(peer.obj)
# rownames(peer.normalized.rpkm.50)  <- rownames(rpkm.matrix)
# colnames(peer.normalized.rpkm.50)  <- colnames(rpkm.matrix)


save(file='RData/GEUVADIS.peer.normalized.rpkm.RData',list=c('peer.normalized.rpkm.25'))



# save(file='RData/GEUVADIS.peer.normalized.rpkm.RData',
#      list=c('peer.normalized.rpkm.25','peer.normalized.rpkm.20','peer.normalized.rpkm.15',
#             'peer.normalized.rpkm.30','peer.normalized.rpkm.50','peer.normalized.rpkm.10'
#             )
#     )



pdf('pca.pdf')
# 
# pca.result.1 <- prcomp(peer.normalized.rpkm.10,center = T,retx = TRUE,scale. = FALSE)
# pca.result.2 <- prcomp(peer.normalized.rpkm.10,center = T,retx = TRUE,scale. = TRUE)
# pairs(pca.result.1$x[,1:5])
# pairs(pca.result.2$x[,1:5])
# 
# 
# 
# 
# 
# 
# pca.result.1 <- prcomp(peer.normalized.rpkm.15,center = T,retx = TRUE,scale. = FALSE)
# pca.result.2 <- prcomp(peer.normalized.rpkm.15,center = T,retx = TRUE,scale. = TRUE)
# pairs(pca.result.1$x[,1:5])
# pairs(pca.result.2$x[,1:5])
# 
# 
# 
# 
# pca.result.1 <- prcomp(peer.normalized.rpkm.20,center = T,retx = TRUE,scale. = FALSE)
# pca.result.2 <- prcomp(peer.normalized.rpkm.20,center = T,retx = TRUE,scale. = TRUE)
# pairs(pca.result.1$x[,1:5])
# pairs(pca.result.2$x[,1:5])
# 
# 
# 
# 
pca.result.1 <- prcomp(peer.normalized.rpkm.25,center = T,retx = TRUE,scale. = FALSE)
pca.result.2 <- prcomp(peer.normalized.rpkm.25,center = T,retx = TRUE,scale. = TRUE)
pairs(pca.result.1$x[,1:5])
pairs(pca.result.2$x[,1:5])
# 
# 
# 
# pca.result.1 <- prcomp(peer.normalized.rpkm.30,center = T,retx = TRUE,scale. = FALSE)
# pca.result.2 <- prcomp(peer.normalized.rpkm.30,center = T,retx = TRUE,scale. = TRUE)
# pairs(pca.result.1$x[,1:5])
# pairs(pca.result.2$x[,1:5])
# 
# 
dev.off()

