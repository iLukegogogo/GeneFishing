require(data.table)
require(plyr)
require(dplyr)
require(foreach)
require(peer)

data                            <-  fread(input = 'raw.data/CHORI/39128by924counts.txt',header=TRUE, stringsAsFactors = FALSE)  %>% as.data.frame 
flag                            <-  grepl(x=colnames(data),pattern='sham')
CHORI.read.cnt.df               <-  cbind(data[,1],data[,flag])
colnames(CHORI.read.cnt.df)     <-  gsub(colnames(CHORI.read.cnt.df),pattern='sham',replacement="")
colnames(CHORI.read.cnt.df)[1]  <-  'gene.id'
rownames(CHORI.read.cnt.df )    <-  CHORI.read.cnt.df$gene.id
CHORI.read.cnt.df$gene.id       <-  NULL
CHORI.read.cnt.matrix           <-  as.matrix(CHORI.read.cnt.df)
cell.line.to.cs.id              <-  read.table('raw.data/CHORI/N426CSIDs.txt',header=T)
colnames(cell.line.to.cs.id)    <-  c('cell.line','cs.id')
rownames(cell.line.to.cs.id)    <-  cell.line.to.cs.id$cell.line
CHORI.read.cnt.matrix           <-  CHORI.read.cnt.matrix[,colnames(CHORI.read.cnt.matrix) %in% cell.line.to.cs.id$cell.line]
colnames(CHORI.read.cnt.matrix) <-  cell.line.to.cs.id[colnames(CHORI.read.cnt.matrix),'cs.id']


meta.df           <- read.table(file = 'raw.data/CHORI/N426shamStatinTotalAlignedCounts.txt',header = T,stringsAsFactors = F)
meta.df           <- meta.df[meta.df$Treatment == 'sham',]
rownames(meta.df) <- meta.df$CSID %>% as.character


lib.depth  <- meta.df[ as.character(colnames(CHORI.read.cnt.matrix)),'totalAligned' ]
rpm.matrix <- foreach(j=1:ncol(CHORI.read.cnt.matrix),.combine='cbind') %do% {
    CHORI.read.cnt.matrix[,j]/(lib.depth[j]/1000000)
}
rownames(rpm.matrix)   <- rownames(CHORI.read.cnt.matrix)
colnames(rpm.matrix)   <- colnames(CHORI.read.cnt.matrix)



ensemble.87.tr.length           <- read.csv(file = 'raw.data/CHORI/ensemble.87.tr.length.txt',header = F,stringsAsFactors = F)
colnames(ensemble.87.tr.length) <- c('gene.id','tr.id','tr.length')
gene.mean.tr.length             <- ddply(ensemble.87.tr.length,.(gene.id),function(x) mean(x$tr.length))
colnames(gene.mean.tr.length)   <- c('gene.id','mean.tr.length')
rownames(gene.mean.tr.length)   <- gene.mean.tr.length$gene.id %>% as.character 


flag                   <- rownames(rpm.matrix) %in% rownames(gene.mean.tr.length)
rpm.matrix             <- rpm.matrix[flag,]
mean.tr.length         <- gene.mean.tr.length[as.character(rownames(rpm.matrix)),'mean.tr.length']
rpkm.matrix            <- foreach(i=1:nrow(rpm.matrix),.combine='rbind') %do% {
    rpm.matrix[i,]/(mean.tr.length[i]/1000)
}
rownames(rpkm.matrix) <- rownames(rpm.matrix)
colnames(rpkm.matrix) <- colnames(rpm.matrix)


colnames(rpkm.matrix) <- meta.df[ as.character(colnames(rpkm.matrix)),'CellLine' ]
rpkm.matrix           <- t(rpkm.matrix)
tmp                   <- apply(rpkm.matrix,MARGIN=2,function(x) median(x))
rpkm.matrix           <- rpkm.matrix[,tmp>=0.1]

load('RData/util_data.RData')
rpkm.matrix <- rpkm.matrix[rownames(PEER.expr.matrix),]

peer.obj <-  PEER()
PEER_setAdd_mean(peer.obj,TRUE)
PEER_setNk(peer.obj,25)
PEER_setPhenoMean(peer.obj,log10(rpkm.matrix+0.01))
PEER_update(peer.obj)
peer.normalized.rpkm.25            <- PEER_getResiduals(peer.obj)
rownames(peer.normalized.rpkm.25)  <- rownames(rpkm.matrix)
colnames(peer.normalized.rpkm.25)  <- colnames(rpkm.matrix)
save(file='RData/CHORI.peer.normalized.rpkm.RData',list=c('peer.normalized.rpkm.25'))


pdf('pca.pdf')
pca.result.1 <- prcomp(peer.normalized.rpkm.25,center = T,retx = TRUE,scale. = FALSE)
pca.result.2 <- prcomp(peer.normalized.rpkm.25,center = T,retx = TRUE,scale. = TRUE)
pairs(pca.result.1$x[,1:5])
pairs(pca.result.2$x[,1:5])
dev.off()


