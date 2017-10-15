require(foreach)
require(plyr)
require(dplyr)
require(ggplot2)
require(doParallel)
require(LICORS)
registerDoParallel(40)
source('code/gene.fishing.package.R')

l                          <- readRDS('RData/simData_realCov_1cluster_10000Genes.RDS')
gene.expr                  <- l$simData
rownames(gene.expr)        <- paste('gene',1:nrow(gene.expr),sep=".")
colnames(gene.expr)        <- paste('sample',1:ncol(gene.expr),sep=".")
rownames(gene.expr)[1:l$size.cluster1] <- paste('bait',1:l$size.cluster1,sep=".")
gene.expr                  <- t(gene.expr)

bait.gene.list <-  foreach(i=1:50) %do% {
    sample(colnames(gene.expr)[1:l$size.cluster1],21)  
}

gene.fishing.rs <- foreach(c1 = bait.gene.list) %do% {
    df <- do.gene.fishing.v3(bait.gene   = c1,
                             expr.matrix = gene.expr,
                             alpha=5,
                             fishing.round=1000,
                             method.name='berkeley.405')  
    df <- df[df$fish.id !='NOTHING',]
    df <- arrange(df,desc(fish.freq))
    df
}
save(file = 'RData/fishing.in.sim.data.RData',list=c('gene.expr','bait.gene.list','gene.fishing.rs'))



