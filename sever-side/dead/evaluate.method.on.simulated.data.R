source('code/gene.fishing.package.R')
require(dplyr)
require(foreach)
require(doParallel)
registerDoParallel(30)


l                          <- readRDS('RData/simData_RandomCov_1Cluster.RDS')
gene.expr                  <- l$simData
rownames(gene.expr)        <- paste('gene',1:nrow(gene.expr),sep=".")
colnames(gene.expr)        <- paste('sample',1:ncol(gene.expr),sep=".")
rownames(gene.expr)[1:100] <- paste('bait',1:100,sep=".")
gene.expr                  <- t(gene.expr)

rs <- foreach(i=1:100) %do% {
    c1  <- sample(colnames(gene.expr),21)

    jordan <- do.gene.fishing.v3(bait.gene = c1,
                                 expr.matrix = gene.expr,       
                                 alpha=5,
                                 fishing.round=1000,
                                 method.name='jordan')

random.walk<- do.gene.fishing.v3(bait.gene = c1,
                                 expr.matrix =gene.expr,       
                                 alpha=5,
                                 fishing.round=1000,
                                 method.name='random.walk')

berkeley.405 <- do.gene.fishing.v3(bait.gene = c1,
                                   expr.matrix = gene.expr,       
                                   alpha=5,
                                   fishing.round=1000,
                                   method.name='berkeley.405')
jordan       <- arrange(jordan,     desc(fish.freq))
random.walk  <- arrange(random.walk,desc(fish.freq))
berkeley.405 <- arrange(berkeley.405,desc(fish.freq))

jordan       <- jordan[jordan$fish.id != 'NOTHING',]
random.walk  <- random.walk[random.walk$fish.id != 'NOTHING',]
berkeley.405 <- berkeley.405[berkeley.405$fish.id != 'NOTHING',]



list(jordan=jordan,random.walk=random.walk,berkeley.405=berkeley.405)
}
save(file='RData/evaluate.method.on.simulated.data.RData',list=c('rs'))

