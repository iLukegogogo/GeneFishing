load('raw.data/GTex/normalized.GTex.RData')
source('code/gene.fishing.package.R')
require(dplyr)
require(foreach)
require(parallel)
require(doParallel)

registerDoParallel(20)
original.bait.gene  <- read.table("./raw.data/c1.csv", header=FALSE,stringsAsFactors=FALSE)[,1]
ep.list <- foreach(tissue = names(normalized.GTex.data)) %do% {
    expr.matrix <- normalized.GTex.data[[tissue]] 
    bait.gene   <- original.bait.gene[original.bait.gene %in% colnames(expr.matrix)]
    ep.vector   <- get.empirical.distance.ratio.distribution(expr.matrix = expr.matrix,bait.gene.number = length(bait.gene),alpha=5)
    p.value.vec <- foreach(j = 1:1000,.combine='c') %dopar%{
        fishing.once(bait.gene = bait.gene,expr.matrix = expr.matrix,alpha=5,ep.vecotr = ep.vector)  
    }
    tissue <- gsub(x = tissue,pattern = '\\s+',replacement = "")
    file.name <- sprintf('Figures/%s.pdf',tissue)
    pdf(file.name)
    hist(ep.vector,  xlab='ratio',                  main=sprintf('histogram of ratio (%s)'  ,tissue))
    hist(p.value.vec,xlab='p-value (21 bait genes)',main=sprintf('histogram of p-value (%s)',tissue))
    dev.off()
    ep.vector
}
names(ep.list) <- names(normalized.GTex.data)


