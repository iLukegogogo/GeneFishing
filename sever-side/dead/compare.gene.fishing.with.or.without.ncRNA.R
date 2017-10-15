source('code/gene.fishing.package.R')
require(dplyr)
require(foreach)
require(doParallel)

registerDoParallel(30)
c1        <- read.table("./raw.data/c1.csv", header=FALSE,stringsAsFactors=FALSE)[,1]
tissue    <- 'Liver'

load('raw.data/GTex/normalized.GTex.with.ncRNA.RData')
tissue.expr.matrix            <- normalized.GTex.data[[tissue]]
with.ncRNA.fish.freq.df       <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = tissue.expr.matrix,alpha=5,fishing.round=1000)


load('raw.data/GTex/normalized.GTex.RData')
tissue.expr.matrix            <- normalized.GTex.data[[tissue]]
without.ncRNA.fish.freq.df    <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = tissue.expr.matrix,alpha=5,fishing.round=1000)
