source('code/gene.fishing.package.R')
require(dplyr)
require(foreach)
require(doParallel)
registerDoParallel(40)

c1        <- read.table("./raw.data/c1.csv", header=FALSE,stringsAsFactors=FALSE)[,1]

load('raw.data/GTex/normalized.GTex.with.ncRNA.median.rpkm.larger.than.0.1.RData')
fish.list <- foreach( tissue  = names(normalized.GTex.data) )%do%{
    tissue.expr.matrix <- normalized.GTex.data[[tissue]]
    if(nrow(tissue.expr.matrix) <100){
        fish <- c('NOTHING')  
        fish
    }else{
        fish.freq.df <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = tissue.expr.matrix,alpha=5,fishing.round=1000)
        fish.freq.df
    }
}
names(fish.list) <- names(normalized.GTex.data)
save(file='RData/fishing.in.GTex.with.ncRNA.median.rpkm.larger.than.0.1.RData',list=c('fish.list'))



# load('raw.data/GTex/normalized.GTex.with.ncRNA.RData')
# fish.list <- foreach( tissue  = names(normalized.GTex.data) )%do%{
#   tissue.expr.matrix <- normalized.GTex.data[[tissue]]
#   if(nrow(tissue.expr.matrix) <100){
#     fish <- c('NOTHING')  
#     fish
#   }else{
#     fish.freq.df <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = tissue.expr.matrix,alpha=5,fishing.round=1000)
#     fish.freq.df
#   }
# }
# names(fish.list) <- names(normalized.GTex.data)
# save(file='RData/fishing.in.GTex.with.ncRNA.at.least.10.sample.rpkm.larger.than.0.1.RData',list=c('fish.list'))
# 
# 
# 
# 
# load('raw.data/GTex/normalized.GTex.RData')
# fish.list <- foreach( tissue  = names(normalized.GTex.data) )%do%{
#   tissue.expr.matrix <- normalized.GTex.data[[tissue]]
#   if(nrow(tissue.expr.matrix) <100){
#     fish <- c('NOTHING')  
#     fish
#   }else{
#     fish.freq.df <- do.gene.fishing.v3(bait.gene = c1,expr.matrix = tissue.expr.matrix,alpha=5,fishing.round=1000)
#     fish.freq.df
#   }
# }
# names(fish.list) <- names(normalized.GTex.data)
# save(file='RData/fishing.in.GTex.without.ncRNA.at.least.10.sample.rpkm.larger.than.0.1.RData',list=c('fish.list'))







