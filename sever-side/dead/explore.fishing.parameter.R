source('code/gene.fishing.package.R')
load('raw.data/GTex/normalized.GTex.RData')
require(dplyr)
require(foreach)
require(doParallel)

registerDoParallel(30)

tissue             <- 'Liver'
c1                 <- read.table("./raw.data/c1.csv", header=FALSE,stringsAsFactors=FALSE)[,1]
tissue.expr.matrix <- normalized.GTex.data[[tissue]]
gene.freq.data     <- foreach(alpha = 5:20) %do%{
    rs <- foreach(i = 1:20) %do% {
        fish.freq.df               <- do.gene.fishing(bait.gene = c1,expr.matrix = tissue.expr.matrix,alpha=alpha,fishing.time=1000)
        fish.freq.df               <- fish.freq.df[fish.freq.df$fish.id != 'NOTHING',]
        fish.freq.df$fish.id       <- as.character(fish.freq.df$fish.id)
        data                       <- rep(0,ncol(tissue.expr.matrix))
        names(data)                <- colnames(tissue.expr.matrix)
        data[fish.freq.df$fish.id] <- fish.freq.df$fish.freq * 1000
        data
    }
    rs
}
save(file = 'RData/explore.fishing.paramater.RData',list = c('gene.freq.data'))








do.mixture.EM <- function(data,n){
    p1 <- 0.1
    p2 <- 0.3
    pi <- 0.4
    N  <- length(data)
    for(i in 1:20){
        A <- foreach(x=data,.combine='c') %do% {
            f <- exp(x*log(p2/p1) + (n-x)*log((1-p2)/(1-p1)))
            pi/(pi + (1-pi)*f)
        }    
      
        B <- foreach(x=data,.combine='c') %do% {
            f <- exp(x*log(p1/p2) + (n-x)*log((1-p1)/(1-p2))) 
            (1-pi)/((1-pi) + pi*f)
        } 
        pi <- sum(A)/N
        p1 <- crossprod(data,A)/sum(A * n)
        p2 <- crossprod(data,B)/sum(B * n)
        print(c(pi,p1,p2))
    }
    B
}

