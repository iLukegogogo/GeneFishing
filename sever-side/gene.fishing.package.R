require(foreach)
require(doParallel)
require(plyr)
require(dplyr)
require(LICORS)



get.spectral.clustering.coordinates <- function(A,force.no.of.cluster=0){ 
    diag(A)  <-   0
    A        <-  abs(A)
    d        <-  apply(A,1,sum)
    I        <-  diag(1,nrow = nrow(A)) 
    L        <-   I -  diag(1/sqrt(d)) %*% A %*% diag(1/sqrt(d))
    

    #L is positive semi-definite and have n non-negative real-valued eigenvalues
    tmp             <-  eigen(L)
    eigen.values    <-  tmp$values[length(tmp$values):1]
    eigen.vectors   <-  tmp$vectors[,length(tmp$values):1]
    
    q1              <-  quantile(eigen.values)[2]
    q3              <-  quantile(eigen.values)[4]
    low             <-  q1 - 3*(q3 - q1)
    up              <-  q1 + 3*(q3 - q1)
    k               <-  sum(eigen.values < low) # pick out eigen-values that are relatively small
    
    if(force.no.of.cluster > 1){
        k <- force.no.of.cluster
    }
    if(k == 1){
        list(eigen.values=eigen.values,no.of.cluster=1)
    }else{
       coordinate.matrix            <- eigen.vectors[,2:(k+1)] # pick out the eigen-vectors associated with the K smallest non-zero eigen-values
       rownames(coordinate.matrix)  <- rownames(A)
       colnames(coordinate.matrix)  <- paste('dim',1:(k),sep="-")
       list(coordinates=coordinate.matrix,eigen.values=eigen.values,no.of.cluster=k)
    }
}



probe.fishability <- function(bait.gene, expr.matrix,alpha=5){
    bait.gene <- bait.gene[bait.gene %in% colnames(expr.matrix)]
    all.gene  <- colnames(expr.matrix)      
    pool.gene <- setdiff(all.gene,bait.gene)
    bait.cluster.purity.matrix <- foreach(k =c(2:10),.combine='rbind') %do%{
        bait.cluster.purity      <- foreach(i= 1:1000,.combine='c') %dopar% {
            bait.gene.expr.matrix          <- expr.matrix[,bait.gene]
            rd.gene                        <- sample(pool.gene,size=length(bait.gene) * alpha)
            bait.and.pool.gene.expr.matrix <- cbind(expr.matrix[,bait.gene],expr.matrix[,rd.gene])
            bait.and.pool.gene.cor.matrix  <- cor(bait.and.pool.gene.expr.matrix,method='spearman')
            rs                             <- get.spectral.clustering.coordinates(bait.and.pool.gene.cor.matrix,k)
            if(rs[['no.of.cluster']] == 1){
                return(1/(1+alpha))
            }
        
            coordinates                    <- rs[['coordinates']]
            kmeans.obj                     <- kmeanspp(data=coordinates,k=ncol(coordinates))
            cluster.vec                    <- kmeans.obj$cluster
            names(cluster.vec)             <- rownames(coordinates)
            cluster.freq.df                <- table(cluster.vec[bait.gene]) %>% as.data.frame
            cluster.freq.df                <- cluster.freq.df[order(cluster.freq.df$Freq,decreasing = T),]
            value                          <- cluster.freq.df$Var1[1] %>% as.character %>% as.integer
            f1                             <- cluster.vec == value
            f2                             <- names(cluster.vec) %in% bait.gene
            sum(f1 & f2)/sum(f1)
        }
        bait.cluster.purity
    }
    
}




do.gene.fishing <- function(bait.gene, expr.matrix,alpha=5,fishing.round=1000){
    bait.gene <- bait.gene[bait.gene %in% colnames(expr.matrix)]
    all.gene  <- colnames(expr.matrix)      
    pool.gene <- setdiff(all.gene,bait.gene)
    tmp       <- foreach(i=1:fishing.round) %dopar%{
        shuffle.pool.gene <- sample(pool.gene,length(pool.gene))
        sub.pool.size     <- length(bait.gene) * alpha
        l <- foreach(j=seq(from=1,to=length(shuffle.pool.gene),by = sub.pool.size)) %do% {
            up <- min(j+sub.pool.size - 1,length(shuffle.pool.gene) )
            shuffle.pool.gene[j:up]          
        }  
        l
    }
    sub.pool.list <- unlist(tmp,recursive = F)
    fish.vec      <- foreach(sub.pool = sub.pool.list,.combine='c') %dopar% {
        bait.gene.expr.matrix          <- expr.matrix[,bait.gene]
        rd.gene                        <- sub.pool
        bait.and.pool.gene.expr.matrix <- cbind(expr.matrix[,bait.gene],expr.matrix[,rd.gene])
        bait.and.pool.gene.cor.matrix  <- cor(bait.and.pool.gene.expr.matrix,method='spearman')
        rs                             <- get.spectral.clustering.coordinates(bait.and.pool.gene.cor.matrix,method.name)
        if(rs[['no.of.cluster']] == 1){
            return('NOTHING')  
        }
        
        coordinates                    <- rs[['coordinates']]
        print(ncol(coordinates))
        kmeans.obj                     <- kmeanspp(data=coordinates,k=ncol(coordinates))
        cluster.vec                    <- kmeans.obj$cluster
        names(cluster.vec)             <- rownames(coordinates)
        cluster.freq.df                <- table(cluster.vec[bait.gene]) %>% as.data.frame
        cluster.freq.df                <- cluster.freq.df[order(cluster.freq.df$Freq,decreasing = T),]
        value                          <- cluster.freq.df$Var1[1] %>% as.character %>% as.integer
        f1                             <- cluster.vec == value
        f2                             <- names(cluster.vec) %in% bait.gene
        if(sum(f1 & !f2) ==0){
            return('NOTHING')
        }
        names(cluster.vec)[f1 & !f2]
    }
    fish.freq.df            <- table(fish.vec) %>% as.data.frame
    colnames(fish.freq.df)  <- c('fish.id','fish.freq')
    fish.freq.df$fish.freq  <- fish.freq.df$fish.freq / fishing.round
    fish.freq.df            <- arrange(fish.freq.df,desc(fish.freq))
    fish.freq.df            <- fish.freq.df[fish.freq.df$fish.id != 'NOTHING',]
    df                      <- data.frame(fish.id=setdiff(pool.gene.bait.gene),fish.freq=0)
    fish.freq.df            <- rbind(fish.freq.df,df)
    fish.freq.df 
}







