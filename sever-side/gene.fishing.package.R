require(foreach)
require(plyr)
require(dplyr)
require(ggplot2)
require(doParallel)
require(LICORS)
registerDoParallel(40)

# possible value of method.name:
# jordan
# random.walk
# berkeley.405

get.spectral.clustering.coordinates <- function(A,method.name){ 
    diag(A)  <-   0
    A        <-  abs(A)
    d        <-  apply(A,1,sum)
    I        <-  diag(1,nrow = nrow(A)) 
    
    if(method.name == 'random.walk' ){
        L    <-   I -  diag(1/d) %*% A 
    }else{
        L    <-   I -  diag(1/sqrt(d)) %*% A %*% diag(1/sqrt(d))
    }

    #L is positive semi-definite and have n non-negative real-valued eigenvalues
    tmp             <-  eigen(L)
    eigen.values    <-  tmp$values[length(tmp$values):1]
    eigen.vectors   <-  tmp$vectors[,length(tmp$values):1]
    
    q1              <-  quantile(eigen.values)[2]
    q3              <-  quantile(eigen.values)[4]
    low             <-  q1 - 3*(q3 - q1)
    up              <-  q1 + 3*(q3 - q1)
    k               <-  sum(eigen.values < low) # pick out eigen-values that are relatively small
    
    if(k == 1){
        list(eigen.values=eigen.values,no.of.cluster=1)
    }else{
       coordinate.matrix          <- eigen.vectors[,2:(k+1)] # pick out the eigen-vectors associated with the K smallest eigen-values
       if(method.name == 'jordan'){
           coordinate.matrix      <- apply(coordinate.matrix,MARGIN=1,function(x) x/ (sum(x*x) %>% sqrt) ) %>% t # re-normalize to an unit sphere
        }
        rownames(coordinate.matrix)  <- rownames(A)
        colnames(coordinate.matrix)  <- paste('dim',1:(k),sep="-")
        list(coordinates=coordinate.matrix,eigen.values=eigen.values,no.of.cluster=k)
    }
    
}



do.gene.fishing.v3 <- function(bait.gene, expr.matrix,alpha=5,fishing.round=1000,method.name='berkeley.405'){
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
    fish.freq.df
}




#######################
# do.gene.fishing.with.permutation <- function(bait.gene, expr.matrix,alpha=5,fishing.round=1000){
#     bait.gene <- bait.gene[bait.gene %in% colnames(expr.matrix)]
#     all.gene  <- colnames(expr.matrix)      
#     pool.gene <- setdiff(all.gene,bait.gene)
#     tmp       <- foreach(i=1:fishing.round) %dopar%{
#         shuffle.pool.gene <- sample(pool.gene,length(pool.gene))
#         sub.pool.size     <- length(bait.gene) * alpha
#         l <- foreach(j=seq(from=1,to=length(shuffle.pool.gene),by = sub.pool.size)) %do% {
#             up <- min(j+sub.pool.size - 1,length(shuffle.pool.gene) )
#             shuffle.pool.gene[j:up]          
#         }  
#         l
#     }
#     sub.pool.list <- unlist(tmp,recursive = F)
#     fish.vec      <- foreach(sub.pool = sub.pool.list,.combine='c') %dopar% {
#         bait.gene.expr.matrix          <- expr.matrix[,bait.gene]
#         rd.gene                        <- sub.pool
#         rd.matrix                      <- expr.matrix[,rd.gene]
#         rd.matrix                      <- apply(rd.matrix,MARGIN=2,sample)
#         bait.and.pool.gene.expr.matrix <- cbind(expr.matrix[,bait.gene],rd.matrix)
#         bait.and.pool.gene.cor.matrix  <- cor(bait.and.pool.gene.expr.matrix,method='spearman')
#         rs                             <- get.spectral.clustering.coordinates(bait.and.pool.gene.cor.matrix)
#         coordinates                    <- rs[['coordinates']]
#         eigen.value                    <- rs[['eigen.value']]
#         
# #         #eigen.value                    <- eigen.value[-1] # remove the first eigen-value, it is always zero
# #         #eigen.value.delta              <- eigen.value[2:length(eigen.value)] - eigen.value[1:(length(eigen.value)-1)]
# #         #if(eigen.value.delta[1] != max(eigen.value.delta)){
# #         #     return('NOTHING')
# #         # }  
#         tmp            <- eigen.value
#         q1             <- quantile(tmp)[2]
#         q3             <- quantile(tmp)[4]
#         low            <- q1 - 3*(q3 - q1)
#         up             <- q1 + 3*(q3 - q1)
#         no.of.cluster  <- sum(eigen.value < low)
#         if( no.of.cluster == 1){
#             return('NOTHING')
#         }  
#         
#         
#          #K-Means clustering with clustering centers       
# #         center.matrix                  <- rbind(
# #           c(median(coordinates[bait.gene,1]),median(coordinates[bait.gene,2])),
# #           c(median(coordinates[rd.gene,  1]),median(coordinates[rd.gene,  2]))                                 
# #         )
# #         kmeans.obj                     <- kmeans(x=coordinates,centers=center.matrix)
#         kmeans.obj                      <- kmeanspp(data=coordinates,k=no.of.cluster)
#         
#         cluster.vec                    <- kmeans.obj$cluster
#         names(cluster.vec)             <- rownames(coordinates)
#         cluster.freq.df                <- table(cluster.vec[bait.gene]) %>% as.data.frame
#         cluster.freq.df                <- cluster.freq.df[order(cluster.freq.df$Freq,decreasing = T),]
#         value                          <- cluster.freq.df$Var1[1] %>% as.character %>% as.integer
#         f1                             <- cluster.vec == value
#         f2                             <- names(cluster.vec) %in% bait.gene
#         if(sum(f1 & !f2) ==0){
#             return('NOTHING')
#         }
#         names(cluster.vec)[f1 & !f2]
#     }
#     fish.freq.df            <- table(fish.vec) %>% as.data.frame
#     colnames(fish.freq.df)  <- c('fish.id','fish.freq')
#     fish.freq.df$fish.freq  <- fish.freq.df$fish.freq / fishing.round
#     fish.freq.df
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# explore.fishing.potential <- function(bait.gene, expr.matrix,alpha=5){
#   bait.gene             <- bait.gene[bait.gene %in% colnames(expr.matrix)]
#   all.gene              <- colnames(expr.matrix)      
#   pool.gene             <- setdiff(all.gene,bait.gene)
#   bait.gene.expr.matrix <- expr.matrix[,bait.gene]
#   pool.gene.expr.matrix <- expr.matrix[,pool.gene]
#   
#   rd.gene                        <- sample(colnames(pool.gene.expr.matrix), length(bait.gene) * alpha  )
#   bait.and.pool.gene.expr.matrix <- cbind(bait.gene.expr.matrix,
#                                           pool.gene.expr.matrix[,rd.gene]
#   )
#   bait.and.pool.gene.cor.matrix  <- cor(bait.and.pool.gene.expr.matrix,method='spearman')
#   rs                             <- get.spectral.clustering.coordinates(bait.and.pool.gene.cor.matrix)
#   coordinates                    <- rs$vectors[,(ncol(bait.and.pool.gene.cor.matrix)-2):(ncol(bait.and.pool.gene.cor.matrix)-1)]
#   rownames(coordinates)          <- rownames(bait.and.pool.gene.cor.matrix)
#   colnames(coordinates)          <- c('eigen2','eigen1')
#   
#   bait.distance                  <- dist(x = coordinates[bait.gene,]) %>% as.matrix 
#   radius                         <- quantile(bait.distance[upper.tri(bait.distance)])[4]/2
#   center.x                       <- median(coordinates[bait.gene,1])
#   center.y                       <- median(coordinates[bait.gene,2])
#   
#   
#   
#   center.matrix                  <- rbind(c(median(coordinates[bait.gene,1]),median(coordinates[bait.gene,2])),
#                                           c(median(coordinates[rd.gene,  1]),median(coordinates[rd.gene,  2]))                                 
#   )
#   kmeans.obj                     <- kmeans(x=coordinates,centers=center.matrix)
#   cluster.vec                    <- kmeans.obj$cluster
#   draw.df                        <- data.frame(eigen1 = coordinates[,2],eigen2=coordinates[,1],cluster = cluster.vec) 
#   rownames(draw.df)              <- rownames(coordinates)
#   draw.df$gene.type                      <- 'pool'
#   draw.df$gene.type[1:length(bait.gene)] <- 'bait'
#   
#   
#   p <- ggplot(draw.df) + 
#     geom_point(aes(x=eigen1 ,y=eigen2 ,color=gene.type,shape=factor(cluster)),size=5) + 
#     geom_path(data = circle(center.x,center.y,radius=radius),aes(x=y,y=x) ) +
#     theme(text = element_text(size=20))
#   p
#   
#   eigen1        <- draw.df$eigen1
#   names(eigen1) <- rownames(draw.df)
#   eigen1        <- sort(eigen1)
#   image(bait.and.pool.gene.cor.matrix[names(eigen1),names(eigen1)])
#   list(p)
# }
# 
# 
# 
# 
# ########### code after this line  is dead!###############
# 
# # circle <- function(center_x=0,center_y=0, radius  = 1,npoints = 100) {
# #   cycle <- seq(0,2*pi,length.out = npoints)
# #   xx    <- center_x + radius * cos(cycle)
# #   yy    <- center_y + radius * sin(cycle)
# #   return(data.frame(x = xx, y = yy))
# # }
# 
# 
# 
# 
# 
# 
# # fishing.once <- function(bait.gene, expr.matrix,alpha=5,ep.vecotr){
# #   bait.gene             <- bait.gene[bait.gene %in% colnames(expr.matrix)]
# #   all.gene              <- colnames(expr.matrix)      
# #   pool.gene             <- setdiff(all.gene,bait.gene)
# #   bait.gene.expr.matrix <- expr.matrix[,bait.gene]
# #   pool.gene.expr.matrix <- expr.matrix[,pool.gene]
# #   
# #   rd.gene                        <- sample(colnames(pool.gene.expr.matrix), length(bait.gene) * alpha  )
# #   bait.and.pool.gene.expr.matrix <- cbind(bait.gene.expr.matrix,
# #                                           pool.gene.expr.matrix[,rd.gene]
# #   )
# #   bait.and.pool.gene.cor.matrix  <- cor(bait.and.pool.gene.expr.matrix,method='spearman')
# #   rs                             <- get.spectral.clustering.coordinates(bait.and.pool.gene.cor.matrix)
# #   coordinates                    <- rs$vectors[,(ncol(bait.and.pool.gene.cor.matrix)-2):(ncol(bait.and.pool.gene.cor.matrix)-1)]
# #   rownames(coordinates)          <- rownames(bait.and.pool.gene.cor.matrix)
# #   colnames(coordinates)          <- c('eigen2','eigen1')
# #   
# #   d1 <- (dist(coordinates[bait.gene,])   %>% as.matrix %>% sum)
# #   d2 <- (dist(coordinates[rd.gene  ,])   %>% as.matrix %>% sum)
# #   d  <- (dist(coordinates)               %>% as.matrix %>% sum)
# #   
# #   
# #   ratio   <- (d1)/d
# #   p.value <- sum(ep.vector <= ratio)/length(ep.vector)
# #   p.value
# # }
# 
# 
# # get.empirical.distance.ratio.distribution <- function(expr.matrix,bait.gene.number,alpha){
# #   rs <- foreach(i=1:100000,.combine='c') %dopar% {
# #     g                     <- sample(x = colnames(expr.matrix),size=bait.gene.number * (1 + alpha))
# #     cor.matrix            <- cor(expr.matrix[,g],method='spearman')
# #     rs                    <- get.spectral.clustering.coordinates(cor.matrix)
# #     coordinates           <- rs$vectors[,(ncol(cor.matrix)-2):(ncol(cor.matrix)-1)]
# #     bait.gene.coordinates <- coordinates[1:bait.gene.number,]
# #     rd.gene.coordinates <- coordinates[(bait.gene.number + 1):nrow(coordinates),]
# #     
# #     d1 <- (dist(bait.gene.coordinates)   %>% as.matrix %>% sum)
# #     d2 <- (dist(  rd.gene.coordinates)   %>% as.matrix %>% sum)
# #     d  <- (dist(coordinates)             %>% as.matrix %>% sum)
# #     
# #     ratio                 <- (d1)/d
# #     ratio    
# #   } 
# #   rs
# # }
# #In this version of gene fishing, we first judege whether the 21 bait gens
# #are close enough to each other before clustering. If they just randomly distribute in the eigen2 vs eigen1 space,
# #we do nothing!
# # do.gene.fishing.v2 <- function(bait.gene, expr.matrix,alpha=5,fishing.round=1000,ep.vector){
# #   bait.gene <- bait.gene[bait.gene %in% colnames(expr.matrix)]
# #   all.gene  <- colnames(expr.matrix)      
# #   pool.gene <- setdiff(all.gene,bait.gene)
# #   tmp       <- foreach(i=1:fishing.round) %dopar%{
# #     shuffle.pool.gene <- sample(pool.gene,length(pool.gene))
# #     sub.pool.size     <- length(bait.gene) * alpha
# #     l                 <- foreach(j=seq(from=1,to=length(shuffle.pool.gene),by = sub.pool.size)) %do% {
# #       up <- min(j+sub.pool.size - 1,length(shuffle.pool.gene) )
# #       shuffle.pool.gene[j:up]          
# #     }  
# #     l
# #   }
# #   sub.pool.list <- unlist(tmp,recursive = F)
# #   fish.vec      <- foreach(sub.pool = sub.pool.list,.combine='c') %dopar% {
# #     bait.gene.expr.matrix <- expr.matrix[,bait.gene]
# #     rd.gene               <- sub.pool
# #     bait.and.pool.gene.expr.matrix <- cbind(expr.matrix[,bait.gene],
# #                                             expr.matrix[,rd.gene]
# #     )
# #     bait.and.pool.gene.cor.matrix  <- cor(bait.and.pool.gene.expr.matrix,method='spearman')
# #     rs                             <- get.spectral.clustering.coordinates(bait.and.pool.gene.cor.matrix)
# #     coordinates                    <- rs$vectors[,(ncol(bait.and.pool.gene.cor.matrix)-2):(ncol(bait.and.pool.gene.cor.matrix)-1)]
# #     rownames(coordinates)          <- c(bait.gene,rd.gene)
# #     ratio                          <- (dist(coordinates[bait.gene,]) %>% as.matrix %>% sum) / (dist(coordinates) %>% as.matrix %>% sum)
# #     p.value                        <- sum(ep.vector<=ratio) / length(ep.vector)
# #     if(p.value>0.05){
# #       return('NOTHING')
# #     }
# #     center.matrix                  <- rbind(c(median(coordinates[bait.gene,1]),median(coordinates[bait.gene,2])),
# #                                             c(median(coordinates[rd.gene,  1]),median(coordinates[rd.gene,  2]))                                 
# #     )
# #     kmeans.obj                     <- kmeans(x=coordinates,centers=center.matrix)
# #     cluster.vec                    <- kmeans.obj$cluster
# #     names(cluster.vec)             <- rownames(coordinates)
# #     cluster.freq.df                <- table(cluster.vec[bait.gene]) %>% as.data.frame
# #     cluster.freq.df                <- cluster.freq.df[order(cluster.freq.df$Freq,decreasing = T),]
# #     value                          <- cluster.freq.df$Var1[1] %>% as.character %>% as.integer
# #     f1                             <- cluster.vec == value
# #     f2                             <- names(cluster.vec) %in% bait.gene
# #     if(sum(f1 & !f2) ==0){
# #       return('NOTHING')
# #     }
# #     names(cluster.vec)[f1 & !f2]
# #     
# #   }
# #   #fish.vec               <- fish.vec[fish.vec != 'NOTHING']
# #   fish.freq.df           <- table(fish.vec) %>% as.data.frame
# #   colnames(fish.freq.df) <- c('fish.id','fish.freq')
# #   fish.freq.df$fish.freq <- fish.freq.df$fish.freq / fishing.round
# #   fish.freq.df
# # }
# 
# 
# 
# 
# do.gene.fishing <- function(bait.gene, expr.matrix,alpha=5,fishing.round=1000){
#   bait.gene <- bait.gene[bait.gene %in% colnames(expr.matrix)]
#   all.gene  <- colnames(expr.matrix)      
#   pool.gene <- setdiff(all.gene,bait.gene)
#   tmp       <- foreach(i=1:fishing.round) %dopar%{
#     shuffle.pool.gene <- sample(pool.gene,length(pool.gene))
#     sub.pool.size     <- length(bait.gene) * alpha
#     l <- foreach(j=seq(from=1,to=length(shuffle.pool.gene),by = sub.pool.size)) %do% {
#       up <- min(j+sub.pool.size - 1,length(shuffle.pool.gene) )
#       shuffle.pool.gene[j:up]          
#     }  
#     l
#   }
#   sub.pool.list <- unlist(tmp,recursive = F)
#   fish.vec <- foreach(sub.pool = sub.pool.list,.combine='c') %dopar% {
#     bait.gene.expr.matrix <- expr.matrix[,bait.gene]
#     rd.gene               <- sub.pool
#     bait.and.pool.gene.expr.matrix <- cbind(expr.matrix[,bait.gene],
#                                             expr.matrix[,rd.gene]
#     )
#     bait.and.pool.gene.cor.matrix  <- cor(bait.and.pool.gene.expr.matrix,method='spearman')
#     rs                             <- get.spectral.clustering.coordinates(bait.and.pool.gene.cor.matrix)
#     coordinates                    <- rs$vectors[,(ncol(bait.and.pool.gene.cor.matrix)-2):(ncol(bait.and.pool.gene.cor.matrix)-1)]
#     rownames(coordinates)          <- c(bait.gene,rd.gene)
#     ratio                          <- (dist(coordinates[bait.gene,]) %>% as.matrix %>% sum) / (dist(coordinates) %>% as.matrix %>% sum)
#     
#     
#     
#     
#     #  K-Means clustering        
#     center.matrix                  <- rbind(
#       c(median(coordinates[bait.gene,1]),median(coordinates[bait.gene,2])),
#       c(median(coordinates[rd.gene,  1]),median(coordinates[rd.gene,  2]))                                 
#     )
#     kmeans.obj                     <- kmeans(x=coordinates,centers=center.matrix)
#     cluster.vec                    <- kmeans.obj$cluster
#     names(cluster.vec)             <- rownames(coordinates)
#     cluster.freq.df                <- table(cluster.vec[bait.gene]) %>% as.data.frame
#     cluster.freq.df                <- cluster.freq.df[order(cluster.freq.df$Freq,decreasing = T),]
#     value                          <- cluster.freq.df$Var1[1] %>% as.character %>% as.integer
#     f1                             <- cluster.vec == value
#     f2                             <- names(cluster.vec) %in% bait.gene
#     if(sum(f1 & !f2) ==0){
#       return('NOTHING')
#     }
#     names(cluster.vec)[f1 & !f2]
#     
#   }
#   #fish.vec               <- fish.vec[fish.vec != 'NOTHING']
#   fish.freq.df            <- table(fish.vec) %>% as.data.frame
#   colnames(fish.freq.df)  <- c('fish.id','fish.freq')
#   fish.freq.df$fish.freq  <- fish.freq.df$fish.freq / fishing.round
#   fish.freq.df
# }
# 



