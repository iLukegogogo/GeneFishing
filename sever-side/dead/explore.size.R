load('RData//util_data.RData')
source('code/do.laplacian.R')
source('code/include.R')
require(LICORS)

PEER.raw.data              <- read.delim("original_data/expression/426sham0covK40addMeanPEERresidualsTranspose.txt", stringsAsFactors=FALSE)
rownames(PEER.raw.data)    <- PEER.raw.data$Gene
PEER.raw.data$Gene         <- NULL
PEER.raw.data.mt           <- t(as.matrix(PEER.raw.data))
PEER.expr.matrix           <- PEER.raw.data.mt[rownames(PEER.expr.matrix),]


c1 <- read.table("Results/c1.csv", quote="\"", stringsAsFactors=FALSE)$V1


fishing.once <- function(pool.matrix,bait.gene.name){
  
  cor.matrix       <- cor(pool.matrix,method='spearman')
  rs               <- do.laplacian(cor.matrix)
  v1               <- rs$vectors[,nrow(cor.matrix) - 1]
  v2               <- rs$vectors[,nrow(cor.matrix) - 2]
  m                <-  cbind(v1,v2)
  names(v1)        <- rownames(cor.matrix)
  names(v2)        <- rownames(cor.matrix)
  center.matrix    <- rbind(c(median(v1[bait.gene.name]),median(v2[bait.gene.name])),
                            c(median(v1),median(v2))
  )
  if(is.na(center.matrix[1])){
    print(v1[bait.gene.name])
  }
  kmeans.obj         <- kmeans(x=m,centers=center.matrix)
  #kmeans.obj         <-  kmeanspp(data = m,k = 2)
  cluster.vec        <-  kmeans.obj$cluster
  names(cluster.vec) <-  rownames(cor.matrix)
  
  cluster.freq.df    <-  table(cluster.vec[bait.gene.name]) %>% as.data.frame
  cluster.freq.df    <-  cluster.freq.df[order(cluster.freq.df$Freq,decreasing = T),]
  value              <-  cluster.freq.df$Var1[1] %>% as.character %>% as.integer
  #print(cluster.freq.df)
  #print(value)
  #print(sum(cluster.vec == value))
  #if(sum(cluster.vec == value) >= 35){
  #  file.name  <- paste('tmp/',runif(1),".pdf",sep="")
  #  pdf(file.name)
  #  plot(m,col=c('red','black')[as.factor(cluster.vec)])
  #  dev.off()
  #}
  
  f1  <- cluster.vec == value
  sum(f1)
}

registerDoParallel(40)
df <- foreach(total.gene.number = seq(from=400,to=5000,by=100),.combine='rbind') %do%{
    vec <- foreach(i=1:1000,.combine='c') %dopar%{
        rd.gene   <- sample(x=colnames(PEER.expr.matrix),size=total.gene.number - 21)  
        pool.gene <- c(c1,rd.gene) %>% unique
        fishing.once(PEER.expr.matrix[,pool.gene],c1)
    }
    data.frame(gene.num=total.gene.number,cluster.size=vec)
  
}
save(file='RData/explore.size.RData',list=c('df'))



