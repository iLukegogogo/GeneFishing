require(plyr)
require(dplyr)
require(foreach)
require(doParallel)
load('normalized.GTex.RData')



setwd('./GTEx_Analysis_v6p_eQTL_covariates')
cmd       <- 'ls |grep txt'
file.list <- system(cmd,intern = T,wait=T)
gender.df <- foreach(file= file.list,.combine='rbind') %do% {
    df          <- read.table(file,header=T,stringsAsFactors = F)
    df$ID       <- NULL
    gender.vec  <- df[nrow(df)-1,] %>% unlist
    individual.id   <- colnames(df)
    data.frame(individual.id=individual.id,gender=gender.vec)
}
gender.df <- ddply(gender.df,.(individual.id),function(x) x$gender[1])
colnames(gender.df) <- c('individual.id','gender')
rownames(gender.df) <- gender.df$individual.id
setwd('../')

pca.gene <- read.table('PCA.gene.txt',header=F,stringsAsFactors = F)[,1] %>% as.character
pdf('pca.pdf')
for(j in 1:43){
  expr.matrix <- normalized.GTex.data[[j]]
  expr.matrix <- expr.matrix[,colnames(expr.matrix) %in% pca.gene]
  if(nrow(expr.matrix) <100){
    next
  }
  pca.result <- prcomp(expr.matrix,center = T,retx = TRUE,scale. = FALSE)
  individual.id <- sapply(rownames(expr.matrix),function(x){
      l <- strsplit(x = x,split='-') %>% unlist  
      paste(l[1],l[2],sep=".")
  })
  gender <- gender.df[individual.id,'gender']
  plot(pca.result$x[,1:2],xlab='PC1',ylab='PC2',main=names(normalized.GTex.data)[j],col=c('red','green')[gender %>% as.factor])
}
dev.off()