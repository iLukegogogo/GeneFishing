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
    #v3               <- rs$vectors[,nrow(cor.matrix) - 3]
    m                <-  cbind(v1,v2)
    #m                  <-  apply(m,MARGIN=1,function(x) x/sqrt(sum(x*x))) %>% t
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
    if(sum(cluster.vec == value) >= 35){
        file.name  <- paste('tmp/',runif(1),".pdf",sep="")
        pdf(file.name)
        plot(m,col=c('red','black')[as.factor(cluster.vec)])
        dev.off()
    }

    f1                 <- cluster.vec == value
    f2                 <- names(cluster.vec) %in% bait.gene.name
    names(cluster.vec)[f1 & !f2] #our real fish
}




registerDoParallel(40)
bait.matrix       <- PEER.expr.matrix[,c1] # the 21 genes as bait
flag              <- colnames(PEER.expr.matrix) %in% c1
rd.gene.name      <- colnames(PEER.expr.matrix)[!flag]



total.gene.number <- 1000
other.gene.number <- total.gene.number - length(c1)
idx.vec           <- seq(from =1,to=length(rd.gene.name),by =other.gene.number) # 

all.fish.vec      <- foreach(k=1:1000,.combine='c') %do%{
    rd.gene.name  <- rd.gene.name[sample(length(rd.gene.name))]
    fish.vec <- foreach(i = idx.vec,.combine='c') %dopar%{
        j              <- min(c(i+other.gene.number-1,length(rd.gene.name)))
        pool.gene.name <- c(rd.gene.name[i:j],c1) %>% unique
        fishing.once(PEER.expr.matrix[,pool.gene.name],c1)
    }
    print(length(fish.vec))
    fish.vec
   
}

fish.freq        <- table(all.fish.vec)
all.fish.vec     <- all.fish.vec %>% unique
cholesterol.gene <- read.table("original_data/cholesterol.gene.list", quote="\"", stringsAsFactors=FALSE)$V1
cholesterol.gene <- cholesterol.gene[cholesterol.gene %in% colnames(PEER.expr.matrix)]


new.bait.matrix <- PEER.expr.matrix[,c(all.fish.vec,cholesterol.gene) %>% unique] 

save(file='RData/new.Peter.gene.fishing.1000.RData',list = c('new.bait.matrix','all.fish.vec','fish.freq'))



