require(plyr)
require(dplyr)
require(foreach)
require(doParallel)
require(data.table)
require(peer)
registerDoParallel(20)

GTex.sample.annotation           <- read.delim("GTEx_Data_V6_Annotations_SampleAttributesDS.txt", stringsAsFactors=FALSE)
GTex.sample.annotation           <- GTex.sample.annotation[,c('SAMPID','SMTSD')]
colnames(GTex.sample.annotation) <- c('sample.id','tissue')
sample.no.per.tissue             <- table(GTex.sample.annotation$tissue) %>% as.data.frame
sample.no.per.tissue             <- sample.no.per.tissue[order(sample.no.per.tissue$Freq,decreasing = T),]
target.tissue                    <- sample.no.per.tissue$Var1[sample.no.per.tissue$Freq>=100] %>% as.character


gene.expr.data             <- fread('GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct',header=TRUE) %>% as.data.frame
gene.expr.data$Description <- NULL
gene.expr.data$Name        <- as.character(gene.expr.data$Name)
gene.expr.data$Name        <- sapply(gene.expr.data$Name, 
                             function(x)  { l <- strsplit(x,split = '\\.') %>% unlist  
                                            l[1]
                             }
                                    )
rownames(gene.expr.data) <- gene.expr.data$Name
gene.expr.data$Name      <- NULL
gene.expr.matrix         <- as.matrix(gene.expr.data) %>% t





#target.gene                      <- read.table('hg38.protein.coding.gene.id.txt',stringsAsFactors = F,header=F)[,1]
#gene.expr.matrix                 <- gene.expr.matrix[,colnames(gene.expr.matrix) %in% target.gene]


normalized.GTex.data    <- foreach(tissue = target.tissue) %dopar% {
    target.sample            <-  GTex.sample.annotation$sample.id[GTex.sample.annotation$tissue == tissue]   %>% as.character
    target.gene.expr.matrix  <-  gene.expr.matrix[rownames(gene.expr.matrix) %in% target.sample,]
    #tmp                      <-  apply(target.gene.expr.matrix,MARGIN=2,function(x) sum(x>=0.1))
    #target.gene.expr.matrix  <-  target.gene.expr.matrix[,tmp>=10]
    tmp                      <-  apply(target.gene.expr.matrix,MARGIN=2,function(x) median(x))
    target.gene.expr.matrix  <-  target.gene.expr.matrix[,tmp>=0.1]

    if(nrow(target.gene.expr.matrix) <= 100){
        matrix(1:10,nrow=2)
    }else{
        peer.obj                 <-  PEER()
        PEER_setAdd_mean(peer.obj,TRUE)
        PEER_setNk(peer.obj,25)
        PEER_setPhenoMean(peer.obj,log10(target.gene.expr.matrix+0.01))
        PEER_update(peer.obj)
        normalized.gene.expr.matrix            <- PEER_getResiduals(peer.obj)
        rownames(normalized.gene.expr.matrix)  <- rownames(target.gene.expr.matrix)
        colnames(normalized.gene.expr.matrix)  <- colnames(target.gene.expr.matrix)
        normalized.gene.expr.matrix
    }
}
names(normalized.GTex.data) <- target.tissue
save(file = 'normalized.GTex.with.ncRNA.median.rpkm.larger.than.0.1.RData',list=c('normalized.GTex.data','GTex.sample.annotation'))




