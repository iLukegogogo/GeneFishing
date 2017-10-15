require(plyr)
require(dplyr)
require(foreach)
require(parallel)
require(doParallel)
require(reshape2)
registerDoParallel(30)


the.56.gene <- read.table(file='the.56.genes.tsv',header=F,stringsAsFactors = F)$V1


remove.dot.suffix <- function(x){
  unlist(strsplit(x = x,split = "\\."))[1]   
  
}
setwd('pairs')
eQTL.file.list <- system('ls',intern = T,wait = T)
df <- foreach(file.name = eQTL.file.list,.combine='rbind') %dopar% {
    eQTL.df    <- read.table(file=file.name,header=T,stringsAsFactors = F)
    tissue    <- gsub(x = file.name,pattern="_Analysis.v6p.signif_snpgene_pairs.txt",replacement = "") 
    rs        <- eQTL.df[,c(1,2,4)]
    rs[,2]    <- sapply(rs[,2] %>% as.character, remove.dot.suffix)
    rs$tissue <- tissue
    rs        <- rs[rs[,2] %in% the.56.gene,]
    rs
}
new.df  <- dcast(df,formula=variant_id + gene_id ~ tissue,value.var='pval_nominal')


p.value.str.vec <- foreach(i=1:nrow(new.df),.combine='c') %dopar% {
    row <- new.df[i,]
    tissue  <- colnames(new.df)[3:ncol(new.df)]
    p.value <- row[,3:ncol(new.df)] %>% unlist %>% c
    p.value.str <- paste(tissue,p.value,sep=":",collapse="@")
    p.value.str
}
new.df               <- new.df[,1:2]
new.df$eQTL.p.value  <- p.value.str.vec
write.table(x=new.df,file='../GTex.eQTL.gene.assiciation.txt',quote=F,col.names=F,row.names=F)



setwd('../')
GTex.eQTL           <- read.table(file='GTex.eQTL.gene.assiciation.txt',header = F,stringsAsFactors = F)
colnames(GTex.eQTL) <- c('V1','target.gene.id','eQTL.p.value')
bed.lines <- foreach(i=1:nrow(GTex.eQTL),.combine='c') %dopar% {
    pos.str        <-   GTex.eQTL$V1[i]
    tmp            <-   strsplit(split="_",x=pos.str) %>% unlist
    new.pos.str    <-   sprintf("chr%s:%s",as.character(tmp[1]),as.character(tmp[2]))
    target.gene.id <-   GTex.eQTL$target.gene.id[i]
    eQTL.p.value   <-   GTex.eQTL$eQTL.p.value[i]
    bed.str        <-   paste("chr",tmp[1],"\t",as.integer(tmp[2])-1,"\t",tmp[2],"\t",new.pos.str,"\t",target.gene.id,"\t",eQTL.p.value,sep="")
    bed.str
}
write.table(x=unique(bed.lines),row.names=F,col.names=F,quote=F,file='GTex.eQTL.bed')
system("bedtools intersect -loj -a GTex.eQTL.bed -b dbSNP147.bed |cut -f 1,2,3,4,5,10,6 > GTex.eQTL.with.rsid.bed")

