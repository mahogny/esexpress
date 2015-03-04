
### todo use esexpress schema
library("DESeq")
library("DESeq2")
source("common.R")


normalizeDeseq <- function(cnt) t(t(cnt)/estimateSizeFactorsForMatrix(cnt))





####################################################################################
## Read ola counts from all the files
readcntola <- function(f){
  dat <- read.table(f,header=TRUE)
  rownames(dat) <- dat[,1]
  return(dat[,-1])
}
  
# cnt_es<-
#   cbind(
#     readcntola("../../data/ola_mES_2i_2.txt"),
#     readcntola("../../data/ola_mES_2i_3.txt"),
#     readcntola("../../data/ola_mES_2i_4.txt"),
#     readcntola("../../data/ola_mES_2i_5.txt"),
#     readcntola("../../data/ola_mES_a2i_2.txt"),
#     readcntola("../../data/ola_mES_a2i_3.txt"),
#     readcntola("../../data/ola_mES_lif_1.txt"),
#     readcntola("../../data/ola_mES_lif_2.txt"),
#     readcntola("../../data/ola_mES_lif_3.txt"))

cnt_es<-read.table("../../data/n_counts_mES.txt",header=TRUE)



#remove dead cells
#removedcells <- read.table("../../data/removed_cells.txt",header=FALSE,sep=" ",stringsAsFactors=FALSE)[,2]
#cnt_es <- round(cnt_es[,-(which(colnames(cnt_es) %in% removedcells))])

#save new table
#write.table(cnt_es,"static/counttable_es.csv")

#remove ERCs and other crap. no longer needed
#cnt_es <- cnt_es[-grep("ERCC",rownames(cnt_es)),]
#cnt_es <- cnt_es[-grep("_",rownames(cnt_es)),]

#normalize. no longer needed
cnt_es_notnorm <- cnt_es
cnt_es <- normalizeDeseq(cnt_es)

genesub<-rownames(cnt_es)#[1:4000]
ids<-which(rownames(cnt_es) %in% genesub)
ds_ola_lif<-cnt_es[ids,grep("lif",colnames(cnt_es))]
ds_ola_2i <-cnt_es[ids,grep("_2i_",colnames(cnt_es))]
ds_ola_a2i<-cnt_es[ids,grep("_a2i_",colnames(cnt_es))]



cellstates <- read.csv("../../data/cell_states.txt",sep=" ",stringsAsFactors=FALSE)[,c(2,3)]
conditions <- merge(data.frame(cell=colnames(cnt_es),stringsAsFactors = FALSE), cellstates)$state

set_2i=c("2i","2C")
set_nanog=c("Nanog_hi","Nanog_med","Nanog_lo")

ds_ola_2i_2C<-cnt_es[ids,conditions=="2C"]
ds_ola_2i_2i<-cnt_es[ids,conditions=="2i"]
ds_ola_nanog_hi <-cnt_es[ids,conditions=="Nanog_hi"]
ds_ola_nanog_med<-cnt_es[ids,conditions=="Nanog_med"]
ds_ola_nanog_lo <-cnt_es[ids,conditions=="Nanog_lo"]


####################################################################################


## Count tables for sandberg
datsand<-read.table("../../data/Sandberg_data.txt",sep=" ", stringsAsFactors = FALSE)[,-(1:7)]
nsand <- normalizeDeseq(datsand)

ids<-which(rownames(nsand) %in% genesub)
ds_s_lblast<-nsand[ids,grep("lateblast",colnames(nsand))]
ds_s_mblast<-nsand[ids,grep("midblast",colnames(nsand))]
ds_s_eblast<-nsand[ids,grep("earlyblast",colnames(nsand))]



####################################################################################
## Count tables
uploadcounts <- function(dataset, cnt){  
  dbGetQuery(con,sprintf("delete from geneexp where dataset='%s';",dataset))
  geneexp<-cbind(dataset,row.names(cnt),encodearrayS(colnames(cnt)),apply(cnt,1,encodearray))
  colnames(geneexp)<-c("dataset","fromgene","fromcell","exp")
  for(i in 1:nrow(geneexp)){
    if(i %% 100 == 0)
      print(i)
    v<-sprintf("insert into geneexp values ('%s','%s',array[%s],array[%s])",geneexp[i,1],geneexp[i,2],geneexp[i,3],geneexp[i,4])
    dbGetQuery(con,v)
  }
}









#unique(dbReadTable(con,"geneexp")[,1])





####################################################################################
## Gene-gene correlations
uploadcorr <- function(dataset, cnt){
  print(sprintf("# genes in cor matrix: %s",ncol(cnt)))
  thec <- cor(t(cnt),method="spearman")
  for(i in 1:nrow(thec)){
    thec[which(is.na(thec))]<-0
  }
  #colnames(linc)<-c("dataset","fromgene","togene","corr","corrp","corrv","corrm")
  togene<-encodearrayS(colnames(thec))
  dbGetQuery(con,sprintf("delete from genecorr WHERE dataset='%s';",dataset))
  for(i in 1:nrow(thec)){
    
    onecor <- encodearray(thec[i,])
    pval <- onecor
    fromgene <- rownames(thec)[i]
    
    if(i %% 100 == 0)
      print(sprintf("%s of %s",i,nrow(thec)))
    
    v<-sprintf("insert into genecorr values ('%s','%s',array[%s],array[%s],array[%s],array[%s],array[%s])",
      dataset,fromgene, togene, onecor, pval, pval, pval)
    dbGetQuery(con,v)
  }
}




## Pull out a subset for correlation
#Reduce to around 5000 genes. it should be the ones with a few cells having some expr
getexpressed <- function(ds){
  foo <- apply(ds_ola_lif>20,1,function(x) length(which(x))>3)    #30 for 600 genes
  ds[which(rownames(ds) %in% names(foo)[which(foo)]),]
}
#red_ds_ola_2i <- getexpressed(ds_ola_2i)
#nrow(red_ds_ola_2i)

#foo <- apply(ds_ola_lif>20,1,function(x) length(which(x))>3)    #30 for 600 genes
#which(rownames(ds_s_eblast) %in% names(foo)[which(foo)])

#check that all genes on basal list are there!  >20 & >3 is ok
#plugenes <- read.table("genelist_pluripotency.txt",stringsAsFactors = FALSE)[,1]
#plugenes[-(which(mapsymid(plugenes) %in% rownames(red_ds_ola_2i)))]



#ncol(es_cnt)  ##check, 704???

#dbGetQuery(con,sprintf("delete from genecorr;"))


#uploadcorr("es_lif",ds_ola_lif)






