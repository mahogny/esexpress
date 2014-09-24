
### todo use esexpress schema
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
  
cnt_es<-
  cbind(
    readcntola("../../data/ola_mES_2i_2.txt"),
    readcntola("../../data/ola_mES_2i_3.txt"),
    readcntola("../../data/ola_mES_2i_4.txt"),
    readcntola("../../data/ola_mES_2i_5.txt"),
    readcntola("../../data/ola_mES_a2i_2.txt"),
    readcntola("../../data/ola_mES_a2i_3.txt"),
    readcntola("../../data/ola_mES_lif_1.txt"),
    readcntola("../../data/ola_mES_lif_2.txt"),
    readcntola("../../data/ola_mES_lif_3.txt"))

#remove ERCs and other crap
cnt_es <- cnt_es[-grep("ERCC",rownames(cnt_es)),]
cnt_es <- cnt_es[-grep("_",rownames(cnt_es)),]

#remove dead cells, normalize
removedcells <- read.table("../../data/removed_cells.txt",header=FALSE,sep=" ",stringsAsFactors=FALSE)[,2]
cnt_es <- cnt_es[,-(which(colnames(cnt_es) %in% removedcells))]
cnt_es_notnorm <- cnt_es
cnt_es <- normalizeDeseq(cnt_es)




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

#cnt_es<-read.table("../../n_counts_mES.txt",header=TRUE)

genesub<-rownames(cnt_es)#[1:4000]

ids<-which(rownames(cnt_es) %in% genesub)
ds_ola_lif<-cnt_es[ids,grep("lif",colnames(cnt_es))]
ds_ola_2i <-cnt_es[ids,grep("_2i_",colnames(cnt_es))]
ds_ola_a2i<-cnt_es[ids,grep("_a2i_",colnames(cnt_es))]

uploadcounts("es_lif",ds_ola_lif)
uploadcounts("es_2i", ds_ola_2i)
uploadcounts("es_a2i",ds_ola_a2i)
#dbReadTable(con,"geneexp")[,2]


####################################################################################
## Count tables for sandberg
datsand<-read.table("../../data/Sandberg_data.txt",sep=" ", stringsAsFactors = FALSE)[,-(1:7)]
nsand <- normalizeDeseq(datsand)


ids<-which(rownames(nsand) %in% genesub)
ds_s_lblast<-nsand[ids,grep("lateblast",colnames(nsand))]
ds_s_mblast<-nsand[ids,grep("midblast",colnames(nsand))]
ds_s_eblast<-nsand[ids,grep("earlyblast",colnames(nsand))]

uploadcounts("sandberg_lateblast",ds_s_lblast)
uploadcounts("sandberg_midblast", ds_s_mblast)
uploadcounts("sandberg_earlyblast",ds_s_eblast)


#unique(dbReadTable(con,"geneexp")[,1])






####################################################################################
## Gene-gene correlations
uploadcorr <- function(dataset, cnt){
  
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
  ds[which(foo),]
}
#red_ds_ola_2i <- getexpressed(ds_ola_2i)
#nrow(red_ds_ola_2i)



#check that all genes on basal list are there!  >20 & >3 is ok
plugenes <- read.table("genelist_pluripotency.txt",stringsAsFactors = FALSE)[,1]
plugenes[-(which(mapsymid(plugenes) %in% rownames(red_ds_ola_2i)))]



#ncol(es_cnt)  ##check, 704???

uploadcorr("es_lif",getexpressed(ds_ola_lif))
uploadcorr("es_2i", getexpressed(ds_ola_2i))
uploadcorr("es_a2i",getexpressed(ds_ola_a2i))
uploadcorr("sandberg_earlyblast",getexpressed(ds_s_eblast))
uploadcorr("sandberg_midblast",  getexpressed(ds_s_mblast))
uploadcorr("sandberg_lateblast", getexpressed(ds_s_lblast))

#uploadcorr("es_lif",ds_ola_lif)




####################################################################################
## Gene info
geneinfo <- read.csv("../../data/ensembl2genename.txt",stringsAsFactors=FALSE)
geneinfo <- cbind(geneinfo,0)
colnames(geneinfo)<-c("geneid","genesym","pvalbiovar")
geneinfo<-as.data.frame(geneinfo,stringsAsFactors=FALSE)
geneinfo <- geneinfo[which(geneinfo$geneid %in% rownames(ds_ola_lif)),]

dbGetQuery(con,"delete from geneinfo;")
dbWriteTable(con,"geneinfo",geneinfo,append=TRUE,row.names=FALSE)
dbReadTable(con,"geneinfo")





#genelist_1 <- sort(unique(read.csv("data/one.csv",stringsAsFactors=FALSE)[,2]))
#write.table(genelist_1,"genelist_pluripotgenes.txt",row.names=FALSE,col.names = FALSE,quote = FALSE)
#grep("Pou",geneinfo$genesym)




###############################################################


rs<-dbSendQuery(con, "select count(*) from genecorr")
test_genecorr <- fetch(rs,n=-1)

rs<-dbSendQuery(con, "select count(*) from geneinfo")
test_geneinfo <- fetch(rs,n=-1)

rs<-dbSendQuery(con, "select count(*) from geneexp")
test_geneexp <- fetch(rs,n=-1)


test_genecorr
test_geneinfo
test_geneexp




## diagram, axis
## diagram, gauss kernel like ola
## pou5f1 need to be included in corr calc
#### check all the genes in the diff list?
## upload other corr as well
## 


library("DESeq")


cellstates <- read.csv("../../data/cell_states.txt",sep=" ",stringsAsFactors=FALSE)[,c(2,3)]

#conditions <- cellstates[,2]
#colnames(cnt_es) == conditions
#ncol(cnt_es)
#colnames(conditions)



conditions <- merge(data.frame(cell=colnames(cnt_es),stringsAsFactors = FALSE), cellstates)$state


compare2 <- function(cnt_es, seta, setb, conditions, ds1name, ds2name){
  for(x in seta)
    conditions[conditions==x] <- "SETA"
  for(x in setb)
    conditions[conditions==x] <- "SETB"
  cds <- newCountDataSet(cnt_es, conditions)
  cds <- estimateSizeFactors( cds )
  cds <- estimateDispersions( cds )
  out <- nbinomTest(cds, "SETA", "SETB")
  out <- out[order(out$padj),]
  out2 <- out[out$padj<0.05,]
  
#  ds1name <- "foo1"
#  ds2name <- "foo2"
  out2 <- cbind(ds1name, ds2name, out2[,c(1,3,4,7,8)])
  colnames(out2) <- c("dataset1","dataset2","geneid","mean1","mean2","pvalue","padj")
  
  #Store twice for simplicity?
  
  dbGetQuery(con,sprintf("delete from diffexp WHERE dataset1='%s' AND dataset2='%s';",ds1name,ds2name))
  dbGetQuery(con,sprintf("delete from diffexp WHERE dataset1='%s' AND dataset2='%s';",ds2name,ds1name))
  dbWriteTable(con,"diffexp",out2,append=TRUE,row.names=FALSE)
  
}

#dbGetQuery(con,"delete from diffexp;")

set_2i=c("2i","2C")
set_nanog=c("Nanog_hi","Nanog_med","Nanog_lo")

compare2(cnt_es, "2i",   "2C",            conditions, "es_2i_2i",  "es_2i_2C")
compare2(cnt_es, set_2i, set_nanog,       conditions, "es_2i",     "es_nanog")

compare2(cnt_es, "a2i", set_2i,           conditions, "es_a2i",    "es_2i")
compare2(cnt_es, "a2i", set_nanog,        conditions, "es_a2i",    "es_nanog")

compare2(cnt_es, "Nanog_hi", "Nanog_med", conditions, "Nanog_hi",  "Nanog_med")
compare2(cnt_es, "Nanog_hi", "Nanog_lo",  conditions, "Nanog_hi",  "Nanog_lo")
compare2(cnt_es, "Nanog_med", "Nanog_lo", conditions, "Nanog_med", "Nanog_lo")

# 2C 2i
# (2C + 2i) vs nanog*
# a2i + (2*)
# vs nanog
# nanog_hi & lov
# nanog_hi & me (all)

