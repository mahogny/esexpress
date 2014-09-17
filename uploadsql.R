source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")


library("DESeq2")
source("common.R")







####################################################################################
## Count tables
uploadcounts <- function(dataset, cnt){
  geneexp<-cbind(dataset,row.names(cnt),encodearrayS(colnames(cnt)),apply(cnt,1,encodearray))
  colnames(geneexp)<-c("dataset","fromgene","fromcell","exp")
  for(i in 1:nrow(geneexp)){
    v<-sprintf("insert into geneexp values ('%s','%s',array[%s],array[%s])",geneexp[i,1],geneexp[i,2],geneexp[i,3],geneexp[i,4])
    dbGetQuery(con,v)
  }
}

cnt_es<-read.table("../../n_counts_mES.txt",header=TRUE)

genesub<-rownames(cnt_es)[1:200]

ids<-which(rownames(cnt_es) %in% genesub)
ds_ola_lif<-cnt_es[ids,grep("lif",colnames(cnt_es))]
ds_ola_2i <-cnt_es[ids,grep("_2i_",colnames(cnt_es))]
ds_ola_a2i<-cnt_es[ids,grep("_a2i_",colnames(cnt_es))]

dbGetQuery(con,"delete from geneexp;")
uploadcounts("es_lif",ds_ola_lif)
uploadcounts("es_2i", ds_ola_2i)
uploadcounts("es_a2i",ds_ola_a2i)
dbReadTable(con,"geneexp")[,2]


####################################################################################
## Count tables for sandberg
datsand<-read.table("../../Sandberg_data.txt",sep=" ", stringsAsFactors = FALSE)[,-(1:7)]
sfsand <- estimateSizeFactorsForMatrix(datsand)
nsand <- t(t(datsand)/sfsand)


ids<-which(rownames(nsand) %in% genesub)
ds_s_lblast<-nsand[ids,grep("lateblast",colnames(nsand))]
ds_s_mblast<-nsand[ids,grep("midblast",colnames(nsand))]
ds_s_eblast<-nsand[ids,grep("earlyblast",colnames(nsand))]

uploadcounts("sandberg_lateblast",ds_s_lblast)
uploadcounts("sandberg_midblast", ds_s_mblast)
uploadcounts("sandberg_earlyblast",ds_s_eblast)


unique(dbReadTable(con,"geneexp")[,1])


####################################################################################
## Gene-gene correlations

uploadcorr <- function(dataset, cnt){
  thec <- cor(t(cnt),method="spearman")
  for(i in 1:nrow(thec)){
    thec[which(is.na(thec))]<-0
  }
  linc<-cbind(dataset,
              row.names(thec),
              encodearrayS(colnames(thec)),
              apply(thec,1,encodearray),
              apply(cnt,1,encodearray),
              apply(cnt,1,encodearray),
              apply(cnt,1,encodearray)
  )
#  print(linc[1,])
  colnames(linc)<-c("dataset","fromgene","togene","corr","corrp","corrv","corrm")
  for(i in 1:nrow(linc)){
    a <- linc[i,]
    #a[which(is.na(a))]<-0
    v<-sprintf("insert into genecorr values ('%s','%s',array[%s],array[%s],array[%s],array[%s],array[%s])",
               a[1],a[2],a[3],a[4],a[5],a[6],a[7])
    dbGetQuery(con,v)
  }
}
dbGetQuery(con,"delete from genecorr;")


uploadcorr("es_lif",ds_ola_lif)




####################################################################################
## Gene info
geneinfo <- read.csv("../../../liora/ensembl2genename.txt",stringsAsFactors=FALSE)
geneinfo <- cbind(geneinfo,0)
colnames(geneinfo)<-c("geneid","genesym","pvalbiovar")
geneinfo<-as.data.frame(geneinfo,stringsAsFactors=FALSE)
geneinfo <- geneinfo[which(geneinfo$geneid %in% rownames(ds_ola_lif)),]

dbGetQuery(con,"delete from geneinfo;")
dbWriteTable(con,"geneinfo",geneinfo,append=TRUE,row.names=FALSE)
dbReadTable(con,"geneinfo")




####################################################################################
## 
####################################################################################
## 
####################################################################################
## 
####################################################################################
## 
####################################################################################
## 
####################################################################################
## 



## hm. no point doing like-search for geneid!
rs<-dbSendQuery(con, "select * from geneinfo where geneid::tsvector @@ 'ENSMUSG00000078480'::tsquery")
fetch(rs,n=-1)


rs<-dbSendQuery(con, "select * from geneinfo where genesym::tsvector @@ 'Rgs6'::tsquery")
fetch(rs,n=-1)
rs<-dbSendQuery(con, "select * from geneinfo where to_tsvector(genesym) @@ to_tsquery('rgs6')")
fetch(rs,n=-1)

## Submits a statement
rs <- dbSendQuery(con, "select * from geneinfo")
fetch(rs,n=-1)

## Submit and execute the query
dbGetQuery(con, "select * from R_packages")


rs<-dbSendQuery(con, "select * from arrtest")
x<-fetch(rs,n=-1)

dbWriteTable(con,"arrtest",x,append=TRUE,row.names=FALSE)

cat(c(1,2,3),sep = ",")
sprintf("{%s}",cat(as.character(c(1,2,3)),sep = ","))
str_join(c(1,2,3))
?str_join

paste("{","aoue","}",sep="")

#foo<-"{1,2,3,4}"
#as.double(strsplit(substr(foo,2,nchar(foo)-1),",")[[1]])




lapply(as.list(x)$a,expandarray)
unlist(lapply(as.list(x)$a,expandarray))  #reshape dimension here, and done. not nice though


x<-rbind(x,x)
apply(x,1,expandarray)
?lapply

## Closes the connection
#dbDisconnect(con)
#dbUnloadDriver(drv)




