### this file is incomplete, because it is way too slow for practical use


source("common.R")
library(gplots)

if(!exists("datasets"))
  datasets<-c("es_lif")

genes <- geneidsym$geneid


dbSendQuery(con, "CREATE TEMPORARY TABLE datasets (dataset TEXT PRIMARY KEY)")
dbWriteTable(con,"datasets",data.frame(dataset=datasets),append=TRUE,row.names=FALSE)



## Read out all count arrays
rs<-dbSendQuery(con, "select * from geneexp WHERE dataset IN (SELECT * from datasets)")
counts <- fetch(rs,n=-1)
dbSendQuery(con, "DROP TABLE geneset")
dbSendQuery(con, "DROP TABLE datasets")

## Merge them into one single count table

#All genes to be included
totgenes <- unique(counts$fromgene)

#All cells
totcells <- c()
totds <- c()
for(i in 1:nrow(counts)){
  herecells <- expandarray(counts$fromcell[i])
  totcells <- c(totcells,
                sprintf("%s:%s",counts$dataset[i],herecells))
  totds <- c(totds,
             rep(counts$dataset[i], length(herecells)))
}

totmatrix <- matrix(0, ncol=length(totgenes), nrow=length(totcells))
colnames(totmatrix)<-mapidsym(totgenes)
rownames(totmatrix)<-totcells
for(i in 1:nrow(counts)){
  genei <- which(counts$fromgene[i]==totgenes)
  celli <- which(counts$dataset[i]==totds)
  totmatrix[celli,genei] <- as.double(expandarray(counts$exp[i]))
  as.double(expandarray(counts$exp[i])
  )}


#### print out the table
for(i in 1:ncol(totmatrix)){
  print(colnames(totmatrix)[i])
  print(" ")
}
print("\n")



dbDisconnect(con)
