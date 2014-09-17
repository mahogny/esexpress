source("common.R")

if(!exists("genes")){
   genes <- c(
     "ENSMUSG00000000126",
     "ENSMUSG00000000028",
     "Zfy2"
     )
   genes <- c(
     "ENSMUSG00000000126",
     "ENSMUSG00000000028"
   )
   
}



rs<-dbSendQuery(con, "select geneid,genesym from geneinfo")
geneidsym <- fetch(rs,n=-1)
mapidsym <- function(n) merge(data.frame(geneid=n),geneidsym,all.x=TRUE)$genesym

genes <- geneidsym$geneid[which(geneidsym$geneid %in% genes | geneidsym$genesym %in% genes)]

dbSendQuery(con, "CREATE TEMPORARY TABLE geneset (fromgene TEXT)")
dbWriteTable(con,"geneset",data.frame(fromgene=genes),append=TRUE,row.names=FALSE)
rs<-dbSendQuery(con, "select * from genecorr where fromgene in (select * from geneset)")
dat <- fetch(rs,n=-1)


togenes <- expandarray(dat$togene[1])
ind <- which(togenes %in% genes)
mat <- matrix(nrow=0,ncol=length(ind))
for(i in 1:nrow(dat)){
  mat<-rbind(mat,as.double(expandarray(dat$corr[i])[ind]))
}
dbSendQuery(con, "DROP TABLE geneset")
colnames(mat)<-mapidsym(togenes[ind])
rownames(mat)<-mapidsym(genes)


cat("{")
cat("\"geneset\":[")
tab<-cbind(togenes[ind],mapidsym(togenes[ind]))
fst=TRUE
for(i in 1:nrow(tab)){
  if(fst)
    fst=FALSE
  else
    cat(",")
  cat("{\"geneid\":\"")
  cat(tab[i,1])
  cat("\",\"genesym\":\"")
  cat(tab[i,2])
  cat("\"}\n")
}
cat("]")
cat("}")




dbDisconnect(con)

