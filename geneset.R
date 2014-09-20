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

   genes <- read.table("genelist_pluripotency.txt",stringsAsFactors = FALSE)[,1]
}

#genes <- c("Vmn1r125","Gm14036")

#Convert to ensemblid
genes <- geneidsym$geneid[which(geneidsym$geneid %in% genes | geneidsym$genesym %in% genes)]

dbSendQuery(con, "CREATE TEMPORARY TABLE geneset (fromgene TEXT);")

dbWriteTable(con,"geneset",data.frame(fromgene=genes),append=TRUE,row.names=FALSE)
rs<-dbSendQuery(con, "select * from genecorr where fromgene in (select * from geneset);")
dat <- fetch(rs,n=-1)


togenes <- expandarray(dat$togene[1])
ind <- which(togenes %in% genes)
mat <- matrix(nrow=0,ncol=length(ind))
if(nrow(dat)>0)
  for(i in 1:nrow(dat))
    mat<-rbind(mat,as.double(expandarray(dat$corr[i])[ind]))
dbSendQuery(con, "DROP TABLE geneset")
colnames(mat)<-mapidsym(togenes[ind])
rownames(mat)<-mapidsym(togenes[ind])   #was genes


#rs<-dbSendQuery(con, "select * from geneinfo NATURAL JOIN geneset;")
#dat <- fetch(rs,n=-1)

cat("{")
cat("\"geneset\":[")
tab<-cbind(togenes[ind],mapidsym(togenes[ind]))

ginfo <- data.frame(geneid=genes, genesym=mapidsym(genes),stringsAsFactors = FALSE)

### notice that gene corr is looked up, but not used here!
fst=TRUE
if(nrow(ginfo)>0)
  for(i in 1:nrow(ginfo)){
    if(fst)
      fst=FALSE
    else
      cat(",")
    cat("{\"geneid\":\"")
    cat(ginfo[i,1])
    cat("\",\"genesym\":\"")
    cat(ginfo[i,2])
    cat("\"}\n")
  }
cat("]")
cat("}")




dbDisconnect(con)

