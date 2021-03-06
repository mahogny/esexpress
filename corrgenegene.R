source("common.R")
library(gplots)


if(!exists("genes")){
   genes <- c(
     "ENSMUSG00000000126",
     "ENSMUSG00000000028"
     )
   
   genes <- read.table("genelist_pluripotency.txt",stringsAsFactors = FALSE)[,1]
}

genes <- ensureensembl(genes)
if(!exists("graphw"))
  graphw<-500
if(!exists("dataset"))
  dataset<-"es_lif"


dbSendQuery(con, "CREATE TEMPORARY TABLE geneset (fromgene TEXT)")
dbWriteTable(con,"geneset",data.frame(fromgene=genes, stringsAsFactors=FALSE),append=TRUE,row.names=FALSE)
rs<-dbSendQuery(con, sprintf("select * from genecorr where fromgene in (select * from geneset) AND dataset='%s'",dataset))
dat <- fetch(rs,n=-1)
dbSendQuery(con, "DROP TABLE geneset")


if(nrow(dat)==0)
  quit(save = "no")

togenes <- expandarray(dat$togene[1])
ind <- which(togenes %in% genes)
mat <- matrix(nrow=0,ncol=length(ind))   ############somewhere here it dies; likely related to dataset
for(i in 1:nrow(dat)){
  mat<-rbind(mat,as.double(expandarray(dat$corr[i])[ind]))
}
colnames(mat)<-mapidsym(togenes[ind])
rownames(mat)<-mapidsym(togenes[ind])





bitmap("|cat",type="pngalpha",width=graphw,height=graphw,units="px")

if(ncol(mat)>1){
  #Set the diagonal to 0, or the output matrix will be saturated
  for(i in 1:ncol(mat)){
    mat[i,i]<-0
  }
  
  my_palette <- colorRamps::matlab.like
  
  heatmap.2(
    density.info = "none",
    mat,
    margins=c(10,10),
    key=TRUE,
    key.title = "Correlation",
    key.xlab="",
    key.ylab="",
    trace="none",
    col=my_palette
  )
} else {
  plot(0,0, xlim=c(0,100),ylim=c(0,100))
  text(50,50,labels = sprintf("Too few computed genes to display (%d)",ncol(mat)))
}
dev.off.wrap()

# tempfile <- tempfile()
# svg(tempfile)
# heatmap(mat)
# dev.off.wrap()
# 
# readBin(tempfile,what="raw")
# writeBin(stdout)


dbDisconnect(con)

