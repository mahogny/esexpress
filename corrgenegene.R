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
rownames(mat)<-mapidsym(togenes[ind])



#Set the diagonal to 0, or the output matrix will be saturated
for(i in 1:ncol(mat)){
  mat[i,i]<-0
}

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

bitmap("|cat",type="pngalpha",width=graphw,height=graphw,units="px")
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

dev.off.wrap()

# tempfile <- tempfile()
# svg(tempfile)
# heatmap(mat)
# dev.off.wrap()
# 
# readBin(tempfile,what="raw")
# writeBin(stdout)


dbDisconnect(con)

