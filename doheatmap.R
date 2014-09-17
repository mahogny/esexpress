source("common.R")


#tmpfile <- tempfile("tmp/","",".svg")

#gene <- "ENSMUSG00000000126"


#tmpfile <- "/tmp/foo1"

#tmpfile <- "/home/mahogny/Dropbox/ebi/olaES/site/esexpress/temp"
#tmpfile <- "|cat"


genes <- c(
  "ENSMUSG00000000126",
  "ENSMUSG00000000028"
  )
genes<-data.frame(fromgene=genes)

dbSendQuery(con, "CREATE TEMPORARY TABLE geneset (fromgene TEXT)")
dbWriteTable(con,"geneset",genes,append=TRUE,row.names=FALSE)
rs<-dbSendQuery(con, "select * from genecorr where fromgene in (select * from geneset)")
dat <- fetch(rs,n=-1)
#dat$fromgene

togenes <- expandarray(dat$togene[1])
ind <- which(togenes %in% genes$fromgene)
mat <- matrix(nrow=0,ncol=length(ind))
for(i in 1:nrow(dat)){
  mat<-rbind(mat,as.double(expandarray(dat$corr[i])[ind]))
}
dbSendQuery(con, "DROP TABLE geneset")

colnames(mat)<-togenes[ind]
rownames(mat)<-genes$fromgene





# 
# apply(togenes,1,funcion(x) which())
# 
# 
# rs<-dbSendQuery(con, sprintf("select * from genecorr where fromgene='%s'",gene))
# dat <- fetch(rs,n=-1)
# 

#?png
#postscript(file = tmpfile, command="cat",type="svg")
#postscript(file = tmpfile, command="cat")
#svg(filename = tmpfile)
#bitmap(file = "%stdout", type="png16m")


bitmap("|cat",type="pngalpha")
#bitmap("|cat",type="svg")
#postscript(file = tmpfile, command="cat")
#svg("|cat")


heatmap(mat)
#pie(rep(1, 24), col = rainbow(24))
#dev.off()
dev.off.wrap()


#print(argv) 

dbDisconnect(con)
#dbUnloadDriver(drv)
# 

