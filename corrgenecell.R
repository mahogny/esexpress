source("common.R")
library(gplots)
#library(MASS)

if(!exists("genes")){
  genes <- c(
    "ENSMUSG00000000126",
    "ENSMUSG00000000028"
  )
  genes <- c("Gnai3","Pbsn","H19","Apoh","Cav2","Tbx2","Zfy2","Wnt9a","Xpo6","Axin2")
  genes <- read.table("genelist_pluripotency.txt",stringsAsFactors = FALSE)[,1]
}
if(!exists("datasets"))
  datasets<-c("mES_serum")
if(!exists("graphw"))
  graphw<-500

genes <- ensureensembl(genes)


dbSendQuery(con, "CREATE TEMPORARY TABLE geneset (fromgene TEXT PRIMARY KEY)")
dbWriteTable(con,"geneset",data.frame(fromgene=genes, stringsAsFactors=FALSE),append=TRUE,row.names=FALSE)

dbSendQuery(con, "CREATE TEMPORARY TABLE datasets (dataset TEXT PRIMARY KEY)")
dbWriteTable(con,"datasets",data.frame(dataset=datasets),append=TRUE,row.names=FALSE)



## Read out all count arrays
rs<-dbSendQuery(con, "select * from geneexp WHERE fromgene IN (SELECT * FROM geneset) AND dataset IN (SELECT * from datasets)")
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



## Transform to log count
totmatrix <- log(totmatrix+1)
totmatrix <- t(totmatrix) #ola preference for graph orientation



##### todo must use genesym on the sides!!


# isoMDS(dist(totmatrix)+1)
# pc <- princomp(t(totmatrix))
# biplot(pc)
# plot(pc)



bitmap("|cat",type="pngalpha",width=graphw,height=graphw,units="px")

#cor(t(totmatrix))

#my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
my_palette <- colorRamps::matlab.like

rowcol <- rep("black",nrow(totmatrix))
rowcol[which(totds == "mES_serum")] <- "#AA0000"    #this would be better in a separate file!
rowcol[which(totds == "mES_a2i")]   <- "#EEAD0E"
rowcol[which(totds == "mES_2i")]    <- "#00008B"
rowcol[which(totds == "sandberg_earlyblast")] <- "#9ACD32"
rowcol[which(totds == "sandberg_midblast")]   <- "#789C31"
rowcol[which(totds == "sandberg_lateblast")]  <- "#556B2F"

rowcol[which(totds == "mES_2i: blastocyst-like")]            <- "#0099CC"
rowcol[which(totds == "mES_2i: 2C-like")]                    <- "#003399"
rowcol[which(totds == "mES_serum: more pluripotent cells")]  <- "#990000"
rowcol[which(totds == "mES_serum: primed cells")]            <- "#CC4D4D"
rowcol[which(totds == "mES_serum: differentiating cells")]   <- "#FF9999"


heatmap.2(
  totmatrix,
  margins=c(10,10),
  key=TRUE, 
  key.title="Expression",
  key.xlab="",
  key.ylab="",
  ColSideColors = rowcol, 
  trace="none", 
  col=my_palette,
  #dendrogram="none",   #optional
  labCol=NA,
  #labCol=NA,  #optional
  cexCol=2)

#hclust(dist(totmatrix))

dev.off.wrap()


dbDisconnect(con)
