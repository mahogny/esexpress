source("uploadsql.R")

####################################################################################
## Gene info
geneinfo <- read.csv("../../data/ensembl2genename.txt",stringsAsFactors=FALSE)
geneinfo <- cbind(geneinfo,0)
colnames(geneinfo)<-c("geneid","genesym","pvalbiovar")
geneinfo<-as.data.frame(geneinfo,stringsAsFactors=FALSE)
geneinfo <- geneinfo[which(geneinfo$geneid %in% rownames(ds_ola_lif)),]


dbGetQuery(con,"delete from geneinfo;")
dbWriteTable(con,c("esexpress","geneinfo"),geneinfo,append=TRUE,row.names=FALSE)
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






#######################################################################################
readgofile <- function(){
  totgoid <- c()
  totgoname <- c()
  totgodef <- c()
  tf <- file("../../data/go-basic.obo",open="r")
  while (length(oneLine <- readLines(tf, n = 1, warn = FALSE)) > 0) {
    if(oneLine=="[Term]"){
      goid   <- readLines(tf, n = 1, warn = FALSE)   #id
      goname <- readLines(tf, n = 1, warn = FALSE)   #name
      readLines(tf, n = 1, warn = FALSE)             #namespace
      godef  <- readLines(tf, n = 1, warn = FALSE)   #def
      
      goid <- str_sub(goid,5)
      goname <- str_sub(goname,7)
      godef <- str_sub(godef,7)
      godef <- str_split(godef,pattern = "\"")[[1]][1]
      
      totgoid <- c(totgoid, goid)
      totgoname <- c(totgoname, goname)
      totgodef <- c(totgodef, godef)
      if(length(totgoid)%%1000==0)
        print(length(totgoid))
    }
  }
  close(tf)
  gotable <- data.frame(goid=totgoid, goname=totgoname, godef=totgodef, stringsAsFactors = FALSE)
  return(gotable)
}
#gotable <- readgofile()
#write.csv(gotable, "../../data/gosummary.csv")

#head(gotable)

dbWriteTable(con,c("esexpress","goinfo"),gotable,append=TRUE,row.names=FALSE)

