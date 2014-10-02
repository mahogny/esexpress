source("uploadsql.r")

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




## diagram, axis
## diagram, gauss kernel like ola
## pou5f1 need to be included in corr calc
#### check all the genes in the diff list?
## upload other corr as well
## 




cellstates <- read.csv("../../data/cell_states.txt",sep=" ",stringsAsFactors=FALSE)[,c(2,3)]

#conditions <- cellstates[,2]
#colnames(cnt_es) == conditions
#ncol(cnt_es)
#colnames(conditions)





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
  dbWriteTable(con,c("esexpress","diffexp"),out2,append=TRUE,row.names=FALSE)
  
}

#dbGetQuery(con,"delete from diffexp;")

conditions <- merge(data.frame(cell=colnames(cnt_es),stringsAsFactors = FALSE), cellstates)$state

set_2i=c("2i","2C")
set_nanog=c("Nanog_hi","Nanog_med","Nanog_lo")

compare2(cnt_es, "2i",   "2C",            conditions, "es_2i: blastocyst-like",  "es_2i: 2C-like")
compare2(cnt_es, set_2i, set_nanog,       conditions, "es_2i",                   "es_nanog") ????????????????

compare2(cnt_es, "a2i", set_2i,           conditions, "es_a2i",    "es_2i")
compare2(cnt_es, "a2i", set_nanog,        conditions, "es_a2i",    "es_nanog") ???????????????????

compare2(cnt_es, "Nanog_hi",  "Nanog_med", conditions, "es_serum: more pluripotent cells",  "es_serum: primed cells")
compare2(cnt_es, "Nanog_hi",  "Nanog_lo",  conditions, "es_serum: more pluripotent cells",  "es_serum: differentiating cells")
compare2(cnt_es, "Nanog_med", "Nanog_lo",  conditions, "es_serum: primed cells",            "es_serum: differentiating cells")

#rs<-dbSendQuery(con, "select distinct dataset1,dataset2 from diffexp")
#fetch(rs,n=-1)

# 
# all 2i -> 2i
# all serum -> serum
# all a2i -> a2i
# 
# 2C -> 2i: 2C-like
# 2i_non2C -> 2i: blastocyst-like
# Nanog-low -> serum: differentiating cells
# Nanog-medium -> serum: primed cells
# Nanog-high -> serum: more pluripotent cells



# 2C 2i
# (2C + 2i) vs nanog*
# a2i + (2*)
# vs nanog
# nanog_hi & lov
# nanog_hi & me (all)



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

head(gotable)

dbWriteTable(con,c("esexpress","goinfo"),gotable,append=TRUE,row.names=FALSE)





################################


uploadgodm <- function(thef, ds1, ds2){
  onedm <- read.csv(thef,sep="\t",stringsAsFactors=FALSE)
  onedm <- cbind(ds1, ds2, onedm)
  colnames(onedm) <- c("dataset1","dataset2","goid","pvalue","tscore")
  dbGetQuery(con,sprintf("delete from godm where dataset1='%s' AND dataset2='%s';",ds1,ds2))
  dbWriteTable(con,c("esexpress","godm"),onedm,append=TRUE,row.names=FALSE)
  onedm <- read.csv(thef,sep="\t",stringsAsFactors=FALSE)
  onedm <- cbind(ds2, ds1, onedm)
  colnames(onedm) <- c("dataset1","dataset2","goid","pvalue","tscore")
  onedm$tscore <- -onedm$tscore
  dbGetQuery(con,sprintf("delete from godm where dataset1='%s' AND dataset2='%s';",ds2,ds1))
  dbWriteTable(con,c("esexpress","godm"),onedm,append=TRUE,row.names=FALSE)
}

uploadgodm("../../data/dm/goPvalue_2i2_a2i2.txt",    "2i2",   "a2i2")
uploadgodm("../../data/dm/goPvalue_serum2_2i2.txt",  "serum2","2i2")
uploadgodm("../../data/dm/goPvalue_serum2_a2i2.txt", "serum2","a2i2")


dbGetQuery(con,sprintf("select * from godm limit 1","2i2"))


##################################

uploadgenedm <- function(thef,ds){
  onedm <- read.csv(thef,sep="\t",stringsAsFactors=FALSE)
  onedm <- cbind(ds, onedm)
  colnames(onedm) <- c("dataset","geneid","genedm")
  #print(head(onedm))
  dbGetQuery(con,sprintf("delete from genedm where dataset='%s';",ds))
  dbWriteTable(con,c("esexpress","genedm"),onedm,append=TRUE,row.names=FALSE)
}
dbGetQuery(con,sprintf("delete from genedm;"))
uploadgenedm("../../data/dm/DM2i2.txt",    "es_2i")
uploadgenedm("../../data/dm/DMa2i2.txt",   "es_a2i")
uploadgenedm("../../data/dm/DMserum2.txt", "es_serum")
