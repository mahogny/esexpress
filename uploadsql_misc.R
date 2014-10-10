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




## diagram, axis
## diagram, gauss kernel like ola
## pou5f1 need to be included in corr calc
#### check all the genes in the diff list?
## upload other corr as well
## 





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
  out2 <- out[which(out$padj<0.05),]
  out2 <- as.data.frame(out2,stringsAsFactors = FALSE)
  
  ds1name <- rep(ds1name,times = nrow(out2))
  ds2name <- rep(ds2name,times = nrow(out2))
  
  out3 <- data.frame(dataset1=ds1name, dataset2=ds2name, geneid=out2$id, 
                   mean1=out2$baseMeanA, mean2=out2$baseMeanB, pvalue=out2$pval, padj=out2$padj, stringsAsFactors = FALSE)

  out4 <- data.frame(dataset1=ds2name, dataset2=ds1name, geneid=out2$id, 
                   mean1=out2$baseMeanB, mean2=out2$baseMeanA, pvalue=out2$pval, padj=out2$padj, stringsAsFactors = FALSE)
#  colnames(out3) <- c("dataset1","dataset2","geneid","mean1","mean2","pvalue","padj")
  
  
  #print(length(which(is.na(out2$geneid))))
  
  dbGetQuery(con,sprintf("delete from diffexp WHERE dataset1='%s' AND dataset2='%s';",ds1name,ds2name))
  dbGetQuery(con,sprintf("delete from diffexp WHERE dataset1='%s' AND dataset2='%s';",ds2name,ds1name))
  
  #Store twice for simplicity
  print("writing table 1...")
  dbWriteTable(con,c("esexpress","diffexp"),out3,append=TRUE,row.names=FALSE)
  print("writing table 2...")
  dbWriteTable(con,c("esexpress","diffexp"),out4,append=TRUE,row.names=FALSE)

  #an alternative is to just query twice...
}

#dbGetQuery(con,"delete from diffexp;")



dbGetQuery(con,sprintf("delete from diffexp;"))

compare2(cnt_es_notnorm, "2i",   "2C",            conditions, "mES_2i: blastocyst-like",  "mES_2i: 2C-like")   #done
compare2(cnt_es_notnorm, set_2i, set_nanog,       conditions, "mES_2i",                   "mES_nanog")  # ????????????????

compare2(cnt_es_notnorm, "a2i", set_2i,           conditions, "mES_a2i",    "mES_2i")
compare2(cnt_es_notnorm, "a2i", set_nanog,        conditions, "mES_a2i",    "mES_nanog") #  ???????????????????

compare2(cnt_es_notnorm, "Nanog_hi",  "Nanog_med", conditions, "mES_serum: more pluripotent cells",  "mES_serum: primed cells")
compare2(cnt_es_notnorm, "Nanog_hi",  "Nanog_lo",  conditions, "mES_serum: more pluripotent cells",  "mES_serum: differentiating cells")
compare2(cnt_es_notnorm, "Nanog_med", "Nanog_lo",  conditions, "mES_serum: primed cells",            "mES_serum: differentiating cells")

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

# what about sandberg?

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

uploadgodm("../../data/dm/goPvalue_2i2_a2i2.txt",    "2i2",   "a2i2")   #### renaming??? ask ola
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
uploadgenedm("../../data/dm/DM2i2.txt",    "mES_2i")
uploadgenedm("../../data/dm/DMa2i2.txt",   "mES_a2i")      #ok now?
uploadgenedm("../../data/dm/DMserum2.txt", "mES_serum")
