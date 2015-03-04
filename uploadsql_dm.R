source("uploadsql.R")


################################


uploadgodm <- function(thef, ds1, ds2){
  onedm <- read.csv(thef,sep="\t",stringsAsFactors=FALSE)
  onedm <- cbind(ds1, ds2, onedm)
#  print(head(onedm))
  colnames(onedm) <- c("dataset1","dataset2","goid","pvalue","tscore")
  dbGetQuery(con,sprintf("delete from godm where dataset1='%s' AND dataset2='%s';",ds1,ds2))
  dbWriteTable(con,c("espresso","godm"),onedm,append=TRUE,row.names=FALSE)
  onedm <- read.csv(thef,sep="\t",stringsAsFactors=FALSE)
  onedm <- cbind(ds2, ds1, onedm)
  colnames(onedm) <- c("dataset1","dataset2","goid","pvalue","tscore")
  onedm$tscore <- -onedm$tscore
  dbGetQuery(con,sprintf("delete from godm where dataset1='%s' AND dataset2='%s';",ds2,ds1))
  dbWriteTable(con,c("espresso","godm"),onedm,append=TRUE,row.names=FALSE)
}

uploadgodm("../../data/dm/goPvalue_2i_a2i_mincount10_rpm.txt",    "mES_2i (2)",   "mES_a2i (2)")   
uploadgodm("../../data/dm/goPvalue_serum_2iAll_mincount10_rpm.txt",  "mES_serum (2)","mES_2i (2)")
uploadgodm("../../data/dm/goPvalue_serum_a2i_mincount10_rpm.txt", "mES_serum (2)","mES_a2i (2)")


#dbGetQuery(con,sprintf("select * from godm limit 1","2i2"))


##################################

uploadgenedm <- function(thef,ds){
  onedm <- read.csv(thef,sep="\t",stringsAsFactors=FALSE)
  onedm <- cbind(ds, onedm)
  colnames(onedm) <- c("dataset","geneid","genedm")
  #print(head(onedm))
  dbGetQuery(con,sprintf("delete from genedm where dataset='%s';",ds))
  dbWriteTable(con,c("espresso","genedm"),onedm,append=TRUE,row.names=FALSE)
}
dbGetQuery(con,sprintf("delete from genedm;"))
uploadgenedm("../../data/dm/DM2i_mincount10_rpm.txt",    "mES_2i")
uploadgenedm("../../data/dm/DMa2i_mincount10_rpm.txt",   "mES_a2i")      
uploadgenedm("../../data/dm/DMserum_mincount10_rpm.txt", "mES_serum")
