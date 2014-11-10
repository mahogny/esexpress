source("uploadsql.R")

toupload <- Sys.getenv("upload")




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
  
  con <- connectes() #this works around any time-out problems
  
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



#dbGetQuery(con,sprintf("delete from diffexp;"))

if(toupload==1){
  compare2(cnt_es_notnorm, "2i",   "2C",            conditions, "mES_2i: blastocyst-like",  "mES_2i: 2C-like")   #done
  compare2(cnt_es_notnorm, set_2i, set_nanog,       conditions, "mES_2i",                   "mES_serum") 
}
if(toupload==2){
  compare2(cnt_es_notnorm, "a2i", set_2i,           conditions, "mES_a2i",    "mES_2i")
  compare2(cnt_es_notnorm, "a2i", set_nanog,        conditions, "mES_a2i",    "mES_serum") 
}
if(toupload==3){
  compare2(cnt_es_notnorm, "Nanog_hi",  "Nanog_med", conditions, "mES_serum: more pluripotent cells",  "mES_serum: primed cells")
  compare2(cnt_es_notnorm, "Nanog_hi",  "Nanog_lo",  conditions, "mES_serum: more pluripotent cells",  "mES_serum: differentiating cells")
  compare2(cnt_es_notnorm, "Nanog_med", "Nanog_lo",  conditions, "mES_serum: primed cells",            "mES_serum: differentiating cells")
}
if(toupload==4){
  compare2(cnt_es_notnorm, "Nanog_med", "Nanog_lo",  conditions, "mES_serum: primed cells",            "mES_serum: differentiating cells")
}
if(toupload==5){
  
  #Concatenate sandberg with olas data
  cnt_es_and_sb <- cbind(datsand, cnt_es_notnorm[rownames(datsand),])
  conditions_es_and_sb <- c(conditions, colnames(datsand))
  conditions_es_and_sb[grep("lateblast",conditions_es_and_sb)] <- "sb_lateblast"
  conditions_es_and_sb[grep("midblast",conditions_es_and_sb)] <- "sb_midblast"
  conditions_es_and_sb[grep("earlyblast",conditions_es_and_sb)] <- "sb_earlyblast"
  
  
  compare2(cnt_es_and_sb, "2i", "sb_earlyblast",  conditions_es_and_sb, "mES_2i",      "sandberg_earlyblast")
  compare2(cnt_es_and_sb, "2i", "sb_midblast",    conditions_es_and_sb, "mES_2i",      "sandberg_midblast")
  compare2(cnt_es_and_sb, "2i", "sb_lateblast",   conditions_es_and_sb, "mES_2i",      "sandberg_lateblast")
}

