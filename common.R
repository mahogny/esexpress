if(!(R.version$major=="3" & R.version$minor=="1.1")){
        cat("Wrong version of R\n")
        quit(save="no")
}


.libPaths(c(.libPaths(),"/net/isilonP/public/rw/homes/w3_rst01/R/x86_64-redhat-linux-gnu-library/3.1"))
.libPaths(c(.libPaths(),"/var/lib/wwwrun/R/x86_64-suse-linux-gnu-library/3.1"))

library("RPostgreSQL")
library("stringr")


auth <- read.table("auth.txt",stringsAsFactors = FALSE)[,1]
dbhost <- auth[1]
dbname <- auth[2]
dbuser <- auth[3]
dbpass <- auth[4]

drv <- dbDriver("PostgreSQL")
#con <- dbConnect(drv, dbname=paste("postgres://",dbhost,"/",dbname,sep=""), user=dbuser, password=dbpass)
con <- dbConnect(drv, host=dbhost, dbname=dbname, user=dbuser, password=dbpass)

dbSendQuery(con,"SET search_path = public, esexpress;")


expandarray <- function(x){
  toret<-tryCatch({
    strsplit(substr(x,2,nchar(x)-1),",")[[1]]
  }, error=function(){
    cat("error in conversion")
    cat(x)
  })
  toret
}
encodearray <- function(x)   do.call(paste,c(as.list(x),sep=","))
encodearrayS <- function(x)  encodearray(sapply(x,function(y) paste("'",y,"'",sep="")))
encodearrayRS <- function(x) encodearray(sapply(x,function(y) paste("\"",y,"\"",sep="")))



dev.off.wrap <- function(){
  dev.off()
  invisible()
}

time_getgeneinfo <- proc.time()[3]
rs<-dbSendQuery(con, "select geneid,genesym from geneinfo")
geneidsym <- fetch(rs,n=-1)
time_getgeneinfo <- proc.time()[3]-time_getgeneinfo

mapidsym <- function(n)
  merge(data.frame(geneid=n),geneidsym,all.x=TRUE)$genesym
mapsymid <- function(n)
  merge(data.frame(genesym=n),geneidsym,all.x=TRUE)$geneid

ensureensembl <- function(genes)
  geneidsym$geneid[which(geneidsym$geneid %in% genes | geneidsym$genesym %in% genes)]

