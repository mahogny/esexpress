source("uploadsql.r")

toupload <- Sys.getenv("upload")


if(toupload==1){
  uploadcorr("es_serum",getexpressed(ds_ola_lif))  
}
if(toupload==2){
  uploadcorr("es_2i", getexpressed(ds_ola_2i))   #did it die here already?
}
if(toupload==3){
  uploadcorr("es_a2i",getexpressed(ds_ola_a2i))
}

if(toupload==4){
  uploadcorr("sandberg_earlyblast",getexpressed(ds_s_eblast))
}
if(toupload==5){
  uploadcorr("sandberg_midblast",  getexpressed(ds_s_mblast))
}
if(toupload==6){
  uploadcorr("sandberg_lateblast", getexpressed(ds_s_lblast))
}


rs<-dbSendQuery(con, "select distinct dataset from genecorr")
fetch(rs,n=-1)
