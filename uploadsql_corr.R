source("uploadsql.R")

toupload <- Sys.getenv("upload")


if(toupload==1){
  uploadcorr("es_serum",getexpressed(ds_ola_lif))  
}
if(toupload==2){
  uploadcorr("es_2i", getexpressed(ds_ola_2i))   
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



if(toupload==7){
  uploadcorr("mES_2i: blastocyst-like",getexpressed(ds_ola_2i_2i))
}
if(toupload==8){
  uploadcorr("mES_2i: 2C-like",  getexpressed(ds_ola_2i_2C))
}
if(toupload==9){
  uploadcorr("mES_serum: more pluripotent cells", getexpressed(ds_ola_nanog_hi))
}
if(toupload==10){
  uploadcorr("mES_serum: primed cells", getexpressed(ds_ola_nanog_med))
}
if(toupload==11){
  uploadcorr("mES_serum: differentiating cells", getexpressed(ds_ola_nanog_lo))
}




rs<-dbSendQuery(con, "select distinct dataset from genecorr")
fetch(rs,n=-1)
