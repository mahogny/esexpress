source("uploadsql.R")

#dbGetQuery(con,sprintf("delete from geneexp;"))
uploadcounts("mES_serum",ds_ola_lif)
uploadcounts("mES_2i", ds_ola_2i)
uploadcounts("mES_a2i",ds_ola_a2i)
#dbReadTable(con,"geneexp")[,2]




uploadcounts("sandberg_lateblast",ds_s_lblast)
uploadcounts("sandberg_midblast", ds_s_mblast)
uploadcounts("sandberg_earlyblast",ds_s_eblast)





uploadcounts("mES_2i: blastocyst-like", ds_ola_2i_2i)
uploadcounts("mES_2i: 2C-like",         ds_ola_2i_2C)
uploadcounts("mES_serum: more pluripotent cells",ds_ola_nanog_hi)
uploadcounts("mES_serum: primed cells",          ds_ola_nanog_med)
uploadcounts("mES_serum: differentiating cells", ds_ola_nanog_lo)
