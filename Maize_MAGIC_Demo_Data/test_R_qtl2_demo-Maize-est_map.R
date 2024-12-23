library(qtl2)
set.seed(0)

height_pheno = "maize_magic_pheno-Height.csv" #only keep column for height from `maize_magic_pheno.csv`
height_json = "maize_magic-HeightOnly.json" #modify to specify above file, along with other demo files

maize = read_cross2(height_json)

local_cores = 4

##Kosambi map function
resultK = est_map(maize, map_function = "kosambi", cores = local_cores)
for (chr in names(resultK)){
	if(chr == "1"){
		result.values = as.numeric(resultK[[chr]])
		output.table = data.frame(marker = names(resultK[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
	}else{
		result.values = as.numeric(resultK[[chr]])
		temp.table = data.frame(marker = names(resultK[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
		output.table = rbind(output.table, temp.table)
	}#end else
}#end for (chr in names(resultK))
}#end for (chr in names(resultK))
write.csv(output.table, "MaizeMAGIC_gmap-Rqtl2_est_map-kosambi.csv", row.names=F, quote=F)

##save again, just in case something else is needed
save.image(file="test_MAGIC-Maize-est_map-KosambiOnly.RData")