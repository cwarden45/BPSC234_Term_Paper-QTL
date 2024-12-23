library(qtl2)
set.seed(0)

height_pheno = "arabmagic_pheno-Height.csv" #only keep column for height from `arabmagic_pheno.csv`
height_json = "arabmagic_tair9-Height.json" #modify to specify above file, along with other demo files
arab = read_cross2(height_json)

local_cores = 4

#default (Haldane)
#result = est_map(arab, cores = local_cores)

#the run-time was longer than the demo QTL analysis (for height-only), so I used the following commands to be able to resume writing code to export the results
##save.image(file="test_MAGIC-est_map.RData")
load(file="test_MAGIC-est_map.RData")

for (chr in names(result)){
	if(chr == "1"){
		result.values = as.numeric(result[[chr]])
		output.table = data.frame(marker = names(result[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
	}else{
		result.values = as.numeric(result[[chr]])
		temp.table = data.frame(marker = names(result[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
		output.table = rbind(output.table, temp.table)
	}#end else
}#end for (chr in names(result))

write.csv(output.table, "arabmagic_gmap-Rqtl2_est_map-haldane.csv", row.names=F, quote=F)

##Kosambi map function
resultK = est_map(arab, map_function = "kosambi", cores = local_cores)
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
write.csv(output.table, "arabmagic_gmap-Rqtl2_est_map-kosambi.csv", row.names=F, quote=F)

##c-f map function
resultCF = est_map(arab, map_function = "c-f", cores = local_cores)
for (chr in names(resultCF)){
	if(chr == "1"){
		result.values = as.numeric(resultCF[[chr]])
		output.table = data.frame(marker = names(resultCF[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
	}else{
		result.values = as.numeric(resultCF[[chr]])
		temp.table = data.frame(marker = names(resultCF[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
		output.table = rbind(output.table, temp.table)
	}#end else
}#end for (chr in names(resultCF))
write.csv(output.table, "arabmagic_gmap-Rqtl2_est_map-CF.csv", row.names=F, quote=F)

##morgan map function
resultM = est_map(arab, map_function = "morgan", cores = local_cores)
for (chr in names(resultM)){
	if(chr == "1"){
		result.values = as.numeric(resultM[[chr]])
		output.table = data.frame(marker = names(resultM[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
	}else{
		result.values = as.numeric(resultM[[chr]])
		temp.table = data.frame(marker = names(resultM[[chr]]),
								chr = rep(chr, length(result.values)),
								pos = result.values)
		output.table = rbind(output.table, temp.table)
	}#end else
}#end for (chr in names(resultM))
write.csv(output.table, "arabmagic_gmap-Rqtl2_est_map-morgan.csv", row.names=F, quote=F)

##save again, just in case something else is needed
save.image(file="test_MAGIC-est_map-OVERALL.RData")