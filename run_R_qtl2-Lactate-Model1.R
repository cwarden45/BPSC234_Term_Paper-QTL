#cross = "A_BYxRM"
#allele1 = "BY"
#allele2 = "RM"
#genotypesIN = "data/genotype_A.tsv"
#phenotypesIN = "data/phenotypes.tsv"
##defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
#yaml = paste(cross,"_Lactate.Model1.yaml",sep="")
##at least for now, skip defining 'phenocovar'
#output.file = paste(cross,"_Lactate.Model1_LODpeaks.txt",sep="")
#output.plot = paste(cross,"_Lactate.Model1_LODall.png",sep="")

cross = "375_M22xBY"
allele1 = "M22"
allele2 = "BY"
genotypesIN = "data/genotype_375.tsv"
phenotypesIN = "data/genotype_375phenotypes.tsv"
#defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
yaml = paste(cross,"_Lactate.Model1.yaml",sep="")
#at least for now, skip defining 'phenocovar'
output.file = paste(cross,"_Lactate.Model1_LODpeaks.txt",sep="")
output.plot = paste(cross,"_Lactate.Model1_LODall.png",sep="")

library(qtl2)

genoOUT = paste(cross,"_geno.csv",sep="")#simply use "A" and "B" notation, instead of varying this for each cross.
phenoOUT = paste(cross,"_pheno_Lactate.csv",sep="")
variantInfo = paste(cross,"_gmap.csv",sep="")# position can either be centiMorgans (cM) or Megabase pairs (Mbp).

##reformat files
JB_geno.table = read.table(genotypesIN, head=T, sep="\t")
JB_pheno.table = read.table(phenotypesIN, head=T, sep="\t")
cross_samples = JB_geno.table$X

pheno_to_test = c("Lactate..1")
QTL2_pheno.table = JB_pheno.table[match(cross_samples, JB_pheno.table$id),]
QTL2_pheno.table = QTL2_pheno.table[,match(c("id",pheno_to_test), names(QTL2_pheno.table))]
write.csv(QTL2_pheno.table, phenoOUT, row.names=F, quote=F)

QTL2_geno.table = JB_geno.table
colnames(QTL2_geno.table)[1]="id"
convert_geno = function(arr){
	arr[arr == 1]="A"
	arr[arr == 2]="B"
	return(arr)
}#end def convert_geno
QTL2_geno.table[,2:ncol(QTL2_geno.table)]=apply(QTL2_geno.table[,2:ncol(QTL2_geno.table)], 2, convert_geno)
print(table(QTL2_geno.table$chrXIV_467219_A_G))#print genotype frequency for variant that I am trying to understand better
write.csv(QTL2_geno.table, genoOUT, row.names=F, quote=F)

extract_chr_and_pos = function(marker){
	#print(marker)
	temp.marker_info=unlist(strsplit(marker, split="_"))
	#print(temp.marker_info)
	temp.chr = temp.marker_info[1]
	temp.pos = temp.marker_info[2]
	return(c(temp.chr, temp.pos))
}#end def extract_chr
marker = colnames(QTL2_geno.table)[2:ncol(QTL2_geno.table)]
marker_info = t(sapply(marker, extract_chr_and_pos))
colnames(marker_info) = c("chr","pos")
QTL2_map.table = data.frame(marker, marker_info)
QTL2_map.table$pos = as.numeric(as.character(QTL2_map.table$pos))
QTL2_map.table$pos = QTL2_map.table$pos / 1000000
write.csv(QTL2_map.table, variantInfo, row.names=F, quote=F)

rm(JB_geno.table)
rm(JB_pheno.table)
rm(QTL2_geno.table)
rm(QTL2_map.table)
rm(QTL2_pheno.table)
rm(marker_info)
rm(marker)
rm(cross_samples)

##run code similar to https://kbroman.org/qtl2/assets/vignettes/user_guide.html
cross_obj = read_cross2(yaml)
map = insert_pseudomarkers(cross_obj$gmap, step=1)
pr = calc_genoprob(cross_obj, map)#includes parameter map_function = c("haldane", "kosambi", "c-f", "morgan")
out = scan1(pr, cross_obj$pheno)

png(output.plot, type="cairo")
ymx = maxlod(out) # overall maximum LOD score
x_details = unlist(map)
marker_plotX = xpos_scan1(map, chr=names(map), thechr=c("chrXIV"), thepos=c(467219)/1000000)#re-coded as chrXIV.chrXIV_467219_A_G_41200
plot(out, map, lodcolumn=1, ylim=c(0, ymx*1.02), las=3, xlab="")
abline(v=marker_plotX, col="orange", lty=3, lwd=3)
legend("top",legend = c("chrXIV_467219_A_G (MKT1)"), col="orange", lty=3, lwd=2,
		xpd=T, inset = -0.1)
dev.off()

peaks = find_peaks(out, map, threshold=4, peakdrop=1.8, drop=1.5)
#bayes = bayes_int(out, map, lodcolumn=1, prob=0.95)
write.table(peaks, output.file, quote=F, sep="\t", row.names=F)
