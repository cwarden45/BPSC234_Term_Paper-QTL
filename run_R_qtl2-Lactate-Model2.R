#cross = "A_BYxRM"
#allele1 = "BY"
#allele2 = "RM"
#genotypesIN = "data/genotype_A.tsv"
#phenotypesIN = "data/phenotypes.tsv"
##defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
#yaml = paste(cross,"_Lactate.Model2.yaml",sep="")
##at least for now, skip defining 'phenocovar'
#output.file1 = paste(cross,"_Lactate.Model2_LODperm.txt",sep="")
#output.file2 = paste(cross,"_Lactate.Model2_LODpeaks.txt",sep="")
#output.file3 = paste(cross,"_Lactate.Model2_LODeffect.txt",sep="")
#output.plot = paste(cross,"_Lactate.Model2_LODall.png",sep="")

cross = "375_M22xBY"
allele1 = "M22"
allele2 = "BY"
genotypesIN = "data/genotype_375.tsv"
phenotypesIN = "data/phenotypes.tsv"
#defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
yaml = paste(cross,"_Lactate.Model2.yaml",sep="")
#at least for now, skip defining 'phenocovar'
output.file = paste(cross,"_Lactate.Model2_LODpeaks.txt",sep="")
output.file1 = paste(cross,"_Lactate.Model2_LODperm.txt",sep="")
output.file2 = paste(cross,"_Lactate.Model2_LODpeaks.txt",sep="")
output.file3 = paste(cross,"_Lactate.Model2_LODeffect.txt",sep="")
output.plot = paste(cross,"_Lactate.Model2_LODall.png",sep="")

library(qtl2)

phenoOUT = paste(cross,"_pheno_ALL.csv",sep="")

##reformat `pheno` files
JB_geno.table = read.table(genotypesIN, head=T, sep="\t")
JB_pheno.table = read.table(phenotypesIN, head=T, sep="\t")
cross_samples = JB_geno.table$X

QTL2_pheno.table = JB_pheno.table[match(cross_samples, JB_pheno.table$id),]
write.csv(QTL2_pheno.table, phenoOUT, row.names=F, quote=F)

rm(JB_geno.table)
rm(JB_pheno.table)
rm(QTL2_pheno.table)
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
#have to determine that "Lactate..1" is lodcolumn=16
plot(out, map, lodcolumn=16, ylim=c(0, ymx*1.02), las=3, xlab="")
abline(v=marker_plotX, col="orange", lty=3, lwd=3)
legend("top",legend = c("chrXIV_467219_A_G (MKT1)"), col="orange", lty=3, lwd=2,
		xpd=T, inset = -0.1)
dev.off()

operm = scan1perm(genoprobs = pr, pheno = cross_obj$pheno, n_perm = 200)#due to increased run-time, change to 200 instead of 1000
print(summary(operm))
write.table(data.frame(operm), output.file1, quote=F, sep="\t", row.names=F)

peaks = find_peaks(out, map, threshold=4, peakdrop=1.8, drop=1.5)
#bayes = bayes_int(out, map, lodcolumn=1, prob=0.95)
write.table(peaks, output.file2, quote=F, sep="\t", row.names=F)

eff = scan1coef(pr[,"chrXIV"], cross_obj$pheno[,"Lactate..1"])#export effects on chromosome of interest
write.table(data.frame(eff), output.file3, quote=F, sep="\t", row.names=F)