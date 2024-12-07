cross = "A_BYxRM"
allele1 = "BY"
allele2 = "RM"
genotypesIN = "data/genotype_A.tsv"
phenotypesIN = "data/phenotypes.tsv"
#defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
yaml = paste(cross,"_Lactate.Model1.yaml",sep="")
#at least for now, skip defining 'phenocovar'
output.table = paste(cross,"_Lactate.lsr_F.NO_PRUNE.pvalue_and_FDR.txt",sep="")
output.plot = paste(cross,"_Lactate.lsr_F.NO_PRUNE.pvalue_vs_FDR.png",sep="")

#cross = "375_M22xBY"
#allele1 = "M22"
#allele2 = "BY"
#genotypesIN = "data/genotype_375.tsv"
#phenotypesIN = "data/phenotypes.tsv"
##defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
#yaml = paste(cross,"_Lactate.Model1.yaml",sep="")
##at least for now, skip defining 'phenocovar'
#output.table = paste(cross,"_Lactate.lsr_F.NO_PRUNE.pvalue_and_FDR.txt",sep="")
#output.plot = paste(cross,"_Lactate.lsr_F.NO_PRUNE.pvalue_vs_FDR.png",sep="")

library(qtl2)
library(mutoss)
library(qvalue)

pheno_to_test = c("Lactate..1")
print("##Reading and Converting Data##")
cross_obj = read_cross2(yaml)

pheno = cross_obj$pheno[,pheno_to_test]

x_details = unlist(cross_obj$gmap)
geno.table = matrix(unlist(cross_obj$geno),
					ncol=nrow(cross_obj$geno[[1]]),
					byrow=TRUE)
rownames(geno.table) = names(x_details)
colnames(geno.table) = rownames(cross_obj$geno[[1]])

geno.table = geno.table - 1 #use 0 and 1, instead of 1 and 2

#use function from https://github.com/cwarden45/COHCAP/blob/master/R/COHCAP.site.R
##There may not be any missing imputed genotypes, but this should allow code to also work if there are only directly observed genotypes
print("##Calculating P-Values##")
LSR_Wald = function(arr, var1){
	y= var1[!is.na(arr)]
	x= arr[!is.na(arr)]
	n = length(x)
	meanX = mean(x)
	meanY = mean(y)
	
	#slope estimate and variance of slope estimate from Chapter 14 of https://www.amazon.com/Mathematical-Statistics-Analysis-Available-Enhanced/dp/0534399428
	SS_XY = sum((x - meanX)*(y - meanY))
	SS_XX = sum((x - meanX)^2)
	Bhat = SS_XY / SS_XX
	
	intercept = meanY - Bhat * meanX
	
	RSS = sum((y - intercept - Bhat * x)^2)
	sample_variance = RSS / (n-2) #for 1 variable
	Var_Bhat = sample_variance / SS_XX
	
	#Wald Test from this class (such as Chapter 8 of https://www.amazon.com/Principles-Statistical-Genomics-Shizhong-Hardcover/dp/B011MCM4RG)
	Wald_Statistic = Bhat^2 / Var_Bhat
	return(pf(Wald_Statistic, df1 = 1, df2 = n-2, lower.tail = FALSE))
}#end def fastLm_wrapperV2

lm.pvalue = apply(geno.table, 1, LSR_Wald, pheno)
print("##Calculating P-Value Adjustment##")
Bonferroni.correction = p.adjust(lm.pvalue, method="bonferroni")
Sidak_obj = sidak(lm.pvalue, alpha=0.05, silent=TRUE)
Sidak.correction = unlist(Sidak_obj$adjPValues)
BH.FDR = p.adjust(lm.pvalue, method="fdr")
q.value = qvalue(p = lm.pvalue)$qvalue

print("##Write Summary Statistic##")
summary.table = data.frame(Marker = names(x_details),
						 lm.pvalue,
						 Bonferroni.correction, Sidak.correction, 
						 BH.FDR, q.value)
write.table(summary.table, output.table, quote=F, sep="\t", row.names=F)

print("##Create Summary Plot##")
extract.chr = function(full.name){
	nameArr = unlist(strsplit(full.name, split="\\."))
	return(nameArr[1])
}#end def extract.chr
marker.chr = sapply(names(x_details), extract.chr)
X_plot = 1:length(x_details)
chr_boundary = tapply(X_plot, marker.chr, max)
chr_midpoint = tapply(X_plot, marker.chr, median)

Bonferroni.threshold = max(summary.table$lm.pvalue[summary.table$Bonferroni.correction < 0.05])
Sidak.threshold = max(summary.table$lm.pvalue[summary.table$Sidak.correction < 0.05])
BH.FDR.threshold = max(summary.table$lm.pvalue[summary.table$BH.FDR< 0.05])
q.value.threshold = max(summary.table$lm.pvalue[summary.table$q.value < 0.05])

png(output.plot, type="cairo")
plot(X_plot, -log10(lm.pvalue), type="l",
	xaxt="n", xlab="",ylab="-Log10(P-Value)")
legend("top", legend=c("Bonferroni","BH FDR < 0.05","Sidak","q-value < 0.05"),
		col=c("red","blue","purple","green"),lwd=2,
		ncol=2, xpd=T, inset = -0.15)
mtext(unique(marker.chr), side=1, at =chr_midpoint, las=2, line=1)
abline(v=chr_boundary, lty=3, col="gray")
abline(h=-log10(Bonferroni.threshold), col="red", lwd=3)
abline(h=-log10(Sidak.threshold), col="purple", lwd=2)
abline(h=-log10(BH.FDR.threshold), col="blue", lwd=2)
abline(h=-log10(q.value.threshold), col="green", lwd=2)
dev.off()