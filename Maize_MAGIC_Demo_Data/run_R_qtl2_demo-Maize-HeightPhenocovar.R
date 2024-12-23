library(qtl2)
set.seed(0)

full_json = "maize_magic-HeightPhenocovar.json" #modify to specify above file, along with other demo files
maize = read_cross2(full_json)

local_cores = 4 #instead of 24

print("### Part 1 ###")
gmap = insert_pseudomarkers(maize$gmap , step=0.2, stepwidth="max")
pmap = interp_map(gmap, maize$gmap , maize$pmap)
pr = calc_genoprob(maize, gmap, error_prob=0.002, cores=local_cores)

### Haley-Knott Regression
print("### Part 2a ###")
out_hk = scan1(pr, maize$pheno, cores=local_cores)

operm_hk = scan1perm(pr, maize$pheno, n_perm=1000, cores=local_cores)
print(summary(operm_hk))
output.file1 = "MaizeMAGIC_Demo-HeightPhenocovar-HK_LODpeaks.txt"
write.table(data.frame(operm_hk), output.file1, quote=F, sep="\t", row.names=F)

#add some plotting and output tables, based upon https://github.com/cwarden45/BPSC234_Term_Paper-QTL/blob/main/run_R_qtl2-Lactate-Model1.R
output.plot1 = "MaizeMAGIC_Demo-HeightPhenocovar-HK_LOD.png"
png(output.plot1, type="cairo")
ymx = maxlod(out_hk) # overall maximum LOD score
plot(out_hk, pmap, lodcolumn="PH", ylim=c(0, ymx*1.02), las=3, xlab="")
abline(h=summary(operm_hk)[1], col="red", lty=3, lwd=3)
dev.off()

### Mixed Linear Model (with overall kinship matrix)
print("### Part 2b ###")
k = calc_kinship(pr, cores=local_cores)
out_lmm = scan1(pr, maize$pheno, k, cores=local_cores)

#skip permutation, based upon example?

#add some plotting, based upon https://github.com/cwarden45/BPSC234_Term_Paper-QTL/blob/main/run_R_qtl2-Lactate-Model1.R
output.plot2 = "MaizeMAGIC_Demo-HeightPhenocovar-LMM_LOD.png"
png(output.plot2, type="cairo")
ymx = maxlod(out_lmm) # overall maximum LOD score
plot(out_lmm, pmap, lodcolumn="PH", ylim=c(0, ymx*1.02), las=3, xlab="")
dev.off()

### Mixed Linear Model (with "Leave One Chromosome Out" kinship matrix)
print("### Part 2c ###")
k_loco = calc_kinship(pr, "loco", cores=local_cores)
out_loco = scan1(pr, maize$pheno, k_loco, cores=local_cores)

#skip permutation, based upon example?

#add some plotting and output tables, based upon https://github.com/cwarden45/BPSC234_Term_Paper-QTL/blob/main/run_R_qtl2-Lactate-Model1.R
output.plot3 = "MaizeMAGIC_Demo-HeightPhenocovar-LMM.LOCO_LOD.png"
png(output.plot3, type="cairo")
ymx = maxlod(out_loco) # overall maximum LOD score
plot(out_loco, pmap, lodcolumn="PH", ylim=c(0, ymx*1.02), las=3, xlab="")
dev.off()

print("### Part 2d ###")
#also create combined plot, based upon https://kbroman.org/qtl2/assets/vignettes/user_guide.html#Performing_a_genome_scan_with_a_linear_mixed_model
output.plotC = "MaizeMAGIC_Demo-HeightPhenocovar-Combined_LOD.png"
png(output.plotC, type="cairo")
color = c("slateblue", "violetred", "green3")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx = max(maxlod(out_hk), maxlod(out_lmm), maxlod(out_loco))
    #change from using main=colnames(maize$pheno)
	plot(out_hk, pmap, lodcolumn="PH", col=color[1], main="PH",
              ylim=c(0, ymx*1.02))
    plot(out_lmm, pmap, lodcolumn="PH", col=color[2], add=TRUE)
    plot(out_loco, pmap, lodcolumn="PH", col=color[3], add=TRUE, lty=2)
    legend("topleft", lwd=2, col=color, c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
dev.off()

print("### Part 3 - SKIP ###")
#similar to Maize, which I believe was the demo code

print("### Part 4 ###")
fl_peak = max(out_hk, pmap, lodcolumn="PH")
fl_pr = pull_genoprobpos(pr, pmap, fl_peak$chr , fl_peak$pos)
fl_fit1 = fit1(fl_pr, maize$pheno[,"PH"])
fl_blup = fit1(fl_pr, maize$pheno[,"PH"], blup=TRUE)
##skip for reasons of time and file size
#save.image(file="MaizeMAGIC_Demo-HeightPhenocovar.RData")