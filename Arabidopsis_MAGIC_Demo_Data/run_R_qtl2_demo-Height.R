library(qtl2)
set.seed(0)

height_pheno = "arabmagic_pheno-Height.csv" #only keep column for height from `arabmagic_pheno.csv`
height_json = "arabmagic_tair9-Height.json" #modify to specify above file, along with other demo files
arab = read_cross2(height_json)

local_cores = 4 #instead of 24

print("### Part 1 ###")
gmap = insert_pseudomarkers(arab$gmap , step=0.2, stepwidth="max")
pmap = interp_map(gmap, arab$gmap , arab$pmap)
pr = calc_genoprob(arab, gmap, error_prob=0.002, cores=local_cores)

### Haley-Knott Regression
print("### Part 2a ###")
out_hk = scan1(pr, arab$pheno, cores=local_cores)

operm_hk = scan1perm(pr, arab$pheno, n_perm=1000, cores=local_cores)
print(summary(operm_hk))
output.file1 = "ArabidopsisMAGIC_Demo-HeightOnly-HK_LODpeaks.txt"
write.table(data.frame(operm_hk), output.file1, quote=F, sep="\t", row.names=F)

#add some plotting and output tables, based upon https://github.com/cwarden45/BPSC234_Term_Paper-QTL/blob/main/run_R_qtl2-Lactate-Model1.R
output.plot1 = "ArabidopsisMAGIC_Demo-HeightOnly-HK_LOD.png"
png(output.plot1, type="cairo")
ymx = maxlod(out_hk) # overall maximum LOD score
plot(out_hk, pmap, lodcolumn=1, ylim=c(0, ymx*1.02), las=3, xlab="")
abline(v=summary(operm_hk)[1], col="red", lty=3, lwd=3)
dev.off()

### Mixed Linear Model (with overall kinship matrix)
print("### Part 2b ###")
k = calc_kinship(pr, cores=local_cores)
out_lmm = scan1(pr, arab$pheno, k, cores=local_cores)

#skip permutation, based upon example?

#add some plotting, based upon https://github.com/cwarden45/BPSC234_Term_Paper-QTL/blob/main/run_R_qtl2-Lactate-Model1.R
output.plot2 = "ArabidopsisMAGIC_Demo-HeightOnly-LMM_LOD.png"
png(output.plot2, type="cairo")
ymx = maxlod(out_lmm) # overall maximum LOD score
plot(out_lmm, pmap, lodcolumn=1, ylim=c(0, ymx*1.02), las=3, xlab="")
dev.off()

### Mixed Linear Model (with "Leave One Chromosome Out" kinship matrix)
print("### Part 2c ###")
k_loco = calc_kinship(pr, "loco", cores=local_cores)
out_loco = scan1(pr, arab$pheno, k_loco, cores=local_cores)

#skip permutation, based upon example?

#add some plotting and output tables, based upon https://github.com/cwarden45/BPSC234_Term_Paper-QTL/blob/main/run_R_qtl2-Lactate-Model1.R
output.plot3 = "ArabidopsisMAGIC_Demo-HeightOnly-LMM.LOCO_LOD.png"
png(output.plot3, type="cairo")
ymx = maxlod(out_loco) # overall maximum LOD score
plot(out_loco, pmap, lodcolumn=1, ylim=c(0, ymx*1.02), las=3, xlab="")
dev.off()

print("### Part 2d ###")
#also create combined plot, based upon https://kbroman.org/qtl2/assets/vignettes/user_guide.html#Performing_a_genome_scan_with_a_linear_mixed_model
output.plotC = "ArabidopsisMAGIC_Demo-HeightOnly-Combined_LOD.png"
png(output.plotC, type="cairo")
color = c("slateblue", "violetred", "green3")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx = max(maxlod(out_hk), maxlod(out_lmm), maxlod(out_loco))
    plot(out_hk, pmap, lodcolumn=1, col=color[1], main=colnames(arab$pheno),
              ylim=c(0, ymx*1.02))
    plot(out_lmm, pmap, lodcolumn=1, col=color[2], add=TRUE)
    plot(out_loco, pmap, lodcolumn=1, col=color[3], add=TRUE, lty=2)
    legend("topleft", lwd=2, col=color, c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
dev.off()

print("### Part 3 - SKIP ###")
##skip because I am not sure if `arab$fruit,` will work in this modified situation.
#snp_pr = genoprob_to_snpprob(pr, arab)
#out_snps = scan1(snp_pr, arab$fruit, cores=local_cores)

print("### Part 4 ###")
fl_peak = max(out_hk, pmap, lodcolumn="height")#use `height` instead of `fruit_length`
fl_pr = pull_genoprobpos(pr, pmap, fl_peak$chr , fl_peak$pos)
fl_fit1 = fit1(fl_pr, arab$pheno[,"height"])#use `height` instead of `fruit_length`
fl_blup = fit1(fl_pr, arab$pheno[,"height"], blup=TRUE)#use `height` instead of `fruit_length`
##skip for reasons of time and file size
#save.image(file="ArabidopsisMAGIC_Demo-HeightOnly.RData")