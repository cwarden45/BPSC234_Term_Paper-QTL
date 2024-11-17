cross = "A_BYxRM"
allele1 = "BY"
allele2 = "RM"
#defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
yaml = paste(cross,"_Lactate.Model1.yaml",sep="")
#at least for now, skip defining 'covar'

#cross = "375_M22xBY"
#allele1 = "M22"
#allele2 = "BY"
##defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
#yaml = paste(cross,"_Lactate.Model1.yaml",sep="")
##at least for now, skip defining 'covar'

#cross = "A_BYxRM"
#allele1 = "BY"
#allele2 = "RM"
##defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
#yaml = paste(cross,"_RANDOM_TRAIT.yaml",sep="")
##at least for now, skip defining 'covar'

#cross = "375_M22xBY"
#allele1 = "M22"
#allele2 = "BY"
##defined cross type as "haploid", based upon the R/qtl2 input file type definitions: https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type
#yaml = paste(cross,"_RANDOM_TRAIT.yaml",sep="")
##at least for now, skip defining 'covar'


library(qtl2)
library(mutoss)
library(qvalue)

pheno_to_test = c("Lactate..1")
print("##Reading Data##")
cross_obj = read_cross2(yaml)

pheno = cross_obj$pheno[,pheno_to_test]

print("##Preprocessing, Similar to QTL##")
map = insert_pseudomarkers(cross_obj$gmap, step=1)
#follow code to use "error_prob=0.002" instead of "error_prob = 0.0001"
pr = calc_genoprob(cross_obj, map, error_prob=0.002)#includes parameter map_function = c("haldane", "kosambi", "c-f", "morgan")

print("##Kinship Estimation##")
kinship = calc_kinship(pr)
print("##Heritability Calculation##")
hsq = est_herit(pheno, kinship)
print(hsq[1])
