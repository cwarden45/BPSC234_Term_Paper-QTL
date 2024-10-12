#setwd("E:\\UCR\\Fall2024\\BPSC_234\\Project\\Bloom_et_al_2019-Supp")

eLife_input.file = "elife-49212-fig4-data1-v2--TOP_LEFT_TABLE.txt"
RData_input.file = "NIHMS544073-supplement-01-reformat.txt"

#for my purposes, I did not add the equivalent of `pm2=c(pm, 'SWH1', 'TOR1', 'WHI2', 'ENA1', 'ENA5', 'KRE33', 'PMR1', 'MAL11')`
#...so I might be able to match the specific trait.

eLife_table = read.table(eLife_input.file, head=T, sep="\t")
RData_table = read.table(RData_input.file, head=T, sep="\t")

genes=unique(RData_table$Gene)

subset_causal = eLife_table[eLife_table$NAME %in% genes,]
print(sort(table(subset_causal$trait),decreasing = TRUE))

#top is Lactate with 6, followed by Caffeine with 4

Lactate_genes = eLife_table$NAME[eLife_table$trait == "Lactate;;1"]
Caffeine_genes = eLife_table$NAME[eLife_table$trait == "Caffeine;15mM;2"]

prev_results.Lactate_genes = RData_table[RData_table$Gene %in% Lactate_genes,]
prev_results.Caffeine_genes = RData_table[RData_table$Gene %in% Caffeine_genes,]

print(prev_results.Lactate_genes$Phenotype)
#[1] "Petitie frequency"                                                                                                                             
#[2] "Biofilm formation"                                                                                                                             
#[3] "High temperature growth"                                                                                                                       
#[4] "High temperature growth"                                                                                                                       
#[5] "High temperature growth, Sporulation efficiency, Gene expression, DNA damage, petite formation, ethanol tolerance, dipropyldopamine resistance"

print(prev_results.Caffeine_genes$Phenotype)
#[1] "High temperature growth, Sporulation efficiency, Gene expression, DNA damage, petite formation, ethanol tolerance, dipropyldopamine resistance"
#[2] "Antifungal drug resistance"  

##I am not sure about how exact the match is between "lactate" and "ethanol tolerance"
#...I believe that I have seen a lot of studies in yeast related to lactose metabolism

#So, I will try to focus on the trait for "Lactate;;1""