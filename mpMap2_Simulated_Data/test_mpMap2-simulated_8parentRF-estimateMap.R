#setwd("C:\\Users\\Charles\\Documents\\UCR\\mpMap2_Simulated_Data")

library(mpMap2)
set.seed(0)

#list of functions : https://rdrr.io/cran/mpMap2/man/

#pedigree : https://rdrr.io/cran/mpMap2/man/pedigree.html
# --> eightParentPedigreeImproperFunnels : https://rdrr.io/cran/mpMap2/man/eightParentPedigreeImproperFunnels.html
# --> eightParentPedigreeRandomFunnels : https://rdrr.io/cran/mpMap2/man/eightParentPedigreeRandomFunnels.html --> use as guide for formatting
# --> eightParentPedigreeSingleFunnel : https://rdrr.io/cran/mpMap2/man/eightParentPedigreeSingleFunnel.html

testPed1 = eightParentPedigreeRandomFunnels(initialPopulationSize = 8, 
											selfingGenerations = 2, nSeeds = 1, intercrossingGenerations = 10)
testPed2 = eightParentPedigreeRandomFunnels(initialPopulationSize = 8, 
											selfingGenerations = 2, nSeeds = 1, intercrossingGenerations = 1)
testPed3 = eightParentPedigreeRandomFunnels(initialPopulationSize = 8, 
											selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 1)
testPed4 = eightParentPedigreeRandomFunnels(initialPopulationSize = 8, 
											selfingGenerations = 10, nSeeds = 45, intercrossingGenerations = 1)
testPed5 = eightParentPedigreeRandomFunnels(initialPopulationSize = 8, 
											selfingGenerations = 8, nSeeds = 40, intercrossingGenerations = 1)

#mpcross : https://rdrr.io/cran/mpMap2/man/mpcross.html
map_truth = qtl::sim.map(len = 1800/10, n.mar = 40000/10, include.x=FALSE) #use MaizeMAGIC as a guide, then divide by 10 for testing purposes
map_truth.table = data.frame(Marker=names(map_truth[[1]]),
							Chr=rep("1",length(map_truth)),
							cM=as.numeric(map_truth[[1]]))
write.csv(map_truth.table, "testPed5_qtl_sim.map-kosambi-TRUTH.csv", quote=F, row.names=F)
sim_cross = simulateMPCross(map = map_truth, pedigree = testPed5, mapFunction = kosambi, seed = 0)

#estimateRF : https://rdrr.io/cran/mpMap2/man/estimateRF.html
mpRF = estimateRF(sim_cross)

#estimateMap : https://rdrr.io/cran/mpMap2/man/estimateMap.html
grouped = formGroups(mpRF, groups = 1)
mapResult = estimateMap(grouped,  mapFunction = rfToKosambi, maxOffset = 20)# demo uses 10, and this article uses 20 : https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.13827
mapResult.table = data.frame(Marker=names(mapResult[[1]]),
							Chr=rep("1",length(mapResult[[1]])),
							cM=as.numeric(mapResult[[1]]))
write.csv(mapResult.table, "testPed5_mpMap2_estimateMap-kosambi-TEST.csv", quote=F, row.names=F)

##there was a warning:
#In evalq((function (..., call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,  :
#  Input data had heterozygotes but was analysed assuming infinite selfing. All heterozygotes were ignored.


#for future considerations, there is an as.mpInterval() function : https://rdrr.io/cran/mpMap2/man/as.mpInterval.html