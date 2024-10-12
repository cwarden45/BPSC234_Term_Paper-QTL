#setwd("E:\\UCR\\Fall2024\\BPSC_234\\Project\\Bloom_et_al_2019-Supp")

eLife_input.file2 = "elife-49212-fig2-data1-v2.xls"
eLife_input.file4 = "elife-49212-fig4-data1-v2--TOP_LEFT_TABLE.txt"

library(readxl)

heritability.table = read_excel(eLife_input.file2)
Lactate_table = data.frame(heritability.table[heritability.table$trait == "Lactate;;1",])
Lactate_table = Lactate_table[order(Lactate_table$within_cross_variance_explained_by_qtl..cross_validated., decreasing = TRUE),]
#consider crosses "3028", "A", and/or "B"

#"MKT1" is the "Lactate" gene related to "ethanol tolerance" from `stats_to_select_trait.R` were 
eLife_table4 = read.table(eLife_input.file4, head=T, sep="\t")
MKT1.table = eLife_table4[eLife_table4$NAME == "MKT1",]
#                    trait     ORF NAME    chr  start    end width strand
#103       Caffeine;15mM;2 YNL085W MKT1 chrXIV 466886 469837  2952      +
#111               EtOH;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#133           Glycerol;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#142          Trehalose;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#149 Magnesium_Chloride;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#150          Raffinose;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#153              YPD;37;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#156            Lactate;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#177            Sucrose;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#214             Xylose;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#216     Copper_Sulfate;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#275          Galactose;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#301            Maltose;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#304   Congo_Red;75ug/mL;3 YNL085W MKT1 chrXIV 466886 469837  2952      +
#397            Mannose;;1 YNL085W MKT1 chrXIV 466886 469837  2952      +
#    pCausalSum.Gene shared.parent   effect    localFDR jointmaxPPC
#103       0.9777199           BYa increase 0.003858364   0.8614411
#111       0.9721269           BYa decrease 0.005368463   0.8020986
#133       0.9525987           BYa decrease 0.010984339   0.7482146
#142       0.9411475           BYa decrease 0.013719745   0.7025041
#149       0.9275014           BYa decrease 0.016102423   0.7179958
#150       0.9273705           BYa decrease 0.016479270   0.6646430
#153       0.9257128           BYa decrease 0.017602359   0.6900191
#156       0.9241950           BYa decrease 0.018713374   0.7204217
#177       0.8837993           BYa decrease 0.027861707   0.6214648
#214       0.8241134           BYa decrease 0.049343658   0.5491691
#216       0.8220957           BYa decrease 0.050533755   0.5484886
#275       0.7143485           BYa decrease 0.089954786   0.3856295
#301       0.6736074           BYa decrease 0.108454242   0.3660910
#304       0.6718980           BYa decrease 0.110616313   0.3251594
#397       0.5301614           RMx decrease 0.177049574   0.2242991
#     whichjointmaxPPC
#103 chrXIV_467219_A_G
#111 chrXIV_467219_A_G
#133 chrXIV_467219_A_G
#142 chrXIV_467219_A_G
#149 chrXIV_467219_A_G
#150 chrXIV_467219_A_G
#153 chrXIV_467219_A_G
#156 chrXIV_467219_A_G
#177 chrXIV_467219_A_G
#214 chrXIV_467219_A_G
#216 chrXIV_467219_A_G
#275 chrXIV_467219_A_G
#301 chrXIV_467219_A_G
#304 chrXIV_466588_T_G
#397 chrXIV_466588_T_G

#The "Lactate;;1" MKT1 is inherited from parent "BYa"

#"BYa" is part of cross "375" (M22xBY) as well as cross "A" (BYxRM)