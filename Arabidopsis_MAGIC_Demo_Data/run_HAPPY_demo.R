source("magic-mod.R") #`magic-mod.R` downloaded from http://mtweb.cs.ucl.ac.uk/mus/www/magic/magic.R
#`happy.preCC.R` was also downloaded from http://mtweb.cs.ucl.ac.uk/mus/www/magic/happy.preCC.R

##I had to modify `magic-mod.R` to use the same filename as the downloaded dependency:
#Line 1: source("happy.preCC.29062015.R") --> source("happy.preCC.R")

##`happy.hbrem` was installed from https://github.com/tavareshugo/happy.hbrem

prepare.database()

phenotypefile = "MAGIC.phenotype.example.12102015.txt" #downloaded from http://mtweb.cs.ucl.ac.uk/mus/www/magic/MAGIC.phenotype.example.12102015.txt
scan.phenotypes(phenotypefile)