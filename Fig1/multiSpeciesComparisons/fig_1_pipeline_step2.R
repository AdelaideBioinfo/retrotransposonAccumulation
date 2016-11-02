
# This script is designed to take processed repeat data from step 1 and map each species' genome to human
# For this script to work the headers from net AXT alignmnets need to be inside a file named refSepcies.queSpecie.axt.txt 
# where refSpecies needs to be the genome name for the reference species and queSpecies needs to be the genome name for the query species.
# eg. for human and chimp "hg19.panTro4.axt.txt"


rm(list = ls())

library(GenomicRanges)


setwd("~/Desktop/retrotransposonAccumulationAnalysis/retrotransposonAccumulation/")

# load output from step one
load("../accesoryFiles/R_objects/PCA_species")
load("../accesoryFiles/R_objects/Rep_info_species")
source("baseScripts/functions.R")

# output dir
R_objectPath <- "../accesoryFiles/R_objects/"



### Our functions will only work in the correct environmnet
referenceDF <- data.frame(spec = "Human", genome = "hg19")
speciesDF <- data.frame(spec = c( "Dog", "Chimp", "Mouse", "Rhesus"), genome = c("canFam3","panTro4", "mm9", "rheMac3"))

# maybe write this as a function/functions? 

for(i in 1:nrow(speciesDF)){
  
  refSpecGenome <- referenceDF$genome[1] 
  queSpecGenome <- speciesDF$genome[i]
  refSpec <- referenceDF$spec[1]
  queSpec <- speciesDF$spec[i]
  
  align <- removeRepAlign(refSpecGenome=refSpecGenome, queSpecGenome=queSpecGenome, refSpec=refSpec, queSpec=queSpec, alignPath = "../accesoryFiles/Data/usable_alignmnet/")
  # now that i can remove the repeats all i have to do is get the bin map 
  align2 <- isolateBinAlign(align=align,refSpec=refSpec,queSpec=queSpec)
  binMap <- buildBinMap(align=align2,refSpec=refSpec, queSpec=queSpec)
  
  
  dat <- list(get(paste(refSpec, "PCA", sep ="")), get(paste(queSpec, "PCA", sep ="")), binMap)
  names(dat) <- c(paste(refSpec, "Ref", sep = ""), paste(queSpec, "Que", sep = ""), "binMap")
  assign(x=paste(refSpec, "Ref_", queSpec, "Que", sep = ""), dat)
}

saveList <- c(paste(referenceDF$spec[1], "Ref_", speciesDF$spec, "Que", sep = ""))
save(list=saveList, file=paste(R_objectPath,"binMaps", sep = ""))

