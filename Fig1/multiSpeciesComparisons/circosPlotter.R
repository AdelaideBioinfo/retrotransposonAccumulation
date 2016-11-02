
# this script is designed to produce circos plots reflecting the 20% tails of PCA and PC2 repeat distributions
# must also have replication domain data available.
# can be downloaded as mentioned in sup material

rm(list = ls())


setwd("~/Desktop/retrotransposonAccumulationAnalysis/retrotransposonAccumulation/")

# load output from step one and two
load("../accesoryFiles/R_objects/PCA_species")
load("../accesoryFiles/R_objects/binMaps")
source("baseScripts/functions.R")

# replication domain path
# rember to keep file names the same as on GEOdatabase
RDpath = "../accesoryFiles/Data/DNN_HMM_repliDomains/"

plotPath = "../plots/Fig1/"
circosAncName = "AncDist.pdf" 
circosNewName = "newDist.pdf"
  
referenceDF <- data.frame(spec = "Human", genome = "hg19")
speciesDF <- data.frame(spec = c( "Dog", "Chimp", "Mouse", "Rhesus"), genome = c("canFam3","panTro4", "mm9", "rheMac3"))


for(s in speciesDF$spec){
  assign(s,mapRemodeler(refSpec="Human", queSpec = s, cutoff = .05))
}

RDs <- read.table(file = paste(RDpath, "GSE53984_GSM923452_Huvec_Rep1_segments.bed", sep = ""), 
                  col.names = c("chr", "start", "end", "value"), colClasses = c("character", "integer", "integer", "character"))
RDs$value[RDs$value == "ERD"] = 1
RDs$value[RDs$value == "LRD"] = -1
RDs$value[RDs$value == "UTZ"] = NA
RDs$value[RDs$value == "DTZ"] = NA
RDs <- RDs[complete.cases(RDs),]    
RDs$value <- as.numeric(RDs$value)
RDs <- RDs[RDs$end - RDs$start + 1 >= 2000000,]

library("circlize")

layout(1)
par(mar = c(1,1,1,1))



col1 = "darkorange"
col2 = "aquamarine4"
pdf(file = paste(plotPath, circosAncName, sep = ""),onefile = T)

circos.initializeWithIdeogram(species = "hg19", plotType = c( "labels"), sort.chr = T)
circos.par("track.height" = .05)
circos.genomicTrackPlotRegion(data = RDs, bg.border = "white", panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0 , "red","white"), border = NA, ...)
})

for(s in nrow(speciesDF):1){
  List <- get(as.character(speciesDF$spec[s]))
  coords = get(paste("HumanRef_", speciesDF$spec[s],"Que", sep = "" ))
  coords <- coords$HumanRef$binInfo[List$ancient_PC$remodeldPC$refNo,]
  datQue <- List$ancient_PC$remodeldPC$quePC
  bedS <- data.frame(coords[,1:3], value = datQue)
  bedS = bedS[order(bedS$value),][c(1:as.integer(.2*nrow(bedS)), (nrow(bedS)-(as.integer(.2*nrow(bedS))-1)):nrow(bedS)),]
  circos.par("track.height" = .08)
  circos.genomicTrackPlotRegion(data = bedS, panel.fun = function(region,value, ...){
    circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, col1, col2), border = NA, ...)
  })
  
}

bedH = data.frame(HumanPCA$binInfo[,1:3], value = (HumanPCA$x$ancient_PC))
bedH = bedH[order(bedH$value),][c(1:as.integer(.2*nrow(bedH)), (nrow(bedH)-(as.integer(.2*nrow(bedH))-1)):nrow(bedH)),]
circos.par("track.height" = .08)
circos.genomicTrackPlotRegion(data = bedH, panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, col1, col2), border = NA, ...)
})
dev.off()



pdf(file = paste(plotPath, circosNewName, sep = ""),onefile = T)
circos.initializeWithIdeogram(species = "hg19", plotType = c( "labels"))
circos.par("track.height" = .05)
circos.genomicTrackPlotRegion(data = RDs, bg.border = "white", panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0 , "red","white"), border = NA, ...)
})

for(s in nrow(speciesDF):1){
  List <- get(as.character(speciesDF$spec[s]))
  coords = get(paste("HumanRef_", speciesDF$spec[s],"Que", sep = "" ))
  coords <- coords$HumanRef$binInfo[List$new_SINE_PC$remodeldPC$refNo,]
  datQue <- List$new_SINE_PC$remodeldPC$quePC
  bedS <- data.frame(coords[,1:3], value = datQue)
  bedS = bedS[order(bedS$value),][c(1:as.integer(.2*nrow(bedS)), (nrow(bedS)-(as.integer(.2*nrow(bedS))-1)):nrow(bedS)),]
  circos.par("track.height" = .08)
  circos.genomicTrackPlotRegion(data = bedS, panel.fun = function(region,value, ...){
    circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, col1, col2), border = NA, ...)
  })
  
}

bedH = data.frame(HumanPCA$binInfo[,1:3], value = (HumanPCA$x$new_SINE_PC))
bedH = bedH[order(bedH$value),][c(1:as.integer(.1*nrow(bedH)), (nrow(bedH)-(as.integer(.1*nrow(bedH))-1)):nrow(bedH)),]
circos.par("track.height" = .08)
circos.genomicTrackPlotRegion(data = bedH, panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, col1, col2), border = NA, ...)
})

plot.new()
legend("center", legend = c("> 80", "< 20"), fill = c(col1, col2), title = "percentile",bty = "n")
legend("bottom", legend = c("ERD"), fill = c("red"), title = "domain", bty = "n")


dev.off()



## pick some new colours tomorrow !!!










