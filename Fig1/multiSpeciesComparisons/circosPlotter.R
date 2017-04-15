
# this script is designed to produce circos plots reflecting the 20% tails of PCA and PC2 repeat distributions
# must also have replication domain data available.
# can be downloaded as mentioned in sup material

rm(list = ls())


setwd("~/Desktop/retrotransposonAccumulationAnalysis/retrotransposonAccumulation/")

# load output from step one and two
load("../accesoryFiles/R_objects/PCA_species")
load("../accesoryFiles/R_objects/binMaps")
source("baseScripts/functions.R")
source("baseScripts/rep_db.R")

# replication domain path
# rember to keep file names the same as on GEOdatabase
RDpath = "../accesoryFiles/Data/DNN_HMM_repliDomains/"

plotPath = "../plots/Fig1/"
circosAncName = "AncDist.pdf" 
circosNewName = "newDist.pdf"
  
referenceDF <- data.frame(spec = "Human", genome = "hg19")
speciesDF <- data.frame(spec = c("Cow","Pig", "Dog", "Chimp", "Mouse", "Rhesus"), genome = c("bosTau7", "susScr2","canFam3","panTro4", "mm9", "rheMac3"))


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
## ERD and LRDs then print

#### use some of the heatmap code and we can forget about over using our shitty ordination system

# we need to think about how to order our species 
# enriched to non-enriched


refSpec = "Human"
for(species in 1:6){
  queSpec = c("Chimp", "Rhesus", "Mouse", "Dog", "Cow", "Pig")[species]
  DataSet <- get(paste(refSpec, "Ref_", queSpec,"Que", sep = ""))

  minFrac = .1

rowChoice <- (1:nrow(DataSet$binMap))[DataSet$binMap$refFrac >= minFrac & DataSet$binMap$queFrac >= minFrac]
aggQue <- aggregate(DataSet[[paste(queSpec, "Que", sep = "")]]$data[DataSet$binMap$queNo,][rowChoice,] * DataSet$binMap$queFrac[rowChoice], 
                    by = list(DataSet$binMap$refNo[rowChoice]), FUN = sum)

aggHuman <- aggregate(DataSet$binMap$refFrac[rowChoice] ,
                      by = list(DataSet$binMap$refNo[rowChoice]), FUN = sum)

specReMap <- aggQue[,2:ncol(aggQue)]/aggHuman$x

# then give each remapped datset the proper name
repStruct <- repFamStruct(queSpec)
df <- DataSet$HumanRef$binInfo[aggQue$Group.1,c("chr", "start", "end")]
for(r in unique(repStruct$repType)){
  df <- cbind(df, rowSums(data.frame(specReMap[,repStruct$TEname[repStruct$repType == r]])))
  }

colnames(df)[4:ncol(df)] <- unique(repStruct$repType)

assign(paste(queSpec,"_humanised", sep = ""), value = df)

}



repStruct <- repFamStruct("Human")
df <- HumanPCA$binInfo[,c("chr", "start", "end")]
for(r in unique(repStruct$repType)){
  df <- cbind(df, rowSums(data.frame(HumanPCA$data[,repStruct$TEname[repStruct$repType == r]])))
}
colnames(df)[4:ncol(df)] <- unique(repStruct$repType)
assign(paste("Human","_humanised", sep = ""), value = df)


RDs <- read.table(file = paste(RDpath, "GSE53984_GSM923452_Huvec_Rep1_segments.bed", sep = ""), 
                  col.names = c("chr", "start", "end", "value"), colClasses = c("character", "integer", "integer", "character"))
RDs$value[RDs$value == "ERD"] = 1
RDs$value[RDs$value == "LRD"] = -1
RDs$value[RDs$value == "UTZ"] = NA
RDs$value[RDs$value == "DTZ"] = NA
RDs <- RDs[complete.cases(RDs),]    
RDs$value <- as.numeric(RDs$value)
RDs <- RDs[RDs$end - RDs$start + 1 >= 2000000,]



specList <- c("Human", "Chimp", "Rhesus", "Mouse", "Dog", "Pig", "Cow")


pdf(file = paste(plotPath, circosNewName, sep = ""),onefile = T)

circos.initializeWithIdeogram(species = "hg19", plotType = c( "labels"))
circos.par("track.height" = .05)
circos.genomicTrackPlotRegion(data = RDs, bg.border = "white", panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0 , "darkorange","white"), border = NA, ...)
})


for(s in length(specList):1){
  dat <- get(paste(specList[s], "humanised", sep = "_"))
  
  bed_new_SINE <- dat[ dat$new_SINE > sort(dat$new_SINE)[nrow(dat) * .85] ,c("chr", "start", "end")]
  bed_new_SINE$value = 1
  bed_new_LINE <- dat[ dat$new_LINE > sort(dat$new_LINE)[nrow(dat) * .85] ,c("chr", "start", "end")]
  bed_new_LINE$value = -1
  
  bed_new <- rbind(bed_new_SINE, bed_new_LINE)
  
  circos.par("track.height" = .08)
  circos.genomicTrackPlotRegion(data = bed_new, panel.fun = function(region,value, ...){
    circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, "aquamarine2", "purple"), border = NA, ...)
  })
  
  
}

dev.off()

pdf(file = paste(plotPath, circosAncName, sep = ""),onefile = T)

circos.initializeWithIdeogram(species = "hg19", plotType = c( "labels"))
circos.par("track.height" = .05)
circos.genomicTrackPlotRegion(data = RDs, bg.border = "white", panel.fun = function(region,value, ...){
  circos.genomicRect(region,value,col = ifelse(value[[1]] > 0 , "darkorange","white"), border = NA, ...)
})


for(s in length(specList):1){
  dat <- get(paste(specList[s], "humanised", sep = "_"))
  
  bed_ancient <- dat[ dat$ancient > sort(dat$ancient)[nrow(dat) * .85] ,c("chr", "start", "end")]
  bed_ancient$value = 1
  bed_old_LINE <- dat[ dat$old_LINE > sort(dat$old_LINE)[nrow(dat) * .85] ,c("chr", "start", "end")]
  bed_old_LINE$value = -1
  
  bed_old <- rbind(bed_ancient, bed_old_LINE)
  
  circos.par("track.height" = .08)
  circos.genomicTrackPlotRegion(data = bed_old, panel.fun = function(region,value, ...){
    circos.genomicRect(region,value,col = ifelse(value[[1]] > 0, "darkblue", "red"), border = NA, ...)
  })
  
}

dev.off()
