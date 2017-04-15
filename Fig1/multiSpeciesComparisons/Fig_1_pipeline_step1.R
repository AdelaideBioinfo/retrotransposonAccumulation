# this can be the fig 1 pipeline 



library(GenomicRanges)
library(rtracklayer)




rm(list = ls())


# this script is designed to import raw repeat data and process it into repeat groups

setwd("~/Desktop/retrotransposonAccumulationAnalysis/retrotransposonAccumulation/")

## set up output dir
R_objectPath <- "../accesoryFiles/R_objects/"

## add the path to rep_db script and functions.R script

source("baseScripts/rep_db.R")
source("baseScripts/functions.R")
options(stringsAsFactors=FALSE)


# pull each species
bin.size = 1000000
speciesDF <- data.frame(spec = c("Human", "Dog", "Chimp", "Mouse", "Rhesus", "Cow", "Pig"), genome = c("hg19", "canFam3","panTro4", "mm9", "rheMac3", "bosTau7", "susScr2"))


for(s in 1:(nrow(speciesDF))){
  rep <- rep_info(spec1=speciesDF$spec[s], genome=speciesDF$genome[s])
  # remove unplaced chromosomes from repeats
  TE.names <- names(rep)
  for(i in TE.names){
    if(length(grep("_", rep[[i]]$genoName)) > 0){
      rep[[i]] <- rep[[i]][-(grep("_", rep[[i]]$genoName)),]
    }
  }
  bins <- binned.genome.reader(genome=speciesDF$genome[s], bin.size=bin.size, keep.rate=.9)
  bins <- bins[[1]]
  # remove unplaced chromosomes from bins
  if(length(grep("_", bins$chr)) > 0){
    bins <- bins[-(grep("_", bins$chr)),]
  }
  bin.sort = binSort(repList = rep, bins=bins,TE.names = names(rep), repType = rep("repeats", length(rep)))
  assign(paste(speciesDF$spec[s], "covCount", sep = "_"), bin.sort$counts)
  assign(paste(speciesDF$spec[s], "covRate", sep = "_"), bin.sort$rates)
  assign(paste(speciesDF$spec[s], "repInfo", sep = "_"), rep)
}



for(s in 1:nrow(speciesDF)){
  bin.rate <- get(paste(speciesDF$spec[s], "covRate", sep = "_"))
  pca <- prcomp(bin.rate[,5:ncol(bin.rate)], scale.=T)
  repStruct <- repFamStruct(speciesDF$spec[s])
  agg <- aggregate(pca$rotation[repStruct$TEname,1:2], by=list(repStruct$repType), FUN = mean)
  rownames(agg) <- agg[,1]
  agg <- agg[,c("PC1", "PC2")]
  sinePC <- colnames(agg)[max(abs(agg["new_SINE",])) == abs(agg["new_SINE",])]
  ancPC <- colnames(agg)[colnames(agg)!=sinePC]

  sumPCA <- summary(pca)
  bin.ratePCA <- list(binInfo = data.frame(bin.rate[,1:4]),
                      x = data.frame(ancient_PC = pca$x[,ancPC] * (agg["ancient",ancPC] / abs(agg["ancient",ancPC])), 
                                     new_SINE_PC = pca$x[,sinePC] * (agg["new_SINE",sinePC] / abs(agg["new_SINE",sinePC]))
                      ),
                      rotation = data.frame(ancient_PC = pca$rotation[,ancPC] * (agg["ancient",ancPC] / abs(agg["ancient",ancPC])), 
                                            new_SINE_PC = pca$rotation[,sinePC] * (agg["new_SINE",sinePC] / abs(agg["new_SINE",sinePC]))
                      ),
                      importance = data.frame(ancient_PC = sumPCA$importance[2,ancPC], 
                                              new_SINE_PC = sumPCA$importance[2,sinePC]),
                      data = data.frame(bin.rate[,5:ncol(bin.rate)])
  )
  assign(paste(speciesDF$spec[s], "PCA", sep = ""), bin.ratePCA)
}

#reuben.biplot(x = PigPCA$x, y = PigPCA$rotation)

# maybe we need to put genes in there and ordinate it that way
# protein coding regions
# reuben.biplot(x=HumanPCA$x, y = HumanPCA$rotation, x.col = 8)
# reuben.biplot(x=DogPCA$x, y = DogPCA$rotation, x.col = 8)
# reuben.biplot(x=ChimpPCA$x, y = ChimpPCA$rotation, x.col = 8)
# reuben.biplot(x=MousePCA$x, y = MousePCA$rotation, x.col = 8)
# reuben.biplot(x=RhesusPCA$x, y = RhesusPCA$rotation, x.col = 8)


objectsPCA <- paste(speciesDF$spec, "PCA", sep = "")
objectsRep <- paste(speciesDF$spec, "_repInfo", sep = "")
save(list=objectsPCA, file=paste(R_objectPath, "PCA_species", sep = ""))
save(list=objectsRep, file=paste(R_objectPath , "Rep_info_species", sep = ""))









#cols <- colorRampPalette(colors = c("blue", "white", "red"))

#heatmap(cor(bin.rate[,5:ncol(bin.rate)]), zlim = c(-1,1), scale = "none", col = cols(40))




