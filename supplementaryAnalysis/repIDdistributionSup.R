### sp the next step is to get all the repeat stat figs sorted



rm(list = ls())
setwd("~/Desktop/retrotransposonAccumulationAnalysis/retrotransposonAccumulation/")
library(GenomicRanges)
library(rtracklayer)


source("baseScripts/functions.R")
source("baseScripts/rep_db.R")

spec1 = c("Human", "Chimp", "Rhesus", "Mouse", "Dog", "Pig", "Cow")
genome = c("hg19", "panTro4", "rheMac3", "mm9", "canFam3", "susScr2", "bosTau7")





#pdf(file = "../plots/supFigs/mismatch/repDist.pdf", onefile = T, height = 10, width = 12.5)

#out <- c(1,3,2,4)
#lay = matrix(c(out, out + 4, out +8, out +12, out +16), nrow = 4, byrow = F)
#layout(lay)

for(s in 1:length(spec1)){
  rep <- rep_info(spec1 = spec1[s], genome = genome[s])
  repGroup <- repFamStruct(spec1 = spec1[s])
  repGroupName<- unique(repGroup$repType)
  
  if(length(unique(repGroup$repType)) == 4){
    pdf(file = paste("../plots/supFigs/mismatch/",spec1[s],"RepDist.pdf", sep = ""), onefile = T, height = 10, width = 6)
    layout(matrix(1:6, ncol = 2, byrow = TRUE))
    par(oma = c(0,0,5,0), c(7,5,0,5))
  } else if(length(unique(repGroup$repType)) == 6){
    pdf(file = paste("../plots/supFigs/mismatch/",spec1[s],"RepDist.pdf", sep = ""), onefile = T, height = 10, width = 6)
    layout(matrix(1:6, ncol = 2, byrow = TRUE))
    par(oma = c(0,0,5,0), c(7,5,0,5))
  }
  
  for(i in 1:length(repGroupName)){
    teFams <- as.character(repGroup$TEname[repGroup$repType == repGroupName[i]])
    plot(c(0,50), c(0,1), type = "n", main = repGroupName[i], xlab = "percentage mismatch", ylab = "cumulative distribution")
    grid()
    legend("bottomright", legend = teFams, fill = 1:length(teFams), cex = .7)
    for(te in 1:length(teFams)){
      print(teFams[te])
      te.dens <- rep[[teFams[te]]][,"milliDiv"] / 10
      lines(ecdf(te.dens), col = te)
    }
  }
  title(main = spec1[s], outer = TRUE, cex.main = 3)
  dev.off()
}



# should we look at element lengths ?
# 

# 



