# lets put it all together 
# claculating the association with 

setwd("~/Desktop/retrotransposonAccumulationAnalysis/retrotransposonAccumulation/")

rm(list = ls())

library(GenomicRanges)
library(rtracklayer)

spec1 <- "Human"
genome = "hg19"


source(file="baseScripts/functions.R")
source(file="baseScripts/rep_db.R")

supFigPath = "../plots/supFigs/"
RDpath = "../accesoryFiles/Data/ConsTimingDomains"

# so lets read in the whole genome and find a way to build the neighbormatrix

# ref gene
web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)

web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
chrom_info <- read.delim(textConnection(txt), header = FALSE)

# pull intergenic regions
refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))

refgene_gap.gr <- filterIntergenic(refgene)

intronKeep.gr <- filterIntron(refgene)


# seq gaps
ses <- browserSession("UCSC")
genome(ses) <- genome
data.types <- c("gaps")
track.name <- c("gap")
table.name <- c("gap")
for(i in 1:length(data.types)){
  dat <- getTable(
    ucscTableQuery(
      ses, 
      track = track.name[i],
      table = table.name[i]
    )
  )
  assign(data.types[i], dat)        
}
gaps.gr <- GRanges(seqnames = Rle(gaps$chrom),
                   ranges = IRanges(start = gaps$chromStart, end = gaps$chromEnd)
)		 


# removing gaps intergenic
int <- intersect(refgene_gap.gr, gaps.gr)
OL <- as.matrix(findOverlaps(refgene_gap.gr, int))
agg <- aggregate(x=width(int)[OL[,2]], by=list(OL[,1]), FUN = sum)
bins_gene_gap <- data.frame(chr = seqnames(refgene_gap.gr), start = start(refgene_gap.gr), end = end(refgene_gap.gr), Known = width(refgene_gap.gr))
bins_gene_gap$Known[agg[,1]] = bins_gene_gap$Known[agg[,1]] - agg[,2]
bins_gene_gap <- (bins_gene_gap[bins_gene_gap$Known > bins_gene_gap$end - bins_gene_gap$start -101,])
bins_gene_gap <- (bins_gene_gap[(bins_gene_gap$end - bins_gene_gap$start) / bins_gene_gap$Known > .95,])

# remove gaps from intron
gap.int <- intersect(intronKeep.gr, gaps.gr)
gapOl <- as.matrix(findOverlaps(intronKeep.gr, gap.int))
bins_intron = as.data.frame(intronKeep.gr)[,1:4]
colnames(bins_intron)[c(1,4)] <- c("chr", "Known")
bins_intron$Knonw[gapOl[,1]] <- bins_intron$Knonw[gapOl[,1]] - width(gap.int)[gapOl[,2]]



# get repeat info
rep = rep_info(spec1=spec1,genome=genome)
# sort into intergenic bins and intronic bins
intergenic_reps <- binSort(rep=rep, bins=bins_gene_gap, TE.names=names(rep), repType = rep("repeats",length(rep)))
intron_reps <- binSort(rep=rep, bins=bins_intron, TE.names=names(rep), repType = rep("repeats",length(rep)))


# do an overlap to work out how many 

#121313687

intergenicJoinedFam <- data.frame(intergenic_reps$counts[,c("chr", "start", "end","Known")],
                                  old_L1 = rowSums(intergenic_reps$counts[,c("L1ME", "L1MD", "L1MC", "L1MB")]),
                                  new_L1 = rowSums(intergenic_reps$counts[,c("L1MA", "L1PB",  "L1PA","L1HS")]),
                                  Alu = rowSums(intergenic_reps$counts[,c("AluS", "AluY", "AluJ")]),
                                  Ancient = rowSums(intergenic_reps$counts[,c("MIR", "L2")])
)

intronJoinedFam <- data.frame(intron_reps$counts[,c("chr", "start", "end","Known")],
                              old_L1 = rowSums(intron_reps$counts[,c("L1ME", "L1MD", "L1MC","L1MB")]),
                              new_L1 = rowSums(intron_reps$counts[,c("L1MA", "L1PB", "L1PA","L1HS")]),
                              Alu = rowSums(intron_reps$counts[,c("AluS", "AluY", "AluJ")]),
                              Ancient = rowSums(intron_reps$counts[,c("MIR", "L2")])
)




joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB  ,rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))





intergenicSizesMean <- NULL
for(l in seq(0,max(intergenicJoinedFam$Known)+100,100)){
  intergenicSizesMean <- c(intergenicSizesMean, 
                           10^(mean(log10((intergenicJoinedFam$end - intergenicJoinedFam$start)[(intergenicJoinedFam$end - intergenicJoinedFam$start) > l]))))
}


intronSizesMean <- NULL
for(l in seq(0,max(intronJoinedFam$Known)+100,100)){
  intronSizesMean <- c(intronSizesMean, 10^(mean(log10((intronJoinedFam$end - intronJoinedFam$start)[(intronJoinedFam$end - intronJoinedFam$start) > l]))))
}

regions = c("intergenic", "intron")
TEcols <- c("red", "purple", "aquamarine3","darkblue")



for(i in 1:4){
  for(j in 1:2){
    pdf(paste(supFigPath,"biasLine/", regions[j],"/",names(joinRep)[i],".pdf", sep = "" ), onefile = T, height = 5, width = 5)
    
    TEfam = names(joinRep)[i]
    TEjoinFam = get(paste(regions[j], "JoinedFam", sep = ""))
    
    portion = TEjoinFam[,TEfam]/TEjoinFam$Known
    cuts <- cut(log10(TEjoinFam$Known), breaks = seq(2,8,by = .05))
    agg <- aggregate(x = portion,by = list(cuts), FUN = sum)
    # divide each number by how many times you see that cut
    summ <- table(cuts)
    
    n = sum(TEjoinFam[,TEfam])
    aggp <- aggregate(x = TEjoinFam$Known,by = list(cuts), FUN = sum)
    p = aggp$x/sum(aggp$x)
    m <- n*p
    SD <- sqrt(n*p*(1-p))
    
    plot(seq(2,8,by = .05)[summ > 0],agg$x/summ[summ > 0], main = "", xlim = c(1.9,7), 
         xlab = "interval length (kb)", ylab = "retrotransposon density", ylim = c(0,.22), pch = 16, col = TEcols[i], xaxt = "n")
    axis(side = 1, at = c(seq(0,8,by = 1)), labels = c((10^seq(0,8,by = 1))/1000))
    lines(seq(2,8,by = .05)[summ > 0],m/aggp$x, lty = 2)
    lines(seq(2,8,by = .05)[summ > 0],(m + (3*SD))/aggp$x)
    lines(seq(2,8,by = .05)[summ > 0],(m - (3*SD))/aggp$x)
    lo <- loess(agg$x/summ[summ > 0] ~ seq(2,8,by = .05)[summ > 0], span = .2)
    pred <- predict(lo,newdata = seq(2,8,by = .05),se = T)
    lines(seq(2,8,by = .05), pred$fit, col = TEcols[i], lwd = 3)
    # now we know what the preferance looks like we just have to apply it to real data. 
    yPoly = c(pred$fit + (3*pred$se.fit), rev(pred$fit - (3*pred$se.fit)))
    xPoly = c(seq(2,8,by = .05), seq(8,2,by = -.05))
    polygon(xPoly[complete.cases(yPoly)],yPoly[complete.cases(yPoly)], density = 20, col = TEcols[i], border = 0)
    
    assign(x = paste(regions[j], TEfam,  "Lo", sep = ""), value = lo)
    
    sizesMean <- get(paste(regions[j], "SizesMean", sep = ""))
    
    assign(x = paste(regions[j], TEfam,  "Pred", sep = ""), value = data.frame(position = (seq(0,max(TEjoinFam$Known)+100,100)/2),
                                                                     proportion = (predict(lo,newdata = log10(sizesMean))))
           )
    
    assign(x = paste(regions[j], TEfam,  "Bias", sep = ""), value = data.frame(position = (seq(0,max(TEjoinFam$Known)+100,100)/2),
                                                                     bias = (predict(lo,newdata = log10(sizesMean)))/(m/aggp$x)[1])
           )
    
    dev.off()
  }
}



intergenicLenChoice = max(intergenicJoinedFam$end - intergenicJoinedFam$start)/2
intronLenChoice = max(intronJoinedFam$end - intronJoinedFam$start)/2



#TEfam = "new_L1"

for(r in 1:2){
  region = c("intergenic", "intron")[r]
  
  for(i in 1:4){
    TEfam = names(joinRep)[i]
    lenChoice = get(paste(region, "LenChoice", sep =  ""))
    posStatBins = get(paste(region,"JoinedFam",sep = ""))
    
    TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,
                                            repBins = posStatBins,repList = joinRep,refgene = refgene,
                                            type = region,repType = "repeats")
    
    rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
    bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
    bpFreq[is.na(bpFreq)] <- 1
    
    rate <- rawCov/bpFreq
    
    
    
    # we can get both stats, 
    cutSite <- unique(as.integer(10^(seq(0,as.integer(log10(lenChoice+1))+1,.05))))
    cuted <- cut(1:(lenChoice+1), breaks = cutSite,right = F,ordered_result = T)
    cutPos<- as.integer(apply(data.frame(cutSite[1:(length(cutSite)-1)], cutSite[2:(length(cutSite))] ), 1, mean))
    Tab <- table(cuted)
    
    aggTEmean <- aggregate(rate, list(as.integer(cuted)), mean, simplify = F)
    aggTEsd <- aggregate(rate, list(as.integer(cuted)), sd, simplify = F)
    
    aggBPfreq <- aggregate(bpFreq, list(as.integer(cuted)), sum, simplify = F)
    aggBPraw <- aggregate(rawCov, list(as.integer(cuted)), sum, simplify = F)
    
    usePos <- cutPos[aggTEsd$Group.1]
    
    biasData <- get(paste(region, TEfam, "Bias", sep =""))
    SPbiasData <- smooth.spline((biasData$position[complete.cases(biasData)]),biasData$bias[complete.cases(biasData)],all.knots = TRUE)
    biasPred <- predict(SPbiasData, usePos)
    
    shape.x <- c(log10(usePos), rev(log10(usePos)))
    shape.y <- c(aggTEmean$x + (2* aggTEsd$x), rev(aggTEmean$x - (2* aggTEsd$x)))
    shape.y.adj <- c((aggTEmean$x + (2* aggTEsd$x))/biasPred$y, rev((aggTEmean$x - (2* aggTEsd$x))/biasPred$y))
    
    sd <- (sqrt(bpFreq * TEs_posStats$p * (1 - TEs_posStats$p))/(bpFreq))
    
    pdf(file = paste(supFigPath, "corrected/", region, "/", TEfam,".pdf", sep = ""), height = 3, width = 10)
    par(mar = c(5,5,5,2))
    plot(log10(1:(lenChoice +1)), rate, col = 3, type = "n", ylim = c(0,.25), 
         main = "", ylab = "retrotransposon\ndensity", 
         xlab = "distance from boundary (log10 bp)")
    lines(log10(usePos), aggTEmean$x, type = "l", col = "grey60", lwd = 3)
    polygon(shape.x, shape.y, density =40,border = 0,col = "grey60")
    
    lines(log10(usePos), (aggBPraw$x/aggBPfreq$x)/biasPred$y, col = TEcols[i],lwd = 3 )
    polygon(shape.x, shape.y.adj, density =40,border = 0,col = TEcols[i],angle = 135)
    
    lines(log10(1:(lenChoice+1)), TEs_posStats$p + 2*sd, col = 1, lwd = 2)
    lines(log10(1:(lenChoice+1)), TEs_posStats$p - 2*sd, col = 1, lwd = 2)  
    legend("topleft",legend = c("uncorrected", "corrected"), fill = c("grey60", TEcols[i]), bty = "n")
    # question still is, are we accuratly capturing the bias effect through our smoothing
    dev.off()
    
    assign(x = paste(TEfam,region,"adjValues", sep = "_"),value = data.frame(position = usePos, rate = aggTEmean$x, adjRate = aggTEmean$x/biasPred$y)) 
    
  }
  
}


