## the next part is to look at replication timing 

library(rtracklayer)

rm(list = ls())

source(file="baseScripts/functions.R")


load("../accesoryFiles/R_objects/PCA_species")
mouseGR <- GRanges(seqnames = MousePCA$binInfo$chr,
                   ranges = IRanges(start = MousePCA$binInfo$start, end = MousePCA$binInfo$end),
                   score = rep(0, nrow(MousePCA$binInfo)))

fileList <- list.files(path = "../accesoryFiles/Data/mm9_RepliChip/")
sampList <-unique(gsub('.{21}$', '',  substring(fileList,21)))

dfScores <- NULL
rtList <- NULL
for(i in seq(1,length(fileList),by = 2)) {
  dat1 <- import(paste("../accesoryFiles/Data/mm9_RepliChip/",fileList[i], sep = ""), format = "BigWig")
  dat2 <- import(paste("../accesoryFiles/Data/mm9_RepliChip/",fileList[i+1], sep = ""), format = "BigWig")
  score(dat1) <- rowMeans(data.frame(score(dat1), score(dat2)))
  score(dat1) <- scale(score(dat1))
  ol <- findOverlaps(mouseGR, dat1)
  aggOl <- aggregate(score(dat1[subjectHits(ol)]), by = list(queryHits(ol)), FUN = mean)
  score(mouseGR[aggOl$Group.1]) <- aggOl$V1
  rtList <- c(rtList, list(dat1))
  dfScores = cbind(dfScores, score(mouseGR))
  score(mouseGR) <- 0
}
names(rtList) <- sampList
colnames(dfScores) <- sampList


smoothScatter(dfScores)
plot(data.frame(dfScores))

bluered <- colorRampPalette(colors = c("blue", "white", "red"))
image(cor(dfScores), col = bluered(40), zlim = c(-1,1))


plot(density(cor(dfScores)), xlim = c(0,1))



head(MousePCA)



# do the same for human and look at the pariwsie correaltions acorss bins




humanGR <- GRanges(seqnames = HumanPCA$binInfo$chr,
                   ranges = IRanges(start = HumanPCA$binInfo$start, end = HumanPCA$binInfo$end),
                   score = rep(0, nrow(HumanPCA$binInfo)))

fileList <- list.files(path = "../accesoryFiles/Data/UW_repliSeq/")
sampList <-(gsub('.{21}$', '',  substring(fileList,19)))

HumandfScores <- NULL
HumanrtList <- NULL
for(i in seq(1,length(fileList),by = 1)) {
  dat1 <- import(paste("../accesoryFiles/Data/UW_repliSeq/",fileList[i], sep = ""), format = "BigWig")
  score(dat1) <- scale(score(dat1))
  ol <- findOverlaps(humanGR, dat1)
  aggOl <- aggregate(score(dat1[subjectHits(ol)]), by = list(queryHits(ol)), FUN = mean)
  score(humanGR[aggOl$Group.1]) <- aggOl$V1
  HumanrtList <- c(HumanrtList, list(dat1))
  HumandfScores = cbind(HumandfScores, score(humanGR))
  score(humanGR) <- 0
}
names(HumanrtList) <- sampList
colnames(HumandfScores) <- sampList


hist(cor(HumandfScores), breaks = 100)

pdf(file = "../plots/supFigs/rtCon/Human.pdf")
reuben.heatmap(cor(HumandfScores), zlim = c(-1,1), scale = "none", col = bluered(20), margins = c(10,10), main = "Human")
dev.off()

pdf(file = "../plots/supFigs/rtCon/Mouse.pdf")
reuben.heatmap(cor(dfScores), zlim = c(-1,1), scale = "none", col = bluered(20), margins = c(10,10), main = "Mouse")
dev.off()




reuben.heatmap(cor(dfScores))







