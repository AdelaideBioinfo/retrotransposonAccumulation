# this script works the same as the our other one does, it just changes the input into the UCSC type

rm(list = ls())

rep_info <- function(spec1, genome){

rep_name <- paste("~/Desktop/retrotransposonAccumulationAnalysis/accesoryFiles/Data/rmskTable/", genome,".fa.out" , sep = "")
rep <- read.table(file = rep_name, header = FALSE, skip = 3, 
                  col.names = c("swScore", "milliDiv", "milliIns", "milliDel", "genoName", "genoStart", "genoEnd", "genoLeft", 
                                "strand", "repName", "repClass", "repStart", "repEnd", "repLeft", "id"),
                  colClasses = c("integer", "double", "double", "double", "character", "integer", "integer","character", 
                                 "character", "character", "character", "character", "integer", "character", "character"),
                  fill = TRUE
)
rep <- rep[rep$id != "",]
rep[rep$strand == "C",c("repStart", "repLeft")] <- rep[rep$strand == "C",c("repLeft", "repStart")]
rep[rep$strand == "C", "strand"] = "-"
rep$genoLeft <- substring(text = rep$genoLeft, 2, nchar(rep$genoLeft) - 1)
rep$repLeft <- substring(text = rep$repLeft, 2, nchar(rep$repLeft) - 1)



# we should make width proportional to the total number of repeats in each family
rep$genoLeft = as.integer(rep$genoLeft)
rep$repStart = as.integer(rep$repStart)
rep$repLeft = as.integer(rep$repLeft)

rep2 <- rep


rep$repFamily <- unlist(lapply(lapply(strsplit(rep$repClass, "/"), rev), getElement, 1))
rep$repClass <- unlist(lapply(strsplit(rep$repClass, "/"),getElement, 1))

rep$bin <- 0
rep[,c("milliDiv", "milliIns", "milliDel")] <- rep[,c("milliDiv", "milliIns", "milliDel")] * 10

rep <- rep[,c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
              "genoName", "genoStart", "genoEnd", "genoLeft",
              "strand", "repName", "repClass", "repFamily", 
              "repStart", "repEnd", "repLeft", "id")]


head(rep)




MIR <- rep[grep("MIR", rep$repFamily),]
L2 <- rep[grep("L2", rep$repFamily),]
L1 <- rep[grep("L1", rep$repFamily),]
L1ME <- L1[grep("ME", L1$repName),]
L1MD <- L1[grep("MD", L1$repName),]
L1MC <- L1[grep("MC", L1$repName),]
L1MB <- L1[grep("MB", L1$repName),]
L1MA <- L1[grep("MA", L1$repName),]

pie(table(L1$repName))

boxplot(L1$milliDiv ~ L1$repName, las = 2, width = table(L1$repName))

if(spec1 == "Pig"){
  
  L1_SS <- rbind(L1[grep("L1_SS", L1$repName),], L1[grep("HAL1_SS", L1$repName),], L1[grep("L1B_SS", L1$repName),])
  
  L1_SS_1 <- L1[grep("L1-1_SS", L1$repName),]
  L1_SS_2 <- L1[grep("L1-2_SS", L1$repName),]
  L1_SS_3 <- L1[grep("L1-3_SS", L1$repName),]
  
  SINE <- rep[rep$repClass == "SINE",]
  
  Pre0_SS <- SINE[grep("Pre0_SS", SINE$repName),]
  PRE1 <- SINE[grep("PRE1", SINE$repName),]
  SINE1_SS <- SINE[grep("SINE1", SINE$repName),]
}

if(spec1 == "Horse"){
  L1_EC <- L1[grep("_EC", L1$repName),]
  
  L1_EC_1 <- L1_EC[grep("L1-1", L1_EC$repName),]
  L1_EC_2 <- L1_EC[grep("L1-2", L1_EC$repName),]
  L1_EC_3 <- L1_EC[grep("L1-3", L1_EC$repName),]
  L1_EC_4 <- L1_EC[grep("L1-4", L1_EC$repName),]
  L1_EC_5 <- L1_EC[grep("L1-5", L1_EC$repName),]
  
  SINE <- rep[rep$repClass == "SINE",]
  ERE <- rep[grep("ERE", rep$repName),]
  
  ERE1 <- ERE[grep("ERE1", ERE$repName),]
  ERE2 <- ERE[grep("ERE2", ERE$repName),]
  ERE3 <- ERE[grep("ERE3", ERE$repName),]
  ERE4 <- ERE[grep("ERE4", ERE$repName),]
  
}

par(mar = c(10,5,5,5))
boxplot(L1_EC$milliDiv ~ L1_EC$repName, las = 2, width = table(L1_EC$repName))





if(spec1 == "Elephant"){
  L1_LA <- L1[grep("LA", L1$repName),]
  L1_LA_1 <- L1_LA[grep("L1-1", L1_LA$repName),]
  L1_LA_2 <- L1_LA[grep("L1-2", L1_LA$repName),]
  L1_LA_3 <- L1_LA[grep("L1-3", L1_LA$repName),]
  L1_LA_4 <- L1_LA[grep("L1-4", L1_LA$repName),]
  L1_LA_5 <- L1_LA[grep("L1-5", L1_LA$repName),]
  L1_LA_6 <- L1_LA[grep("L1-6", L1_LA$repName),]
  
  
  tRNA_RTE <- rep[rep$repFamily == "tRNA-RTE",]
  AFRO_LA <- tRNA_RTE[tRNA_RTE$repName == "AFRO_LA",]
  AFROSINE1 <- rbind(tRNA_RTE[tRNA_RTE$repName == "AFROSINE",], tRNA_RTE[tRNA_RTE$repName == "AFROSINE1B",])
  AFROSINE2 <- tRNA_RTE[tRNA_RTE$repName == "AFROSINE2",]
  AFROSINE3 <- tRNA_RTE[tRNA_RTE$repName == "AFROSINE3",]
  
  RTE1 <- rep[rep$repFamily == "Core-RTE",]
  
  a.i <- IRanges(start = RTE1$repStart, end = RTE1$repEnd)
  lines(coverage(a.i[RTE1$strand == "-" & RTE1$repName == "RTE1-N1b_LA"]))
  lines(coverage(a.i[RTE1$repName == "RTE1-N1b_LA"]))
  
  plot(coverage(a.i))
  
  
  b.i <- IRanges(start = RTEbovB$repStart, end = RTEbovB$repEnd)
  plot(coverage(b.i))
  
  c.i <- IRanges(start = RTE_X$repStart, end = RTE_X$repEnd)
  plot(coverage(c.i))
  
  d.i <- IRanges(start = AFRO_LA$repStart, end = AFRO_LA$repEnd)
  plot(coverage(d.i))
  
  e.i <- IRanges(start = AFROSINE3$repStart, end = AFROSINE3$repEnd)
  plot(coverage(e.i))
  
  
  RTE_BovB <- rep[rep$repFamily == "RTE-BovB",]
  RTE_X <- rep[rep$repFamily == "RTE-X",]
  
}




if(genome == "oryCun2"){
  
  CSINE1 = rep$repFamily[rep$repName == "CSINE1"]
  CSINE2 = rep$repFamily[rep$repName == "CSINE2"]
  CSINE2B = rep$repFamily[rep$repName == "CSINE2B"] 
  CSINE3A = rep$repFamily[rep$repName == "CSINE3A"] 

  L1A_OC = rep$repFamily[grep("L1A_O", rep$repName)] 
  L1A2_OC = rep$repFamily[grep("L1A2_O", rep$repName)] 
  L1C_OC = rep$repFamily[grep("L1C_O", rep$repName)] 
  L1MA = rep$repFamily[grep("L1MA", rep$repName)] 

}








# maybe lets just go with horse cause elphant is too shit
# the element coverage suggests there is a lot of misannotation




pie(table(Core_RTE$repName))


boxplot(L1$perDiv ~ L1$repName, las = 2)





head(tRNA_RTE[tRNA_RTE$strand == "+" & tRNA_RTE$genoEnd - tRNA_RTE$genoStart > 200 & tRNA_RTE$repName == "AFROSINE",])







rep$repFamily[rep$repClass =="SINE/MIR"] <- "MIR"
rep$repFamily[rep$repClass =="LINE/L2" ] <- "L2"
rep$repGroup[rep$repFamily == "L2" | rep$repFamily == "MIR"] <- "ancient"


rep$repFamily[grep("L1ME", rep$repName)] <- "L1ME"
rep$repFamily[grep("L1MD", rep$repName)] <- "L1MD"
rep$repFamily[grep("L1MC", rep$repName)] <- "L1MC"
rep$repFamily[grep("L1MB", rep$repName)] <- "L1MB"
rep$repGroup[rep$repFamily == "L1ME"



