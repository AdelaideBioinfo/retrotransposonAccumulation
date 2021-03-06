# species info 

# write a function which can decleare the variables needed for the analysis 
# tell it which species to get the info, it will make a bunch of objects 
# these will be available via TE.names

# need to add more species 

#cow 
#elephant
#horse
#pig
#rabbit




rep_info <- function(spec1, genome){
  
  # path to each species repeat masker files downloaded from UCSC and saved as "genome"_rmsk.txt
  if(!(spec1 == "Pig" | spec1 == "Rabbit" | spec1 == "Horse")){
    rep_name <- paste( "~/Desktop/retrotransposonAccumulationAnalysis/accesoryFiles/Data/rmskTable/", genome,"_rmsk.txt" , sep = "")
    
    
    rep <- read.table(file= rep_name,
                      colClasses= c("integer", "integer", "integer", "integer", "integer",  # some stats
                                    "character", "integer", "integer",  # repeat coordinates
                                    "integer",  # genome left
                                    "character",   # strand
                                    "factor", "factor", "factor",     # classifications
                                    "integer", "integer", "integer", "integer"),        # more coordinates
                      col.names= c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                                   "genoName", "genoStart", "genoEnd", "genoLeft",
                                   "strand", "repName", "repClass", "repFamily", 
                                   "repStart", "repEnd", "repLeft", "id"),
                      comment.char="#",
                      header = FALSE, fill = TRUE
    )
  }else{
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
    
    rep$repFamily <- unlist(lapply(lapply(strsplit(rep$repClass, "/"), rev), getElement, 1))
    rep$repClass <- unlist(lapply(strsplit(rep$repClass, "/"),getElement, 1))
    rep$bin <- 0
    rep[,c("milliDiv", "milliIns", "milliDel")] <- rep[,c("milliDiv", "milliIns", "milliDel")] * 10
    rep <- rep[,c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                  "genoName", "genoStart", "genoEnd", "genoLeft",
                  "strand", "repName", "repClass", "repFamily", 
                  "repStart", "repEnd", "repLeft", "id")]
  }
  
  MIR <- rep[grep("MIR", rep$repFamily),]
  L2 <- rep[grep("L2", rep$repFamily),]
  L1 <- rep[grep("L1", rep$repFamily),]
  L1ME <- L1[grep("ME", L1$repName),]
  L1MD <- L1[grep("MD", L1$repName),]
  L1MC <- L1[grep("MC", L1$repName),]
  L1MB <- L1[grep("MB", L1$repName),]
  L1MA <- L1[grep("MA", L1$repName),]
  
  if(spec1 == "Human"){
    AluJ <- rep[grep("AluJ", rep$repName),]
    AluS <- rep[grep("AluS", rep$repName),]
    AluY <- rep[grep("AluY", rep$repName),]
    L1PB <- L1[grep("PB", L1$repName),]
    L1PA <- L1[grep("PA", L1$repName),]
    L1HS <- L1[grep("HS", L1$repName),] 
  }
  
  if(spec1 == "Mouse"){
    ALU <- rep[grep("Alu", rep$repFamily),]
    PB <- ALU[grep("PB", ALU$repName),]
    B1 <- rbind(ALU[grep("B1_", ALU$repName),], ALU[grep("B1F", ALU$repName),])
    B2all <- rep[grep(pattern="B2", rep$repFamily),]
    B2 <- B2all[grep("B2", B2all$repName),]
    B3 <- B2all[grep("B3", B2all$repName),]
    B4all <- rep[grep(pattern="B4", rep$repFamily),]
    B4 <- B4all[-(grep(pattern="RSINE", B4all$repName)),]
    L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),]
    Lx <- L1new[grep("Lx", L1new$repName),]
    L1Md <- L1new[grep("L1Md", L1new$repName),]
    L1_Mus <- L1new[grep("L1_Mus", L1new$repName),]
    L1_Mur <- L1new[grep("L1_Mur", L1new$repName),]
    L1_Mm <- L1new[grep("L1_Mm", L1new$repName),]
  }
  
  if(spec1 == "Chimp"){
    Alu <- rep[grep("Alu", rep$repFamily),]
    AluJ <- Alu[grep("AluJ", Alu$repName),]
    AluS <- Alu[grep("AluS", Alu$repName),]
    AluY <- Alu[grep("AluY", Alu$repName),]
    L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),] 
    L1PB <- L1new[grep("PB", L1new$repName),]
    L1PA <- L1new[grep("PA", L1new$repName),]
    L1Pt <- L1new[grep("Pt", L1new$repName),] 
  }
  
  if(spec1 == "Rhesus"){
    Alu <- rep[grep("Alu", rep$repFamily),]
    AluJ <- Alu[grep("AluJ", Alu$repName),]
    AluS <- Alu[grep("AluS", Alu$repName),]
    AluY <- Alu[grep("AluY", Alu$repName),]
    L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),] 
    L1PB <- L1new[grep("PB", L1new$repName),]
    L1PA <- L1new[grep("PA", L1new$repName),]
    L1RS <- rbind(L1new[grep("_RS", L1new$repName),], L1new[grep("PREC", L1new$repName),])
  }
  
  if(spec1 == "Dog"){
    Lys <- rep[grep("tRNA-Lys", rep$repFamily),]
    SINEC_c <- Lys[grep("SINEC_c", Lys$repName),]
    SINEC_b <- Lys[grep("SINEC_b", Lys$repName),]
    SINEC_a <- Lys[grep("SINEC_a", Lys$repName),]
    SINEC_old <- Lys[grep("SINEC_old", Lys$repName),]
    SINEC_Cf <- Lys[grep("SINEC_Cf", Lys$repName),]
    L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),]
    L1_Carn <- L1new[grep("L1_Carn", L1new$repName), ]
    L1_Canid <- L1new[grep("L1_Canid", L1new$repName), ]
    L1_Canis <- L1new[grep("L1_Canis", L1new$repName), ]
    L1_Cf <- L1new[grep("L1_Cf", L1new$repName), ]
  }
  
  if(spec1 == "Cow"){
    BovB = rep[rep$repFamily == "RTE-BovB",]
    BovA = rep[rep$repFamily == "BovA",]
    tRNA_Glu = rep[rep$repFamily == "tRNA-Glu",]
    L1_BT <- rbind(L1[grep("BT", L1$repName),], L1[grep("Art", L1$repName),])
  }
  
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
  
  if(spec1 == "Rabbit"){
    
    CSINE1 = rep[rep$repName == "CSINE1",]
    CSINE2 = rep[rep$repName == "CSINE2",]
    CSINE2B = rep[rep$repName == "CSINE2B",] 
    CSINE3A = rep[rep$repName == "CSINE3A",] 
    
    L1A_OC = rep[grep("L1A_O", rep$repName),] 
    L1A2_OC = rep[grep("L1A2_O", rep$repName),] 
    L1C_OC = rep[grep("L1C_O", rep$repName),] 

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
  
  
  # just get some more TE clases over the genome so we can get a more refined look at whats happening
  if(spec1 == "Human"){
    TE.names <- c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1HS", "L2", "MIR")
  }
  if(spec1 == "Mouse"){
    TE.names <- c("B1", "B2", "B3", "B4", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "Lx", "L1_Mur", "L1_Mus","L1Md" ,"L1_Mm", "L2", "MIR")
  }
  if(spec1 == "Chimp"){
    TE.names <- c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1Pt", "L2", "MIR")
  }
  if(spec1 == "Rhesus"){
    TE.names <- c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1RS", "L2", "MIR")
  }
  if(spec1 == "Dog"){
    TE.names <- c("SINEC_old", "SINEC_c", "SINEC_b", "SINEC_a", "SINEC_Cf", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_Carn", "L1_Canid", "L1_Canis", "L1_Cf", "L2", "MIR")
  }
  if(spec1 == "Cow"){
    TE.names <- c("tRNA_Glu", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_BT", "L2", "MIR", "BovA", "BovB")
  }
  if(spec1 == "Pig"){
    TE.names <- c("Pre0_SS", "PRE1", "SINE1_SS", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_SS", "L1_SS_1", "L1_SS_2", "L1_SS_3", "L2", "MIR")
  }
  if(spec1 == "Rabbit"){
    TE.names <- c("CSINE1", "CSINE2", "CSINE2B", "CSINE3A", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1A_OC", "L1A2_OC", "L1C_OC", "L2", "MIR")
  }
  if(spec1 == "Horse"){
    TE.names <- c("ERE1", "ERE2", "ERE3", "ERE4", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_EC_1", "L1_EC_2", "L1_EC_3", "L1_EC_4", "L1_EC_5", "L2", "MIR")
  }
  
  sorted.TE <- NULL
  for(i in TE.names){
    sorted.TE = c(sorted.TE, list(get(i)))
  }
  names(sorted.TE) <- TE.names
  
  return(sorted.TE)
  
}

#### a function where we can get tables that describe the family structure for each TE

repFamStruct <- function(spec1){
  if(spec1 == "Human"){
    repFamStruct <- data.frame(TEname = c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1HS", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Mouse"){
    repFamStruct <- data.frame(TEname = c("B1", "B2", "B3", "B4", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "Lx", "L1_Mur", "L1_Mus","L1Md" ,"L1_Mm", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Chimp"){
    repFamStruct <- data.frame(TEname = c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1Pt", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Rhesus"){
    repFamStruct <- data.frame(TEname = c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1RS", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")    
                               )
  }
  if(spec1 == "Dog"){
    repFamStruct <- data.frame(TEname = c("SINEC_old", "SINEC_c", "SINEC_b", "SINEC_a", "SINEC_Cf", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_Carn", "L1_Canid", "L1_Canis", "L1_Cf", "L2", "MIR"),
                              repType = c("new_SINE", "new_SINE", "new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
                              )
  }
  if(spec1 == "Cow"){
    repFamStruct <- data.frame(TEname = c("tRNA_Glu", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_BT", "L2", "MIR", "BovA", "BovB"),
                               repType = c("new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE","new_LINE" ,"new_LINE", "ancient", "ancient", "new_SINE_RTE", "new_LINE_RTE")
    )
  }
  if(spec1 == "Pig"){
    repFamStruct <- data.frame(TEname = c("Pre0_SS", "PRE1", "SINE1_SS", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_SS", "L1_SS_1", "L1_SS_2", "L1_SS_3", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE","new_LINE" ,"new_LINE" ,"new_LINE" ,"new_LINE" ,"new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Rabbit"){
    repFamStruct <- data.frame(TEname =  c("CSINE1", "CSINE2", "CSINE2B", "CSINE3A", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1A_OC", "L1A2_OC", "L1C_OC", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE","new_LINE" ,"new_LINE" ,"new_LINE" ,"new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Horse"){
    repFamStruct <- data.frame(TEname =  c("ERE1", "ERE2", "ERE3", "ERE4", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_EC_1", "L1_EC_2", "L1_EC_3", "L1_EC_4", "L1_EC_5", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE","new_LINE" ,"new_LINE" ,"new_LINE" ,"new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
    )
  }
  
  return(repFamStruct)
}




