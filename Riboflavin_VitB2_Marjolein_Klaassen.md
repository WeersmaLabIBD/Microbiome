**Effects of vitamin B2 (Riboflavin) supplementation on the gut microbiome of CD Patients**

**Authors: Marjolein Klaassen**

**Date: ** 

After having found no significant changes in the relative abundances of gut bacterial taxonomies and pathways after three weeks of riboflavin supplementation (vitamin B2) in CD patients, we decided to test whether the ratio between the relative abundance of F. prausnitzii/E. coli changes after three weeks. Below, one can appreciate the R-code for this. 
Bio-informatical tools MetaPhlAn and HUMAnN2 were used to determine the presence/absence and relative abundances of microbial taxonomies and pathways, from the metagenomic data. Since these tools produce relative abundances (and inferring absolute abundances from the number of reads seems inaccurate), I have calculated the ratios based on relative abundances. However, since we are calculating ratios, I believe this is no problem.

 
# Setting the working directory

setwd("~/Documents/IBD Weersma/Vitamin B2/Working directory")

# Importing microbiome taxonomy data 
db_VitB2 = read.csv("./metaphlanmerged.txt", header = T, sep = "\t", stringsAsFactors = F)
db_VitB2 = as.data.frame(db_VitB2)

db_Clin_Jul = read.csv("Rise-UpDB.csv", header = T, sep = ",", stringsAsFactors = F)
db_Clin_Jul = db_Clin_Jul[,c("Pnumber", "Age", "Sex", "Lab1Calprotectin", "Lab2Calprotectin", "MontrealL", "HBI1", "HBI2", "Colectomy", "PPI")]


# Making the first column (i.e. microbiome features) the rownames. 
rownames(db_VitB2) = db_VitB2$ID
# Removing the first column. 
db_VitB2$ID = NULL


# Only selecting the taxonomical level we are interested in (all). 
# The 'filterMetaGenomeDF' function is written by Ranko Gacesa in his R-scripts for microbiome data. Firstly, this has
# to run. 

#_____________________________ R-script part of Ranko Gacesa _____________________________
filterMetaGenomeDF <- function(inDF,presPerc = 0.1,minMRelAb = 0.5,minMedRelAb=0.0, rescaleTaxa=F,verbose=T,
                               keepDomains=c('Bacteria'),
                               keepLevels=c('T','S','G','F','O','C','P')) {
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  tColsNMG = grep('k__',colnames(inDF),invert = T) # colums with microbiome
  # replaces NAs in microbiome with 0s
  for (c in tCols) {
    inDF[,c][is.na(inDF[,c])] <- 0.0
  }
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    nrnZ = as.numeric(sum(inDF[,c]!=0.0))
    if ( (nrnZ/as.numeric(nrow(inDF))) < presPerc) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; nr non Zero:',nrnZ,'=',nrnZ/as.numeric(nrow(inDF)),'>> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = mean(inDF[,c])
    if ( mn < minMRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; mean rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }
  
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = median(inDF[,c])
    if ( mn < minMedRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; median rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }  
  
  # keep domains
  # -----------------------------
  if (keepDomains!='All' & keepDomains !="") {
    toKeep = grep('k__',colnames(inDF),invert = T)
    for (d in keepDomains) {
      #d print(d)
      toKeep = c(toKeep,grep(paste('k__',d,sep=''),colnames(inDF)))
    }
    inDF <- inDF[,toKeep]
  }
  
  #TODO
  
  
  # remove taxonomic levels
  # -----------------------------
  inDFnonTaxa <- as.data.frame(inDF[,grep('k__',colnames(inDF),invert=T)])
  colnames(inDFnonTaxa) <- colnames(inDF)[grep('k__',colnames(inDF),invert=T)]
  inDF2 <- inDF[,grep('k__',colnames(inDF),invert=F)]
  # pick strains (T)
  taxaTCols <- grep('t__',colnames(inDF2))
  taxaT <- inDF2[,taxaTCols]
  if (length(taxaT) > 0) {inDF2 <- inDF2[,-taxaTCols]}
  # pick species (S)
  taxaSCols <- grep('s__',colnames(inDF2))
  taxaS <- inDF2[,taxaSCols]
  if (length(taxaS) > 0) {inDF2 <- inDF2[,-taxaSCols]}
  # pick genera (G)
  taxaGCols <- grep('g__',colnames(inDF2))
  taxaG <- inDF2[,taxaGCols]
  if (length(taxaG) > 0) {inDF2 <- inDF2[,-taxaGCols]}
  # pick families (F)
  taxaFCols <- grep('f__',colnames(inDF2))
  taxaF <- inDF2[,taxaFCols]
  if (length(taxaF) > 0) {inDF2 <- inDF2[,-taxaFCols]}
  # pick orders (O)
  taxaOCols <- grep('o__',colnames(inDF2))
  taxaO <- inDF2[,taxaOCols]
  if (length(taxaO) > 0) {inDF2 <- inDF2[,-taxaOCols]}
  # pick classes (C)
  taxaCCols <- grep('c__',colnames(inDF2))
  taxaC <- inDF2[,taxaCCols]
  if (length(taxaC) > 0) {inDF2 <- inDF2[,-taxaCCols]}
  # pick phyla (P)
  taxaPCols <- grep('p__',colnames(inDF2))
  taxaPColsKeep <- grep('p__',colnames(inDF2),invert = T)
  taxaPColsKeepN <- colnames(inDF2)[grep('p__',colnames(inDF2),invert = T)]
  taxaP <- inDF2[,taxaPCols]
  if (length(taxaP) > 0) {inDF2 <- as.data.frame(inDF2[,taxaPColsKeep])}
  colnames(inDF2) <- taxaPColsKeepN
  taxaK <- inDF2
  # pick 
  oDF <- inDFnonTaxa
  if (verbose) {print ('Keeping following taxonomic levels:'); print(keepLevels)}
  if ('K' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaK),'Kingdoms'))}
    oDF <- cbind(oDF,taxaK)}
  if ('P' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaP),'Phyla'))} 
    oDF <- cbind(oDF,taxaP)}
  if ('C' %in% keepLevels) {
    if (verbose) {print(paste(' -> kept',ncol(taxaC),'Classes'))} 
    oDF <- cbind(oDF,taxaC)}
  if ('O' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaO),'Orders'))}
    oDF <- cbind(oDF,taxaO)}
  if ('F' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaF),'Families'))}
    oDF <- cbind(oDF,taxaF)}
  if ('G' %in% keepLevels) {if (verbose){ print(paste(' -> kept',ncol(taxaG),'Genera'))}
    oDF <- cbind(oDF,taxaG)}
  if ('S' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaS),'Species'))}
    oDF <- cbind(oDF,taxaS)}
  if ('T' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaT),'Strains'))}
    oDF <- cbind(oDF,taxaT)}
  if (verbose) {print ('data processing done, returning Dataframe')}
  oDF
}


#_____________________________ end of R-script of Ranko Gacesa _____________________________

tTrans = t.data.frame(db_VitB2)
db_VitB2_King = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = 0.0000001,minMedRelAb = -1,keepLevels = c("K"))
db_VitB2_King = as.data.frame(db_VitB2_King)
db_VitB2_Phylum = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = 0.0000001,minMedRelAb = -1,keepLevels = c("P"))
db_VitB2_Phylum = as.data.frame(db_VitB2_Phylum)
db_VitB2_Class = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = 0.0000001,minMedRelAb = -1,keepLevels = c("C"))
db_VitB2_Class = as.data.frame(db_VitB2_Class)
db_VitB2_Order = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = 0.0000001,minMedRelAb = -1,keepLevels = c("O"))
db_VitB2_Order = as.data.frame(db_VitB2_Order)
db_VitB2_Family = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = 0.0000001,minMedRelAb = -1,keepLevels = c("F"))
db_VitB2_Family = as.data.frame(db_VitB2_Family)
db_VitB2_Genus = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = 0.0000001,minMedRelAb = -1,keepLevels = c("G"))
db_VitB2_Genus = as.data.frame(db_VitB2_Genus)
db_VitB2_Species = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = 0.0000001,minMedRelAb = -1,keepLevels = c("S"))
db_VitB2_Species = as.data.frame(db_VitB2_Species)



for (c in c(1:nrow(db_VitB2_Species))) {
  db_VitB2_Species[c,] <- asin(sqrt(db_VitB2_Species[c,]/100.0))
}



###############################################################################
################# Sharp disease activity cut-off ##############################
###############################################################################

# Making t=1 and t=4 time groups 
db_VitB2_Species = as.data.frame(t(db_VitB2_Species))
t1 = db_VitB2_Species[,c(grep("M1", colnames(db_VitB2_Species)))]
rownames(t1) = rownames(db_VitB2_Species)
t4 = db_VitB2_Species[,c(grep("M4", colnames(db_VitB2_Species)))]
rownames(t4) = rownames(db_VitB2_Species)
t4 = as.data.frame(t4)

# Removing samples that have not both T1 and T4 (we only want to include samples that have both time points)
# These might change, for we have 3 T4's that have been send to the Broad, but have not came back from sequencing.  
colnames(t1) = gsub("_M1_metaphlan", "", colnames(t1))
colnames(t4) = gsub("_M4_metaphlan", "", colnames(t4))

InBothSamples = intersect(colnames(t1), colnames(t4))
t1 = t1[,InBothSamples]
t4 = t4[,InBothSamples]

tt1 = as.data.frame(t(t1))
tt1["Pnumber"] = rownames(tt1)
tt1 = tt1[,c(ncol(tt1), 1:ncol(tt1)-1)]

t0_Clin = merge(db_Clin_Jul, tt1, by="Pnumber", all = FALSE)
t0_Clin = t0_Clin[,c(1:4, 6,7, 9:ncol(t0_Clin))]
t0_Clin = t0_Clin[t0_Clin$Colectomy == "No",]

t0_Clin["DiseaseActivityAtBaseLine"] = NA
t0_Clin = t0_Clin[,c(ncol(t0_Clin), 1:ncol(t0_Clin)-1)]

#for (i in 1:nrow(t0_Clin)){
#  if(t0_Clin$Lab1Calprotectin[i] >= 200){
#    t0_Clin$DiseaseActivityAtBaseLine[i] = "Active_baseline"
#  } else if (t0_Clin$Lab1Calprotectin[i] <= 60 ){
#    t0_Clin$DiseaseActivityAtBaseLine[i] = "Remission_baseline"
#  } else if (t0_Clin$Lab1Calprotectin[i] < 200 | t0_Clin$Lab1Calprotectin[i] > 60){
#    t0_Clin$DiseaseActivityAtBaseLine[i] = "Exclude"
#  } else
#    t0_Clin$Pnumber[i] = t0_Clin$Pnumber[i]
#}

for (i in 1:nrow(t0_Clin)){
  if(t0_Clin$Lab1Calprotectin[i] > 200){ #we intentionally refer to calpotection from T1
    t0_Clin$DiseaseActivityAtBaseLine[i] = "Active_baseline"
  } else 
    t0_Clin$DiseaseActivityAtBaseLine[i] = "Remission_baseline"
}

t0_Clin["Pnumber_tp"] = t0_Clin$Pnumber
t0_Clin = t0_Clin[,c(ncol(t0_Clin), 1:ncol(t0_Clin)-1)]
t0_Clin$Pnumber_tp = paste0("T0_", t0_Clin$Pnumber_tp)
t0_Clin["MergeQC"] = t0_Clin$Pnumber
t0_Clin$MergeQC = paste0(t0_Clin$Pnumber, "_M1")
t0_Clin = t0_Clin[,c(ncol(t0_Clin), 1:ncol(t0_Clin)-1)]
write.table(t0_Clin, "PatientsMetagenomics_Riboflavin.csv", sep = "\t", quote = F, row.names = F)


tt4 = as.data.frame(t(t4))
tt4["Pnumber"] = rownames(tt4)
tt4 = tt4[,c(ncol(tt4), 1:ncol(tt4)-1)]

t3_Clin = merge (db_Clin_Jul, tt4, by="Pnumber", all= FALSE)
t3_Clin = t3_Clin[,c(1:3, 5,6, 8:ncol(t3_Clin))]
t3_Clin = t3_Clin[t3_Clin$Colectomy == "No",]

t3_Clin["Pnumber_tp"] = t3_Clin$Pnumber
t3_Clin = t3_Clin[,c(ncol(t3_Clin), 1:ncol(t3_Clin)-1)]
t3_Clin$Pnumber_tp = paste0("T3_", t3_Clin$Pnumber_tp)

t3_Clin["DiseaseActivityAtBaseLine"] = NA
t3_Clin = t3_Clin[,c(ncol(t3_Clin), 1:ncol(t3_Clin)-1)]

#for (i in 1:nrow(t3_Clin)){
#  if(t0_Clin$Lab1Calprotectin[i] >= 200){
#    t3_Clin$DiseaseActivityAtBaseLine[i] = "Active_baseline"
#  } else if (t0_Clin$Lab1Calprotectin[i] <= 60){
#    t3_Clin$DiseaseActivityAtBaseLine[i] = "Remission_baseline"
#  } else if (t0_Clin$Lab1Calprotectin[i] < 200 | t0_Clin$Lab1Calprotectin[i] > 60){
#    t3_Clin$DiseaseActivityAtBaseLine[i] = "Exclude"
#  } else
#    t3_Clin$Pnumber[i] = t3_Clin$Pnumber[i]
#}


for (i in 1:nrow(t3_Clin)){
  if(t0_Clin$Lab1Calprotectin[i] > 200){ #we intentionally refer to calpotection from T1
    t3_Clin$DiseaseActivityAtBaseLine[i] = "Active_baseline"
  } else 
    t3_Clin$DiseaseActivityAtBaseLine[i] = "Remission_baseline"
}

t3_Clin["MergeQC"] = t3_Clin$Pnumber
t3_Clin$MergeQC = paste0(t3_Clin$Pnumber, "_M4")
t3_Clin = t3_Clin[,c(ncol(t3_Clin), 1:ncol(t3_Clin)-1)]

#table(t0_Clin$DiseaseActivityAtBaseLine)
#Active at baseline n=29, remission at baseline n=38. 


t0_Clin["time_point"] = "M1"
t0_Clin = t0_Clin[,c(ncol(t0_Clin), 1:ncol(t0_Clin)-1)]

t3_Clin["time_point"] = "M3"
t3_Clin = t3_Clin[,c(ncol(t3_Clin), 1:ncol(t3_Clin)-1)]

colnames(t0_Clin)[8] = "LabCalprotectin"
colnames(t3_Clin)[8] = "LabCalprotectin"

colnames(t0_Clin)[10] = "HBI"
colnames(t3_Clin)[10] = "HBI"

# --> Creating dataframes of patients who have active disease at baseline, and a dataframe for patients 
# who have inactive disease at baseline. 
Active_Rib_M1 = t0_Clin[t0_Clin$DiseaseActivityAtBaseLine == "Active_baseline",]
ActiveRib_M3_wantM1 = t3_Clin[t3_Clin$DiseaseActivityAtBaseLine == "Active_baseline",]

Rem_Rib_M1 = t0_Clin[t0_Clin$DiseaseActivityAtBaseLine == "Remission_baseline",]
Rem_Rib_M3_wantM1 = t3_Clin[t3_Clin$DiseaseActivityAtBaseLine == "Remission_baseline",]

Act_basegroup_Species = rbind(Active_Rib_M1, ActiveRib_M3_wantM1)
Rem_basegroup_Species=  rbind(Rem_Rib_M1, Rem_Rib_M3_wantM1)

# Importing Readdepth (quality report)
# Because we would like to correct for read depth. 
QC = read.csv("QualityReportRiboflavin.csv", header = T, sep = ",", stringsAsFactors = F)
QC = QC[,c(2,8)]
QC["MergeQC"] = QC$Sample


Rem_Species_fin = merge(QC, Rem_basegroup_Species, by ="MergeQC", all = F)
Rem_Species_fin = Rem_Species_fin[,c(5, 6, 4, 3, 7:ncol(Rem_Species_fin))]
Rem_Species_fin$Total.Reads = gsub(",", "", Rem_Species_fin$Total.Reads)
Rem_Species_fin$Total.Reads = as.numeric(Rem_Species_fin$Total.Reads)

# First, tidy the column names 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = ":", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = " ", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = "-", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = ")", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = "/", replacement = "_")
MetaCycVTTidy = make.names(colnames(Rem_Species_fin), unique = TRUE)
colnames(Rem_Species_fin) = MetaCycVTTidy 

Rem_Species_fin["AbundancesFprau"] = NA
Rem_Species_fin["AbundancesEcoli"] = NA
Rem_Species_fin = Rem_Species_fin[,c(177,176,1:175)]


for (i in 1:nrow(Rem_Species_fin)){
  if (Rem_Species_fin$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii[i]==0){
    Rem_Species_fin$AbundancesFprau[i] = "no"
  } else 
    Rem_Species_fin$AbundancesFprau[i] = "yes"
}

Rem_Species_fin = Rem_Species_fin[,c(144, 1:143, 145:177)]

# Getting column name number based on column name
#which(colnames(Rem_Species_fin)=="k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii" )
#144
#which(colnames(Rem_Species_fin)=="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli" )
#174


Rem_Species_fin = Rem_Species_fin[,c(174, 2, 1, 3:173, 175:177)]

for (i in 1:nrow(Rem_Species_fin)){
  if (Rem_Species_fin$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli[i]==0){
    Rem_Species_fin$AbundancesEcoli[i] = "no"
  } else 
    Rem_Species_fin$AbundancesEcoli[i] = "yes"
}

Rem_Species_T0Prau = Rem_Species_fin[Rem_Species_fin$time_point=="M1",]
table(Rem_Species_T0Prau$AbundancesFprau)

Rem_Species_T3Prau = Rem_Species_fin[Rem_Species_fin$time_point=="M3",]
table(Rem_Species_T3Prau$AbundancesFprau)


Rem_Species_T0Coli = Rem_Species_fin[Rem_Species_fin$time_point=="M1",]
Rem_Species_T0Coli = Rem_Species_T0Coli[,c(2,5,1:3, 6:177)]
table(Rem_Species_T0Prau$AbundancesEcoli)

Rem_Species_T3Coli = Rem_Species_fin[Rem_Species_fin$time_point=="M3",]
Rem_Species_T3Coli = Rem_Species_T3Coli[,c(2,5,1:3, 6:177)]
table(Rem_Species_T3Prau$AbundancesEcoli)


Rem_KNHH = Rem_Species_fin[,c(5, 7, 4, 2, 1, 3)]
write.table(Rem_KNHH, "RemmissionRiboflavinAbundancesPrauColi.tsv", sep = "\t", quote = F, row.names = F)


# Dus, nu hebben we van elk sample (per tijdspunt), of het wel of niet E.coli en F.prau heeft. Nu gaan we alleen de samples
# houden, die allebei hebben (F.prau en E.coli). Daarna pas alleen de samples die allebei de tijdspunten hebben. 

Rem_Species_both = Rem_Species_fin[which(Rem_Species_fin$AbundancesEcoli=="yes" & Rem_Species_fin$AbundancesFprau=="yes"),]











###########################
######### ACT #############
###########################
Act_Species_fin = merge(QC, Act_basegroup_Species, by ="MergeQC", all = F)
Act_Species_fin = Act_Species_fin[,c(5, 6, 4, 3, 7:ncol(Act_Species_fin))]
Act_Species_fin$Total.Reads = gsub(",", "", Act_Species_fin$Total.Reads)
Act_Species_fin$Total.Reads = as.numeric(Act_Species_fin$Total.Reads)

Act_Species_fin["Ratio_Prau_Coli"] = NA
Act_Species_fin = Act_Species_fin[,c(176, 1:175)]

# First, tidy the column names 
names(Act_Species_fin) = gsub(x = names(Act_Species_fin), pattern = ":", replacement = "_") 
names(Act_Species_fin) = gsub(x = names(Act_Species_fin), pattern = " ", replacement = "_") 
names(Act_Species_fin) = gsub(x = names(Act_Species_fin), pattern = "-", replacement = "_") 
names(Act_Species_fin) = gsub(x = names(Act_Species_fin), pattern = ")", replacement = "_") 
names(Act_Species_fin) = gsub(x = names(Act_Species_fin), pattern = "/", replacement = "_")
MetaCycVTTidy = make.names(colnames(Act_Species_fin), unique = TRUE)
colnames(Act_Species_fin) = MetaCycVTTidy 


Act_Species_fin["AbundancesFprau"] = NA
Act_Species_fin["AbundancesEcoli"] = NA
Act_Species_fin = Act_Species_fin[,c(178,177,1:176)]


for (i in 1:nrow(Act_Species_fin)){
  if (Act_Species_fin$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii[i]==0){
    Act_Species_fin$AbundancesFprau[i] = "no"
  } else 
    Act_Species_fin$AbundancesFprau[i] = "yes"
}

Act_Species_fin = Act_Species_fin[,c(145, 1:144, 146:178)]

# Getting column name number based on column name
#which(colnames(Act_Species_fin)=="k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii" )
#144
#which(colnames(Act_Species_fin)=="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli" )
#175


Act_Species_fin = Act_Species_fin[,c(175, 2, 1, 3:174, 176:178)]

for (i in 1:nrow(Act_Species_fin)){
  if (Act_Species_fin$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli[i]==0){
    Act_Species_fin$AbundancesEcoli[i] = "no"
  } else 
    Act_Species_fin$AbundancesEcoli[i] = "yes"
}

Act_Species_T0Prau = Act_Species_fin[Act_Species_fin$time_point=="M1",]
table(Act_Species_T0Prau$AbundancesFprau)

Act_Species_T3Prau = Act_Species_fin[Act_Species_fin$time_point=="M3",]
table(Act_Species_T3Prau$AbundancesFprau)


Act_Species_T0Coli = Act_Species_fin[Act_Species_fin$time_point=="M1",]
Act_Species_T0Coli = Act_Species_T0Coli[,c(2,6,1,3:5, 7:178)]
table(Act_Species_T0Prau$AbundancesEcoli)

Act_Species_T3Coli = Act_Species_fin[Act_Species_fin$time_point=="M3",]
Act_Species_T3Coli = Act_Species_T3Coli[,c(2,6,1,3:5, 7:178)]
table(Act_Species_T3Prau$AbundancesEcoli)



Act_KNHH = Act_Species_fin[,c(6, 8, 4, 2, 1, 3)]
write.table(Act_KNHH, "ActiveRiboflavinAbundancesPrauColi.tsv", sep = "\t", quote = F, row.names = F)


# Dus, nu hebben we van elk sample (per tijdspunt), of het wel of niet E.coli en F.prau heeft. Nu gaan we alleen de samples
# houden, die allebei hebben (F.prau en E.coli). Daarna pas alleen de samples die allebei de tijdspunten hebben. 

Act_Species_both = Act_Species_fin[which(Act_Species_fin$AbundancesEcoli=="yes" & Act_Species_fin$AbundancesFprau=="yes"),]



Ratio_PrauColi_rem = Rem_Species_both[Rem_Species_both$Pnumber_tp %in% c("T0_p003", "T3_p003", "T0_p004", "T3_p004", "T0_p010", "T3_p010", "T0_p018", "T3_p018", "T0_p023", "T3_p023", "T0_p037", "T3_p037", "T0_p041", "T3_p041", "T0_p058", "T3_p058", "T0_p058", "T3_p058", "T0_p063", "T3_p063"),]

Ratio_PrauColi_act = Act_Species_both[Act_Species_both$Pnumber_tp %in% c("T0_p006", "T3_p006", "T0_p024", "T3_p024", "T0_p028", "T3_p028", "T0_p030", "T3_p030", "T0_p044", "T3_p044", "T0_p047", "T3_p047", "T0_p061", "T3_p061", "T0_p065", "T3_p065", "T0_p066", "T3_p066", "T0_p070", "T3_p070", "T0_p074", "T3_p074", "T0_p077", "T3_p077", "T0_p078", "T3_p078", "T0_p079", "T3_p079"),]



# Ratio F.prau/E.coli berekenen. Ik deel f.prau door e.coli, dus als er een groot getal uitkomt, dan is er relatief veel 
# f.prau in vergelijking tot e.coli. Komt er een laag getal uit, dan is er relatief veel e.coli itt tot F.prau. 
for (i in 1:nrow(Ratio_PrauColi_act)){
  if (Ratio_PrauColi_act$DiseaseActivityAtBaseLine == "Active_baseline"){
    Ratio_PrauColi_act$Ratio_Prau_Coli[i] = (Ratio_PrauColi_act$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii[i])/(Ratio_PrauColi_act$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli[i])
  } else 
    Ratio_PrauColi_act$DiseaseActivityAtBaseLine[i] = Ratio_PrauColi_act$DiseaseActivityAtBaseLine[i]
}


## BMI
Act_M1_Ratio = Ratio_PrauColi_act[Ratio_PrauColi_act$time_point=="M1",]
Act_M3_Ratio = Ratio_PrauColi_act[Ratio_PrauColi_act$time_point=="M3",]

Act_m1 = Act_M1_Ratio$Ratio_Prau_Coli
Act_m3 = Act_M3_Ratio$Ratio_Prau_Coli

Act_Rat = data.frame(Y=c(Act_m1, Act_m3), Site=factor(rep(c("Act_m1", "Act_m3"), times=c(length(Act_m1), length(Act_m3)))))
y = Act_Rat$Y
x = Act_Rat$Site
wilcox.test(y ~ x, data=Act_Rat) 

#p=0.6347
mean(Act_M1_Ratio$Ratio_Prau_Coli)
mean(Act_M3_Ratio$Ratio_Prau_Coli)












###### Remission ########
# Ratio F.prau/E.coli berekenen. Ik deel f.prau door e.coli, dus als er een groot getal uitkomt, dan is er relatief veel 
# f.prau in vergelijking tot e.coli. Komt er een laag getal uit, dan is er relatief veel e.coli itt tot F.prau. 
Ratio_PrauColi_rem["Ratio_Prau_Coli"] = NA
Ratio_PrauColi_rem = Ratio_PrauColi_rem[,c(178, 1:177)]

for (i in 1:nrow(Ratio_PrauColi_rem)){
  if (Ratio_PrauColi_rem$DiseaseActivityAtBaseLine == "Remission_baseline"){
    Ratio_PrauColi_rem$Ratio_Prau_Coli[i] = (Ratio_PrauColi_rem$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii[i])/(Ratio_PrauColi_rem$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli[i])
  } else 
    Ratio_PrauColi_rem$DiseaseActivityAtBaseLine[i] = Ratio_PrauColi_rem$DiseaseActivityAtBaseLine[i]
}


## BMI
Rem_M1_Ratio = Ratio_PrauColi_rem[Ratio_PrauColi_rem$time_point=="M1",]
Rem_M3_Ratio = Ratio_PrauColi_rem[Ratio_PrauColi_rem$time_point=="M3",]

Rem_m1 = Rem_M1_Ratio$Ratio_Prau_Coli
Rem_m3 = Rem_M3_Ratio$Ratio_Prau_Coli

Rem_Rat = data.frame(Y=c(Rem_m1, Rem_m3), Site=factor(rep(c("Rem_m1", "Rem_m3"), times=c(length(Rem_m1), length(Rem_m3)))))
y = Rem_Rat$Y
x = Rem_Rat$Site
wilcox.test(y ~ x, data=Rem_Rat) 

#p=0.3865
mean(Rem_M1_Ratio$Ratio_Prau_Coli)
mean(Rem_M3_Ratio$Ratio_Prau_Coli)


