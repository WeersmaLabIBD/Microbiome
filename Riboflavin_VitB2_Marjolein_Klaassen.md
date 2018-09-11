**Effects of vitamin B2 (Riboflavin) supplementation on the gut microbiome of CD Patients**

**Authors: Marjolein Klaassen**

**Date: 11-09-2018** 

After having found no significant changes in the relative abundances of gut bacterial taxonomies and pathways after three weeks of riboflavin supplementation (vitamin B2) in CD patients, we decided to test whether the ratio between the relative abundance of F. prausnitzii/E. coli changes after three weeks. Below, one can appreciate the R-code for this. 
Bio-informatical tools MetaPhlAn and HUMAnN2 were used to determine the presence/absence and relative abundances of microbial taxonomies and pathways, from the metagenomic data. Since these tools produce relative abundances (and inferring absolute abundances from the number of reads seems inaccurate), I have calculated the ratios based on relative abundances. However, since we are calculating ratios, I believe this is no problem.

 
**Setting the working directory**
```
setwd("~/Documents/IBD Weersma/Vitamin B2/Working directory")
```

**Importing microbiome taxonomy data**
```
db_VitB2 = read.csv("./metaphlanmerged.txt", header = T, sep = "\t", stringsAsFactors = F)
db_VitB2 = as.data.frame(db_VitB2)
```

**Importing the database from Julius von Martels and Arno Bourgonje**
```
db_Clin_Jul = read.csv("Rise-UpDB.csv", header = T, sep = ",", stringsAsFactors = F)
db_Clin_Jul = db_Clin_Jul[,c("Pnumber", "Age", "Sex", "Lab1Calprotectin", "Lab2Calprotectin", "MontrealL", "HBI1", "HBI2", "Colectomy", "PPI")]
```

**Making the sampleID's the rownames for clearity.** 
```
rownames(db_VitB2) = db_VitB2$ID
```

**Removing the first column.** 
```
db_VitB2$ID = NULL
```

**Only selecting the taxonomical levels we are interested in.The 'filterMetaGenomeDF' function is written by Ranko Gacesa in his R-scripts for microbiome data.** 
```

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

```

# Creating databases with different levels of taxonomies. 
```
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
```

**Performing arc sine square root transformations on the relative abundances of the database containing species level**
```
for (c in c(1:nrow(db_VitB2_Species))) {
  db_VitB2_Species[c,] <- asin(sqrt(db_VitB2_Species[c,]/100.0))
}
```

**I only want to include patients in the analyses, when they have a microbiome sample taken both at T=0weeks and T=3 weeks. Therefore, I first want to create dataframes with all the samples at T=0 and one with all the samples at T=3, to test which sampleIDs are present in both databases.** 
```
db_VitB2_Species = as.data.frame(t(db_VitB2_Species))
t1 = db_VitB2_Species[,c(grep("M1", colnames(db_VitB2_Species)))]
rownames(t1) = rownames(db_VitB2_Species)
t4 = db_VitB2_Species[,c(grep("M4", colnames(db_VitB2_Species)))]
rownames(t4) = rownames(db_VitB2_Species)
t4 = as.data.frame(t4)
```

**Removing samples that have not both T1 and T4 (we only want to include samples that have both time points.**
```
colnames(t1) = gsub("_M1_metaphlan", "", colnames(t1))
colnames(t4) = gsub("_M4_metaphlan", "", colnames(t4))

InBothSamples = intersect(colnames(t1), colnames(t4))
t1 = t1[,InBothSamples]
t4 = t4[,InBothSamples]
```
**Now, I want to merge the clinical data (clinical database of Julius) with the taxonomy data (t1), at T=0 weeks.** 
```
tt1 = as.data.frame(t(t1))
tt1["Pnumber"] = rownames(tt1)
tt1 = tt1[,c(ncol(tt1), 1:ncol(tt1)-1)]

t0_Clin = merge(db_Clin_Jul, tt1, by="Pnumber", all = FALSE)
t0_Clin = t0_Clin[,c(1:4, 6,7, 9:ncol(t0_Clin))]
```
**Removing patients who have had a colectomy (i.e. stoma patients)**
```
t0_Clin = t0_Clin[t0_Clin$Colectomy == "No",]
```
**Creating a new column defining disease activity based on calprotectin level at baseline.** 
```
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
```
```
for (i in 1:nrow(t0_Clin)){
  if(t0_Clin$Lab1Calprotectin[i] > 200){ #we intentionally refer to calpotection from T1
    t0_Clin$DiseaseActivityAtBaseLine[i] = "Active_baseline"
  } else 
    t0_Clin$DiseaseActivityAtBaseLine[i] = "Remission_baseline"
}
```

**I'm creating an extra column with sampleID information, so that later analyses will be easier.**
```
t0_Clin["Pnumber_tp"] = t0_Clin$Pnumber
t0_Clin = t0_Clin[,c(ncol(t0_Clin), 1:ncol(t0_Clin)-1)]
t0_Clin$Pnumber_tp = paste0("T0_", t0_Clin$Pnumber_tp)
t0_Clin["MergeQC"] = t0_Clin$Pnumber
t0_Clin$MergeQC = paste0(t0_Clin$Pnumber, "_M1")
t0_Clin = t0_Clin[,c(ncol(t0_Clin), 1:ncol(t0_Clin)-1)]
```

**I want to do the exact same for the samples at T=3, i.e. to merge the clinical data (clinical database of Julius) with the taxonomy data (t1).**
```
tt4 = as.data.frame(t(t4))
tt4["Pnumber"] = rownames(tt4)
tt4 = tt4[,c(ncol(tt4), 1:ncol(tt4)-1)]

t3_Clin = merge (db_Clin_Jul, tt4, by="Pnumber", all= FALSE)
t3_Clin = t3_Clin[,c(1:3, 5,6, 8:ncol(t3_Clin))]
```
**Removing all patients who had a colectomy, i.e. stoma patients**
```
t3_Clin = t3_Clin[t3_Clin$Colectomy == "No",]
```

**Adding a column of SampleID information**
```
t3_Clin["Pnumber_tp"] = t3_Clin$Pnumber
t3_Clin = t3_Clin[,c(ncol(t3_Clin), 1:ncol(t3_Clin)-1)]
t3_Clin$Pnumber_tp = paste0("T3_", t3_Clin$Pnumber_tp)
```
**Creating a column defining disease activity at baseline, based at baseline calprotectin levels** 
```
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
```
**Adding extra column of SampleID information, for later merging with quality of reads**
```
t3_Clin["MergeQC"] = t3_Clin$Pnumber
t3_Clin$MergeQC = paste0(t3_Clin$Pnumber, "_M4")
t3_Clin = t3_Clin[,c(ncol(t3_Clin), 1:ncol(t3_Clin)-1)]
```
```
#table(t0_Clin$DiseaseActivityAtBaseLine)
#Active at baseline n=29, remission at baseline n=38. 
```

**To merge the T=0 and T=3 samples, all column names need to be matched exactly. This is what I am doing here.** 
```
t0_Clin["time_point"] = "M1"
t0_Clin = t0_Clin[,c(ncol(t0_Clin), 1:ncol(t0_Clin)-1)]

t3_Clin["time_point"] = "M3"
t3_Clin = t3_Clin[,c(ncol(t3_Clin), 1:ncol(t3_Clin)-1)]

colnames(t0_Clin)[8] = "LabCalprotectin"
colnames(t3_Clin)[8] = "LabCalprotectin"

colnames(t0_Clin)[10] = "HBI"
colnames(t3_Clin)[10] = "HBI"
```

**--> Creating dataframes of patients who have active disease at baseline, and a dataframe for patients who have inactive disease at baseline.** 
```
Active_Rib_M1 = t0_Clin[t0_Clin$DiseaseActivityAtBaseLine == "Active_baseline",]
ActiveRib_M3_wantM1 = t3_Clin[t3_Clin$DiseaseActivityAtBaseLine == "Active_baseline",]

Rem_Rib_M1 = t0_Clin[t0_Clin$DiseaseActivityAtBaseLine == "Remission_baseline",]
Rem_Rib_M3_wantM1 = t3_Clin[t3_Clin$DiseaseActivityAtBaseLine == "Remission_baseline",]

Act_basegroup_Species = rbind(Active_Rib_M1, ActiveRib_M3_wantM1)
Rem_basegroup_Species=  rbind(Rem_Rib_M1, Rem_Rib_M3_wantM1)
```

**Importing Readdepth (quality report). Because we would like to correct for read depth.** 
```
QC = read.csv("QualityReportRiboflavin.csv", header = T, sep = ",", stringsAsFactors = F)
QC = QC[,c(2,8)]
QC["MergeQC"] = QC$Sample
```
**Merging previous data column with quality of reads**
```
Rem_Species_fin = merge(QC, Rem_basegroup_Species, by ="MergeQC", all = F)
Rem_Species_fin = Rem_Species_fin[,c(5, 6, 4, 3, 7:ncol(Rem_Species_fin))]
Rem_Species_fin$Total.Reads = gsub(",", "", Rem_Species_fin$Total.Reads)
Rem_Species_fin$Total.Reads = as.numeric(Rem_Species_fin$Total.Reads)
```

**Removing all noise from column names**
```
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = ":", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = " ", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = "-", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = ")", replacement = "_") 
names(Rem_Species_fin) = gsub(x = names(Rem_Species_fin), pattern = "/", replacement = "_")
MetaCycVTTidy = make.names(colnames(Rem_Species_fin), unique = TRUE)
colnames(Rem_Species_fin) = MetaCycVTTidy 
```

**Creating new columns that will categorically describe abundance of F.prau and E.coli (yes/no)**
```
Rem_Species_fin["AbundancesFprau"] = NA
Rem_Species_fin["AbundancesEcoli"] = NA
Rem_Species_fin = Rem_Species_fin[,c(177,176,1:175)]
```
**Loop saying that if F.prau =0, then the categorical abundance column =0. Otherwise, yes**
```
for (i in 1:nrow(Rem_Species_fin)){
  if (Rem_Species_fin$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii[i]==0){
    Rem_Species_fin$AbundancesFprau[i] = "no"
  } else 
    Rem_Species_fin$AbundancesFprau[i] = "yes"
}
Rem_Species_fin = Rem_Species_fin[,c(144, 1:143, 145:177)]
```
```
# Getting column name number based on column name
#which(colnames(Rem_Species_fin)=="k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii" )
#144
#which(colnames(Rem_Species_fin)=="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli" )
#174


Rem_Species_fin = Rem_Species_fin[,c(174, 2, 1, 3:173, 175:177)]
```
**Loop saying that if E.coli =0, then the categorical abundance column =0. Otherwise, yes**
```
for (i in 1:nrow(Rem_Species_fin)){
  if (Rem_Species_fin$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli[i]==0){
    Rem_Species_fin$AbundancesEcoli[i] = "no"
  } else 
    Rem_Species_fin$AbundancesEcoli[i] = "yes"
}
```




**Ratio F.prau/E.coli berekenen. Ik deel f.prau door e.coli, dus als er een groot getal uitkomt, dan is er relatief veel f.prau in vergelijking tot e.coli. Komt er een laag getal uit, dan is er relatief veel e.coli itt tot F.prau.** 
```
Rem_Species_fin["Ratio_Prau_Coli"] = NA
Rem_Species_fin = Rem_Species_fin[,c(178, 1:177)]

for (i in 1:nrow(Rem_Species_fin)){
  if (Rem_Species_fin$DiseaseActivityAtBaseLine == "Remission_baseline"){
    Rem_Species_fin$Ratio_Prau_Coli[i] = (Rem_Species_fin$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii[i])/(Rem_Species_fin$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli[i])
  } else 
    Rem_Species_fin$DiseaseActivityAtBaseLine[i] = Rem_Species_fin$DiseaseActivityAtBaseLine[i]
}
```

**Wilcoxon test on ratio F.prau/E.coli before and after riboflavin**
```Rem_M1_Ratio = Rem_Species_fin[Rem_Species_fin$time_point=="M1",]
Rem_M3_Ratio = Rem_Species_fin[Rem_Species_fin$time_point=="M3",]

Rem_m1 = Rem_M1_Ratio$Ratio_Prau_Coli
Rem_m3 = Rem_M3_Ratio$Ratio_Prau_Coli

Rem_Rat = data.frame(Y=c(Rem_m1, Rem_m3), Site=factor(rep(c("Rem_m1", "Rem_m3"), times=c(length(Rem_m1), length(Rem_m3)))))
y = Rem_Rat$Y
x = Rem_Rat$Site
wilcox.test(y ~ x, data=Rem_Rat) 

#p=0.3865
mean(Rem_M1_Ratio$Ratio_Prau_Coli)
mean(Rem_M3_Ratio$Ratio_Prau_Coli)
```

