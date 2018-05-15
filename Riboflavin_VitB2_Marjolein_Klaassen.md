**Effects of vitamin B2 (Riboflavin) supplementation on the gut microbiome of CD Patients**

**Authors: Marjolein Klaassen and Ranko Gacesa**

**Date: 05-04-2018** 


Taxonomy
-------------  

Setting working directory
``` 
setwd("~/Documents/IBD Weersma/Vitamin B2/Working directory")
``` 

Importing microbiome taxonomy data
``` 
db_VitB2 = read.csv("metaphlanmerged.txt", header = T, sep = "\t", stringsAsFactors = F)
db_VitB2 = as.data.frame(db_VitB2)
``` 

Making the first column (i.e. microbiome features) the rownames. 
``` 
rownames(db_VitB2) = db_VitB2$ID
``` 
Removing the first column. 
``` 
db_VitB2$ID = NULL
``` 
Only selecting the taxonomical level we are interested in (all). 
The 'filterMetaGenomeDF' function is written by Ranko Gacesa in his R-scripts for microbiome data. Firstly, this has
to run. 
``` 
#_____________________________ R-script part of Ranko Gacesa _____________________________
filterMetaGenomeDF <- function(inDF,presPerc = 0.1,minMRelAb = 0.01,minMedRelAb=0.0, rescaleTaxa=F,verbose=T,
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
``` 

Now, I will subset all the microbiome data in the different taxonomical files.
``` 
tTrans = t(db_VitB2)
db_VitB2_King = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = -1,minMedRelAb = -1,keepLevels = c("K"))
db_VitB2_King = as.data.frame(db_VitB2_King)
db_VitB2_Phylum = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = -1,minMedRelAb = -1,keepLevels = c("P"))
db_VitB2_Phylum = as.data.frame(db_VitB2_Phylum)
db_VitB2_Class = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = -1,minMedRelAb = -1,keepLevels = c("C"))
db_VitB2_Class = as.data.frame(db_VitB2_Class)
db_VitB2_Order = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = -1,minMedRelAb = -1,keepLevels = c("O"))
db_VitB2_Order = as.data.frame(db_VitB2_Order)
db_VitB2_Family = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = -1,minMedRelAb = -1,keepLevels = c("F"))
db_VitB2_Family = as.data.frame(db_VitB2_Family)
db_VitB2_Genus = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = -1,minMedRelAb = -1,keepLevels = c("G"))
db_VitB2_Genus = as.data.frame(db_VitB2_Genus)
db_VitB2_Species = filterMetaGenomeDF(inDF=tTrans,presPerc = 0.10, minMRelAb = -1,minMedRelAb = -1,keepLevels = c("S"))
db_VitB2_Species = as.data.frame(db_VitB2_Species)

# Testing whether data is indeed relative abundance (just to be sure)
rowSums(db_VitB2_King, na.rm = FALSE, dims = 1)
rowSums(db_VitB2_Phylum, na.rm = FALSE, dims = 1)
rowSums(db_VitB2_Class, na.rm = FALSE, dims = 1)
rowSums(db_VitB2_Order, na.rm = FALSE, dims = 1)
rowSums(db_VitB2_Family, na.rm = FALSE, dims = 1)
rowSums(db_VitB2_Genus, na.rm = FALSE, dims = 1)
rowSums(db_VitB2_Species, na.rm = FALSE, dims = 1)
``` 

**SPECIES**
Normalize data (so that it will be normally distributed) via arcsine square root transformation (similar to MaAsLin)
``` 
db_VitB2_Species = t(db_VitB2_Species)
for (c in c(1:nrow(db_VitB2_Species))) {
  db_VitB2_Species[c,] <- asin(sqrt(db_VitB2_Species[c,]))
}
db_VitB2_Species[is.na(db_VitB2_Species)] <- 0.0
``` 

Making t=1 and t=4 time groups 
``` 
t1 = db_VitB2_Species[,c(1, grep("M1", colnames(db_VitB2_Species)))]
rownames(t1) = rownames(db_VitB2_Species)
t1 = as.data.frame(t1)
t4 = db_VitB2_Species[,c(1, grep("M4", colnames(db_VitB2_Species)))]
rownames(t4) = rownames(db_VitB2_Species)
t4 = as.data.frame(t4)
``` 

Removing samples that have not both T1 and T4 (we only want to include samples that have both time points)
These might change, for we have 3 T4's that have been send to the Broad, but have not came back from sequencing.  
``` 
colnames(t1) = gsub("_M1_metaphlan", "", colnames(t1))
colnames(t4) = gsub("_M4_metaphlan", "", colnames(t4))

InBothSamples = intersect(colnames(t1), colnames(t4))
t1 = t1[,InBothSamples]
t4 = t4[,InBothSamples]
```







**Adonis Analyses (taxa)**
I will check if there is a larger proportion of explained variance by M1/M3 difference than by interindividual difference (Pnumber). I will split patients who are active at baseline and patients who are in remission at baseline. 


```
tt1 = as.data.frame(t(t1))
tt1["Pnumber"] = rownames(tt1)
tt1 = tt1[,c(164, 1:163)]
```

```
t0_Clin = merge(db_Clin_Jul, tt1, by="Pnumber", all = FALSE)
t0_Clin = t0_Clin[,c(1:4, 6,7, 9:173)]
```

```
t0_Clin["DiseaseActivityAtBaseLine"] = NA
t0_Clin = t0_Clin[,c(1, 172, 2:171)]

for (i in 1:nrow(t0_Clin)){
  if(t0_Clin$Lab1Calprotectin[i] > 200){
    t0_Clin$DiseaseActivityAtBaseLine[i] = "Active_baseline"
  } else 
    t0_Clin$DiseaseActivityAtBaseLine[i] = "Remission_baseline"
}

t0_Clin["Pnumber_tp"] = t0_Clin$Pnumber
t0_Clin = t0_Clin[,c(173, 1:172)]
t0_Clin$Pnumber_tp = paste0("T0_", t0_Clin$Pnumber_tp)
```

```
tt4 = as.data.frame(t(t4))
tt4["Pnumber"] = rownames(tt4)
tt4 = tt1[,c(164, 1:163)]

t3_Clin = merge (db_Clin_Jul, tt4, by="Pnumber", all= FALSE)
t3_Clin = t3_Clin[,c(1:3, 5, 6, 8:173)]

t3_Clin["Pnumber_tp"] = t3_Clin$Pnumber
t3_Clin = t3_Clin[,c(172, 1:171)]
t3_Clin$Pnumber_tp = paste0("T3_", t3_Clin$Pnumber_tp)
```
```
t3_Clin["DiseaseActivityAtBaseLine"] = NA
t3_Clin = t3_Clin[,c(1, 173, 2:172)]

for (i in 1:nrow(t3_Clin)){
  if(t0_Clin$Lab1Calprotectin[i] > 200){
    t3_Clin$DiseaseActivityAtBaseLine[i] = "Active_baseline"
  } else 
    t3_Clin$DiseaseActivityAtBaseLine[i] = "Remission_baseline"
}
```

Active at baseline n=29, remission at baseline n=38. 

Adding phenotype M1 or M3. 
```
t0_Clin["time_point"] = "M1"
t0_Clin = t0_Clin[,c(1, 174, 2:173)]

t3_Clin["time_point"] = "M3"
t3_Clin = t3_Clin[,c(1, 174, 3, 2, 4:173)]
```

Making similar columnames in T0 and T3. 
```
colnames(t0_Clin)[7] = "LabCalprotectin"
colnames(t3_Clin)[7] = "LabCalprotectin"

colnames(t0_Clin)[9] = "HBI"
colnames(t3_Clin)[9] = "HBI"
```

Binding databases. 
```
Active_Rib_M1 = t0_Clin[t0_Clin$DiseaseActivityAtBaseLine == "Active_baseline",]
ActiveRib_M3_wantM1 = t3_Clin[t3_Clin$DiseaseActivityAtBaseLine == "Active_baseline",]

Rem_Rib_M1 = t0_Clin[t0_Clin$DiseaseActivityAtBaseLine == "Remission_baseline",]
Rem_Rib_M3_wantM1 = t3_Clin[t3_Clin$DiseaseActivityAtBaseLine == "Remission_baseline",]

Act_AdonisRiboSpecies = rbind(Active_Rib_M1, ActiveRib_M3_wantM1)
Rem_AdonisRiboSpecies = rbind(Rem_Rib_M1, Rem_Rib_M3_wantM1)
```

Adonis based on remission
```
Rem_species = Rem_AdonisRiboSpecies[,c(12:174)]

#### phenotype data
Phenotype_dat = Rem_AdonisRiboSpecies[,c(2, 3, 5, 6, 8)]


my_results <- matrix(ncol = 3, nrow=ncol(Phenotype_dat))       

library(vegan)

#For each column in factor_table (factor)  
for (i in 1:ncol(Phenotype_dat)) {
  
  #Create a table for complete cases
  final_factor_table = Phenotype_dat[complete.cases(Phenotype_dat[,i]),]
  
  filter_tax_table = Rem_species[rownames(final_factor_table),]
  
  ad = adonis(formula = filter_tax_table ~ final_factor_table[,i] , data = final_factor_table, permutations = 1000, method = "bray")
  aov_table = ad$aov.tab
  
  my_results[i,1]=aov_table[1,1]
  my_results[i,2]=aov_table[1,5]
  my_results[i,3]=aov_table[1,6]
  
}


rownames(my_results) = colnames(Phenotype_dat)
colnames(my_results) = c("Df", "R2", "Pr(>F)")

Adonis_Remission_Species <- as.data.frame(my_results)

write.table(Adonis_Remission_Species, file="~/remspeciesadonis.txt", quote=F, sep = "\t")



### P-value correction

library(stats)

p_correction_fdr <- as.data.frame(p.adjust(Adonis_Remission_Species$`Pr(>F)`, method = "fdr"))
rownames(p_correction_fdr) <- rownames(Adonis_Remission_Species)

p_correction_bonferroni <- as.data.frame(p.adjust(Adonis_Remission_Species$`Pr(>F)`, method = "bonferroni"))
rownames(p_correction_bonferroni) <- rownames(Adonis_Remission_Species)


#Merge with adonis_results table

final_adonis_results <- merge(Adonis_Remission_Species, p_correction_fdr, by="row.names")
rownames(final_adonis_results) <- final_adonis_results[,1]
final_adonis_results <- final_adonis_results[,-1]

final_adonis_results <- merge(final_adonis_results, p_correction_bonferroni, by="row.names")
rownames(final_adonis_results) <- final_adonis_results[,1]
final_adonis_results <- final_adonis_results[,-1]


colnames(final_adonis_results)[4] <- "FDR_p_value"
colnames(final_adonis_results)[5] <- "Bonferroni_p_value"

write.table(final_adonis_results, file="~/final_adonis_results.txt", sep= "\t", quote = F)
```
