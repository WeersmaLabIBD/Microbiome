Healthy cotwins share gut microbiome signatures with their inflammatory bowel disease twins and unrelated patients
--------------

**Author of code: Marjolein A.Y. Klaassen, Ranko Gacesa, Eelco C. Brand**
**Year: 2020**
**Aim: Study gut microbiomes of IBD-discordant and -concordant twin pairs, which offers the unique opportunity to assess individuals at risk of developing IBD, namely healthy cotwins from IBD-discordant twin pairs.**

--------------

**1a. Filtering steps (taxa)** 

```
# ================================================================================
# ================================================================================
# function for filtering microbiome (metaphlan) results
# - takes dataframe (any)
# - replaces NAs with 0.0
# - filters OUT taxa present in less then presPerc samples
# - filters OUT taxa with mean releative abundance < minMRelAb
# - filters OUT taxa with median relative abundance < minMedRelAb
# - filters OUT taxonomic levels not in keepLevels ()
#   -> keepLevels should be vector of following: T = strain, S = species, G = Genera, F = families
#                                                O = Orders, C = classes, P = phyla, K = Kingdoms
#   -> example: keepLevels=c('S','G') keeps only species and genera
# - filters OUT domains not in keepDomains: 
#   -> input is vector of following: Eukaryota, Bacteria, Viruses, Archaea
#   -> example: keepDomains=c('B') keeps only bacteria, removes viruses, archaea and eukarya
# - if rescaleTaxa = True, rescales relative abundancies after filtering
# returns modified dataframe
# NOTES:
# - assumes metaphlan encoding (k__<kingdom>.o__<order> ... ); it will mess stuff
# up if non-metagenome rows are encoded like this!
# DEFAULTS: 
# - removes non-bacteria
# - keeps all except kingdoms
# ================================================================================

#TODO: implement keepDomains
filterMetaGenomeDF <- function(inDF,presPerc = 0.1,minMRelAb = 0.01,minMedRelAb=0.0, rescaleTaxa=F,verbose=T,
                               keepDomains=c('Bacteria','Archaea'),
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
      if (verbose) {
        print (paste('col',c,': ',colnames(inDF)[c],'; nr non Zero:',nrnZ,'=',nrnZ/as.numeric(nrow(inDF)),'>> Column removed!'))
      }
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('^[dgtspcfko]__',colnames(inDF)) # colums with microbiome
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
    if (rescaleTaxa) {taxaK <- taxaK/rowSums(taxaK)}
    oDF <- cbind(oDF,taxaK)}
  if ('P' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaP),'Phyla'))} 
    if (rescaleTaxa) {taxaP <- taxaP/rowSums(taxaP)}
    oDF <- cbind(oDF,taxaP)}
  if ('C' %in% keepLevels) {
    if (verbose) {print(paste(' -> kept',ncol(taxaC),'Classes'))} 
    if (rescaleTaxa) {taxaC <- taxaC/rowSums(taxaC)}
    oDF <- cbind(oDF,taxaC)}
  if ('O' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaO),'Orders'))}
    if (rescaleTaxa) {taxaO <- taxaO/rowSums(taxaO)}
    oDF <- cbind(oDF,taxaO)}
  if ('F' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaF),'Families'))}
    if (rescaleTaxa) {taxaF <- taxaF/rowSums(taxaF)}
    oDF <- cbind(oDF,taxaF)}
  if ('G' %in% keepLevels) {if (verbose){ print(paste(' -> kept',ncol(taxaG),'Genera'))}
    if (rescaleTaxa) {taxaG <- taxaG/rowSums(taxaG)}
    oDF <- cbind(oDF,taxaG)}
  if ('S' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaS),'Species'))}
    if (rescaleTaxa) {taxaS <- taxaS/rowSums(taxaS)}
    oDF <- cbind(oDF,taxaS)}
  if ('T' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaT),'Strains'))}
    if (rescaleTaxa) {taxaT <- taxaT/rowSums(taxaT)}
    oDF <- cbind(oDF,taxaT)}
  
  if (verbose) {print ('data processing done, returning Dataframe')}
  oDF
}
# ================================================================================

```

**1b. Filtering steps (pathways)** 

```

filterHumannDF <- function(inDF,presPerc = 0.05,minMRelAb = 0.001,minMedRelAb=0.0,minSum=90.0, rescale=T,verbose=T,type='MetaCyc') {
  
  nonPWYpwys <- c("ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis",
                  "CHLOROPHYLL-SYN: chlorophyllide a biosynthesis I (aerobic, light-dependent)",
                  "GLYCOLYSIS-E-D: superpathway of glycolysis and Entner-Doudoroff",
                  "GLYCOLYSIS-TCA-GLYOX-BYPASS: superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
                  "GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)",
                  "GLYOXYLATE-BYPASS: glyoxylate cycle",
                  "HEME-BIOSYNTHESIS-II: heme biosynthesis I (aerobic)",
                  "MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS: protein N-glycosylation (eukaryotic, high mannose)",
                  "NAD-BIOSYNTHESIS-II: NAD salvage pathway II",                  
                  "REDCITCYC: TCA cycle VIII (helicobacter)",
                  "TCA-GLYOX-BYPASS: superpathway of glyoxylate bypass and TCA",
                  "TCA: TCA cycle I (prokaryotic)")
  
  colnames(inDF)[colnames(inDF) %in% nonPWYpwys] <- paste0('PWY_',colnames(inDF)[colnames(inDF) %in% nonPWYpwys])
  
  if (type=='MetaCyc') {
    nonPWYdf <- as.data.frame(inDF[,-grep('PWY',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^EC_',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    nonPWYdf <- as.data.frame(inDF[,-grep('RXN',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^PF[01]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^GO',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^K[012]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^K[012]',colnames(inDF))] ])
  }
  colnames(nonPWYdf) <- cnsNonPWYdf
  if (type=='MetaCyc') {
    yesPWYdf <- as.data.frame(inDF[,grep('PWY',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    yesPWYdf <- as.data.frame(inDF[,grep('^EC_',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    yesPWYdf <- as.data.frame(inDF[,grep('RXN',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    yesPWYdf <- as.data.frame(inDF[,grep('^PF[01]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    yesPWYdf <- as.data.frame(inDF[,grep('^GO',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    yesPWYdf <- as.data.frame(inDF[,grep('^K[012]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^K[012]',colnames(inDF))] ])
  }
  
  # replaces NAs with 0s
  for (c in colnames(yesPWYdf)) {
    yesPWYdf[,c][is.na(yesPWYdf[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    nrnZ = as.numeric(sum(yesPWYdf[,c]!=0.0))
    if (nrnZ/as.numeric(nrow(yesPWYdf)) < presPerc) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = mean(yesPWYdf[,c])
    if ( mn < minMRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = median(yesPWYdf[,c])
    if ( mn < minMedRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # do final rescale
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  inDF <- cbind.data.frame(nonPWYdf,yesPWYdf)
  if (verbose) {print ('> DONE')}
  inDF
}

```

**2. PCoA analyses**

```

library("vegan")
Beta = vegdist(taxa_file, pwy_file, method="bray")
PCoA=cmdscale(Beta, k = 5) # R uses cmdscale() to calculate classical multi-dimensional scaling, a synonym for principal coordinates analysis.

ggplot(PCoA_taxa, PCoA_pwy, aes(x=V1,y=V2, geom= "blank", color=df$ColorsDots)) + geom_point()  + ggtitle("") + theme_classic() + geom_point (data=centroids2, size=8, alpha=0.9, col=centroids2$colors) + labs(x="PCoA1", y="PCoA2") +  scale_color_identity ("", breaks=c("green3", "black", "blue2", "orangered"), labels=c("healthy controls (n=495)", "healthy twins (n=38)", "IBD twins (n=61)", "IBD controls (n=99)"), guide="legend") 

## p-values between Bray-Curtis distances per group 

Microbiome_data_frame = dataframe[,c()]
Phenotype_data_frame = dataframe[,c()]
#
my_results <- matrix(ncol = 3, nrow=ncol(Phenotype_data_frame))    
my_results <- matrix(ncol = 3, nrow=ncol(Phenotype_data_frame))       
#
library(vegan)
#For each column in Phenotype_data_frame (factor)  
for (i in 1:ncol(Phenotype_dat)) {
  #Create a table for complete cases
  final_factor_table <- Phenotype_data_frame[complete.cases(Phenotype_data_frame[,i]),]
  filter_tax_table <- Microbiome_data_frame[rownames(final_factor_table),]
  ad <- adonis(formula = filter_tax_table ~ final_factor_table[,i] , data = final_factor_table, permutations = 10000, method = "bray")
  aov_table <- ad$aov.tab
  my_results[i,1]=aov_table[1,1]
  my_results[i,2]=aov_table[1,5]
  my_results[i,3]=aov_table[1,6]
}
#
rownames(my_results) = colnames(Phenotype_dat)
colnames(my_results) = c("Df", "R2", "Pr(>F)")
#
Adonis_my_results <- as.data.frame(my_results)

# P-value correction
library(stats)
p_correction_fdr <- as.data.frame(p.adjust(Adonis_my_results$`Pr(>F)`, method = "fdr"))
rownames(p_correction_fdr) <- rownames(Adonis_my_results)
p_correction_bonferroni <- as.data.frame(p.adjust(Adonis_my_results$`Pr(>F)`, method = "bonferroni"))
rownames(p_correction_bonferroni) <- rownames(Adonis_my_results)
#Merge with adonis_results table
final_adonis_results <- merge(Adonis_my_results, p_correction_fdr, by="row.names")
rownames(final_adonis_results) <- final_adonis_results[,1]
final_adonis_results <- final_adonis_results[,-1]
final_adonis_results <- merge(final_adonis_results, p_correction_bonferroni, by="row.names")
rownames(final_adonis_results) <- final_adonis_results[,1]
final_adonis_results <- final_adonis_results[,-1]
#
colnames(final_adonis_results)[4] <- "FDR_p_value"
colnames(final_adonis_results)[5] <- "Bonferroni_p_value"
#

```

**3.Similarity in gut microbiome composition**

```

## When looking between twin pairs (intertwinpairs), we firstly had to create pairs that were not allowed to originate from the same twin pair

library("vegan")
Beta = vegdist(db_twins, method="bray")
betaMat <- as.matrix(Beta)
PCoA=cmdscale(Beta, k = 5) # R uses cmdscale() to calculate classical multi-dimensional scaling, a synonym for principal coordinates analysis.
PCoA_df = as.data.frame(PCoA)

resDF <- NULL
for (tID in unique(c(PCoA_df$Twin_pair_number)) ) {
  # iterate over twins
  print(tID)
  babyDF <- PCoA_df[PCoA_df$Twin_pair_number==tID, ]
  # test if twin is actually two samples
  if (nrow(babyDF) != 2) {
    print(paste0('WARNING:',tID,' has ',nrow(babyDF),' samples!'))
  } else {
    # calculate distance between twins
    twinDist <- betaMat[ babyDF$Record.Id[1],babyDF$Record.Id[2]  ]
    # make one result
    # collect disease type
    # NOTE: DOES NOT ACCOUNT FOR HEALTHY-HEALTHY PAIRS!!!!
    diseaseTypeL <- unique(babyDF$Disease_status)
    if ("Ulcerative Colitis" %in% diseaseTypeL) {
      disT <- "UC"
    } else if ("IBD unspecified" %in% diseaseTypeL) {
      disT <- "IBDU"
    } else {
      disT <- "CD"
    }
    oneResult <- data.frame(ID1=babyDF$Record.Id[1],
                            ID2=babyDF$Record.Id[2],
                            twinConcordancy=babyDF$Concordancy,
                            twinZigosity=babyDF$zyg_def,
                            diseaseType=disT,
                            bcPairwise=twinDist)
    # collect one result
    resDF <- rbind.data.frame(resDF,oneResult)
  }
}

# When creating random pairs within the healthy controls, IBD controls and unrelated twin pairs, the following code was used:
library("vegan")
Beta = vegdist(db_Twin_Beta, method="bray")
betaMat <- as.matrix(Beta)
library(reshape2)
betaMat2 <- betaMat
betaMat2[lower.tri(betaMat2,diag = T)] <- 0
flatDM <- subset(melt(betaMat2),value !=0)

randomPairs <- sample(x=row.names(flatDM),size = 1000,replace = F)


```


**4. MaAsLin2 analyses**

```

### a. Comparing gut microbiomes between healthy twins and their IBD cotwins 
library(Maaslin2)
Maaslin2(taxa_maas/pwy_maas, meta_taxa_maas, "MaAsLin a", fixed_effects = c("cohort","zygosity","IBD_type_chart_review", "Signs_of_disease_activity", "Sex", "Age_years_calculated", "Antibiotics_past_3_months", "PPI_current", "BMI_survey", "DiseaseLocation_UCCD"), random_effects = c("Twin_pair_number"), normalization = "NONE",transform = "AST", min_prevalence = 0.25)

### b. Comparing gut microbiomes between IBD twins and unrelated IBD controls 
library(Maaslin2)
Maaslin2(taxa_maas/pwy_maas, meta_taxa_maas, "MaAsLin b", fixed_effects = c("cohort","zygosity","IBD_type_chart_review", "Sex", "Age_years_calculated", "Antibiotics_past_3_months", "PPI_current", "BMI_survey", "DiseaseLocation_UCCD"), normalization = "NONE",transform = "AST", min_prevalence = 0.25)

### c. Comparing gut microbiomes between IBD twins and unrelated healthy controls
library(Maaslin2)
Maaslin2(taxa_maas/pwy_maas, meta_taxa_maas, "MaAsLin c", fixed_effects = c("cohort","zygosity","IBD_type_chart_review", "Sex", "Age_years_calculated", "Antibiotics_past_3_months", "PPI_current", "BMI_survey", "DiseaseLocation_UCCD"), normalization = "NONE",transform = "AST", min_prevalence = 0.25)

### d. Comparing gut microbiomes between healthy twins and unrelated healthy controls
library(Maaslin2)
Maaslin2(taxa_maas/pwy_maas, meta_taxa_maas, "MaAsLin d", fixed_effects = c("cohort","zygosity","IBD_type_chart_review", "Sex", "Age_years_calculated", "Antibiotics_past_3_months", "PPI_current", "BMI_survey", "DiseaseLocation_UCCD"), normalization = "NONE",transform = "AST", min_prevalence = 0.25)

### e. Comparing gut microbiomes between healthy controls and unrelated healthy controls
library(Maaslin2)
Maaslin2(taxa_maas/pwy_maas, meta_taxa_maas, "MaAsLin d", fixed_effects = c("cohort","IBD_type_chart_review", "Sex", "Age_years_calculated", "Antibiotics_past_3_months", "PPI_current", "BMI_survey", "DiseaseLocation_UCCD"), normalization = "NONE",transform = "AST", min_prevalence = 0.25)

```
