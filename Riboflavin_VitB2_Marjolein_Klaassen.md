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

res =  NULL
``` 

Now, I will compare T1 with T4, using a paired Wilcoxon. This was written together with Ranko Gacesa. 
``` 
for (c in c(1:nrow(t1)) ) {
  #for (c in c(1:10)) {
  taxName <- strsplit(rownames(t1)[c],'\\|')[[1]][length(strsplit(rownames(t1)[c],'\\|')[[1]])]
  print (paste('testing',taxName))
  tt1 <-as.numeric(as.vector(t1[c,2:ncol(t1)]))
  tt4 <- as.numeric(as.vector(t4[c,2:ncol(t4)]))
  tt1df <- as.data.frame(tt1)
  tt1df$Time <- "T1"
  tt1df$PairNR <- c(1:nrow(tt1df))
  colnames(tt1df) <- c("Value","Time","PairNR")
  tt4df <- as.data.frame(tt4)
  tt4df$Time <- "T4"
  tt4df$PairNR <- c(1:nrow(tt4df))
  colnames(tt4df) <- c("Value","Time","PairNR")
  toPlotDF <- rbind.data.frame(tt1df,tt4df)
  wc <- wilcox.test(tt1,tt4,alternative = "two.sided",paired = T)
  print(wc)
  # plot groups
  g <- ggplot(data=toPlotDF,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_boxplot(outlier.alpha = 0.0,width=0.1,alpha=0.5) + 
    geom_jitter(width = 0.25,height = 0.01) + ylab(paste("Rel. Abundance of",taxName)) 
  #ggsave(filename = paste('p_results/plot_',c,'.png', sep=''))
  # plot changes, drop 0s
  toPlotDF2 <- toPlotDF
  g2 <- ggplot(data=toPlotDF2,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_line(aes(group=PairNR),linetype="longdash",col="darkgray") + geom_point()
  #ggsave(g2,filename = paste('p_results/plot_',c,'_','pairs','.png', sep=''))
  
  if (!is.null(res)) {
    res2 <- data.frame("N"=c, "taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
    res <- rbind.data.frame(res2,res)
  } else {
    res <- data.frame("N"=c,"taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
  }
}
# p-value adjust
res$FDR <- p.adjust(res$p.value,method = "fdr")

write.table(res,file = 'VitB2_Results_Species.csv',row.names = F,sep=",")
``` 


**Genus** 
``` 
# Normalize data (so that it will be normally distributed) via arcsine square root transformation (similar to MaAsLin)
db_VitB2_Genus = t(db_VitB2_Genus)
for (c in c(1:nrow(db_VitB2_Genus))) {
  db_VitB2_Genus[c,] <- asin(sqrt(db_VitB2_Genus[c,]))
}
db_VitB2_Genus[is.na(db_VitB2_Genus)] <- 0.0


# Making t=1 and t=4 time groups 
t1 = db_VitB2_Genus[,c(1, grep("M1", colnames(db_VitB2_Genus)))]
rownames(t1) = rownames(db_VitB2_Genus)
t1 = as.data.frame(t1)
t4 = db_VitB2_Genus[,c(1, grep("M4", colnames(db_VitB2_Genus)))]
rownames(t4) = rownames(db_VitB2_Genus)
t4 = as.data.frame(t4)

# Removing samples that have not both T1 and T4 (we only want to include samples that have both time points)
# These might change, for we have 3 T4's that have been send to the Broad, but have not came back from sequencing.  
colnames(t1) = gsub("_M1_metaphlan", "", colnames(t1))
colnames(t4) = gsub("_M4_metaphlan", "", colnames(t4))

InBothSamples = intersect(colnames(t1), colnames(t4))
t1 = t1[,InBothSamples]
t4 = t4[,InBothSamples]


# Filter out species that are present in less than <10% of the samples.
# Q: <10% of which samples? Only T1? 

res =  NULL
# Comparing (wilcoxon, paired): between t1 and t4 (species) 
for (c in c(1:nrow(t1)) ) {
  #for (c in c(1:10)) {
  taxName <- strsplit(rownames(t1)[c],'\\|')[[1]][length(strsplit(rownames(t1)[c],'\\|')[[1]])]
  print (paste('testing',taxName))
  tt1 <-as.numeric(as.vector(t1[c,2:ncol(t1)]))
  tt4 <- as.numeric(as.vector(t4[c,2:ncol(t4)]))
  tt1df <- as.data.frame(tt1)
  tt1df$Time <- "T1"
  tt1df$PairNR <- c(1:nrow(tt1df))
  colnames(tt1df) <- c("Value","Time","PairNR")
  tt4df <- as.data.frame(tt4)
  tt4df$Time <- "T4"
  tt4df$PairNR <- c(1:nrow(tt4df))
  colnames(tt4df) <- c("Value","Time","PairNR")
  toPlotDF <- rbind.data.frame(tt1df,tt4df)
  wc <- wilcox.test(tt1,tt4,alternative = "two.sided",paired = T)
  print(wc)
  # plot groups
  g <- ggplot(data=toPlotDF,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_boxplot(outlier.alpha = 0.0,width=0.1,alpha=0.5) + 
    geom_jitter(width = 0.25,height = 0.01) + ylab(paste("Rel. Abundance of",taxName)) 
  #ggsave(filename = paste('p_results/plot_',c,'.png', sep=''))
  # plot changes, drop 0s
  toPlotDF2 <- toPlotDF
  g2 <- ggplot(data=toPlotDF2,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_line(aes(group=PairNR),linetype="longdash",col="darkgray") + geom_point()
  #ggsave(g2,filename = paste('p_results/plot_',c,'_','pairs','.png', sep=''))
  
  if (!is.null(res)) {
    res2 <- data.frame("N"=c, "taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
    res <- rbind.data.frame(res2,res)
  } else {
    res <- data.frame("N"=c,"taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
  }
}
# p-value adjust
res$FDR <- p.adjust(res$p.value,method = "fdr")

write.table(res,file = 'VitB2_Results_Genus.csv',row.names = F,sep=",")
``` 

**Family**
Normalize data (so that it will be normally distributed) via arcsine square root transformation (similar to MaAsLin)
``` 
db_VitB2_Family = t(db_VitB2_Family)
for (c in c(1:nrow(db_VitB2_Family))) {
  db_VitB2_Family[c,] <- asin(sqrt(db_VitB2_Family[c,]))
}
db_VitB2_Family[is.na(db_VitB2_Family)] <- 0.0
``` 

Making t=1 and t=4 time groups 
``` 
t1 = db_VitB2_Family[,c(1, grep("M1", colnames(db_VitB2_Family)))]
rownames(t1) = rownames(db_VitB2_Family)
t1 = as.data.frame(t1)
t4 = db_VitB2_Family[,c(1, grep("M4", colnames(db_VitB2_Family)))]
rownames(t4) = rownames(db_VitB2_Family)
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

res =  NULL
``` 

Comparing (wilcoxon, paired): between t1 and t4 (species) 
``` 
for (c in c(1:nrow(t1)) ) {
  #for (c in c(1:10)) {
  taxName <- strsplit(rownames(t1)[c],'\\|')[[1]][length(strsplit(rownames(t1)[c],'\\|')[[1]])]
  print (paste('testing',taxName))
  tt1 <-as.numeric(as.vector(t1[c,2:ncol(t1)]))
  tt4 <- as.numeric(as.vector(t4[c,2:ncol(t4)]))
  tt1df <- as.data.frame(tt1)
  tt1df$Time <- "T1"
  tt1df$PairNR <- c(1:nrow(tt1df))
  colnames(tt1df) <- c("Value","Time","PairNR")
  tt4df <- as.data.frame(tt4)
  tt4df$Time <- "T4"
  tt4df$PairNR <- c(1:nrow(tt4df))
  colnames(tt4df) <- c("Value","Time","PairNR")
  toPlotDF <- rbind.data.frame(tt1df,tt4df)
  wc <- wilcox.test(tt1,tt4,alternative = "two.sided",paired = T)
  print(wc)
  # plot groups
  g <- ggplot(data=toPlotDF,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_boxplot(outlier.alpha = 0.0,width=0.1,alpha=0.5) + 
    geom_jitter(width = 0.25,height = 0.01) + ylab(paste("Rel. Abundance of",taxName)) 
  #ggsave(filename = paste('p_results/plot_',c,'.png', sep=''))
  # plot changes, drop 0s
  toPlotDF2 <- toPlotDF
  g2 <- ggplot(data=toPlotDF2,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_line(aes(group=PairNR),linetype="longdash",col="darkgray") + geom_point()
  #ggsave(g2,filename = paste('p_results/plot_',c,'_','pairs','.png', sep=''))
  
  if (!is.null(res)) {
    res2 <- data.frame("N"=c, "taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
    res <- rbind.data.frame(res2,res)
  } else {
    res <- data.frame("N"=c,"taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
  }
}
# p-value adjust
res$FDR <- p.adjust(res$p.value,method = "fdr")

write.table(res,file = 'VitB2_Results_Family.csv',row.names = F,sep=",")
``` 

**Order**
``` 
# Normalize data (so that it will be normally distributed) via arcsine square root transformation (similar to MaAsLin)
db_VitB2_Order = t(db_VitB2_Order)
for (c in c(1:nrow(db_VitB2_Order))) {
  db_VitB2_Order[c,] <- asin(sqrt(db_VitB2_Order[c,]))
}
db_VitB2_Order[is.na(db_VitB2_Order)] <- 0.0


# Making t=1 and t=4 time groups 
t1 = db_VitB2_Order[,c(1, grep("M1", colnames(db_VitB2_Order)))]
rownames(t1) = rownames(db_VitB2_Order)
t1 = as.data.frame(t1)
t4 = db_VitB2_Order[,c(1, grep("M4", colnames(db_VitB2_Order)))]
rownames(t4) = rownames(db_VitB2_Order)
t4 = as.data.frame(t4)

# Removing samples that have not both T1 and T4 (we only want to include samples that have both time points)
# These might change, for we have 3 T4's that have been send to the Broad, but have not came back from sequencing.  
colnames(t1) = gsub("_M1_metaphlan", "", colnames(t1))
colnames(t4) = gsub("_M4_metaphlan", "", colnames(t4))

InBothSamples = intersect(colnames(t1), colnames(t4))
t1 = t1[,InBothSamples]
t4 = t4[,InBothSamples]


# Filter out species that are present in less than <10% of the samples.
# Q: <10% of which samples? Only T1? 

res =  NULL
# Comparing (wilcoxon, paired): between t1 and t4 (species) 
for (c in c(1:nrow(t1)) ) {
  #for (c in c(1:10)) {
  taxName <- strsplit(rownames(t1)[c],'\\|')[[1]][length(strsplit(rownames(t1)[c],'\\|')[[1]])]
  print (paste('testing',taxName))
  tt1 <-as.numeric(as.vector(t1[c,2:ncol(t1)]))
  tt4 <- as.numeric(as.vector(t4[c,2:ncol(t4)]))
  tt1df <- as.data.frame(tt1)
  tt1df$Time <- "T1"
  tt1df$PairNR <- c(1:nrow(tt1df))
  colnames(tt1df) <- c("Value","Time","PairNR")
  tt4df <- as.data.frame(tt4)
  tt4df$Time <- "T4"
  tt4df$PairNR <- c(1:nrow(tt4df))
  colnames(tt4df) <- c("Value","Time","PairNR")
  toPlotDF <- rbind.data.frame(tt1df,tt4df)
  wc <- wilcox.test(tt1,tt4,alternative = "two.sided",paired = T)
  print(wc)
  # plot groups
  g <- ggplot(data=toPlotDF,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_boxplot(outlier.alpha = 0.0,width=0.1,alpha=0.5) + 
    geom_jitter(width = 0.25,height = 0.01) + ylab(paste("Rel. Abundance of",taxName)) 
  #ggsave(filename = paste('p_results/plot_',c,'.png', sep=''))
  # plot changes, drop 0s
  toPlotDF2 <- toPlotDF
  g2 <- ggplot(data=toPlotDF2,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_line(aes(group=PairNR),linetype="longdash",col="darkgray") + geom_point()
  #ggsave(g2,filename = paste('p_results/plot_',c,'_','pairs','.png', sep=''))
  
  if (!is.null(res)) {
    res2 <- data.frame("N"=c, "taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
    res <- rbind.data.frame(res2,res)
  } else {
    res <- data.frame("N"=c,"taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
  }
}
# p-value adjust
res$FDR <- p.adjust(res$p.value,method = "fdr")

write.table(res,file = 'VitB2_Results_Order.csv',row.names = F,sep=",")
``` 


**Class**
``` 
# Normalize data (so that it will be normally distributed) via arcsine square root transformation (similar to MaAsLin)
db_VitB2_Class = t(db_VitB2_Class)
for (c in c(1:nrow(db_VitB2_Class))) {
  db_VitB2_Class[c,] <- asin(sqrt(db_VitB2_Class[c,]))
}
db_VitB2_Class[is.na(db_VitB2_Class)] <- 0.0


# Making t=1 and t=4 time groups 
t1 = db_VitB2_Class[,c(1, grep("M1", colnames(db_VitB2_Class)))]
rownames(t1) = rownames(db_VitB2_Class)
t1 = as.data.frame(t1)
t4 = db_VitB2_Class[,c(1, grep("M4", colnames(db_VitB2_Class)))]
rownames(t4) = rownames(db_VitB2_Class)
t4 = as.data.frame(t4)

# Removing samples that have not both T1 and T4 (we only want to include samples that have both time points)
# These might change, for we have 3 T4's that have been send to the Broad, but have not came back from sequencing.  
colnames(t1) = gsub("_M1_metaphlan", "", colnames(t1))
colnames(t4) = gsub("_M4_metaphlan", "", colnames(t4))

InBothSamples = intersect(colnames(t1), colnames(t4))
t1 = t1[,InBothSamples]
t4 = t4[,InBothSamples]


# Filter out species that are present in less than <10% of the samples.
# Q: <10% of which samples? Only T1? 

res =  NULL
# Comparing (wilcoxon, paired): between t1 and t4 (species) 
for (c in c(1:nrow(t1)) ) {
  #for (c in c(1:10)) {
  taxName <- strsplit(rownames(t1)[c],'\\|')[[1]][length(strsplit(rownames(t1)[c],'\\|')[[1]])]
  print (paste('testing',taxName))
  tt1 <-as.numeric(as.vector(t1[c,2:ncol(t1)]))
  tt4 <- as.numeric(as.vector(t4[c,2:ncol(t4)]))
  tt1df <- as.data.frame(tt1)
  tt1df$Time <- "T1"
  tt1df$PairNR <- c(1:nrow(tt1df))
  colnames(tt1df) <- c("Value","Time","PairNR")
  tt4df <- as.data.frame(tt4)
  tt4df$Time <- "T4"
  tt4df$PairNR <- c(1:nrow(tt4df))
  colnames(tt4df) <- c("Value","Time","PairNR")
  toPlotDF <- rbind.data.frame(tt1df,tt4df)
  wc <- wilcox.test(tt1,tt4,alternative = "two.sided",paired = T)
  print(wc)
  # plot groups
  g <- ggplot(data=toPlotDF,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_boxplot(outlier.alpha = 0.0,width=0.1,alpha=0.5) + 
    geom_jitter(width = 0.25,height = 0.01) + ylab(paste("Rel. Abundance of",taxName)) 
  #ggsave(filename = paste('p_results/plot_',c,'.png', sep=''))
  # plot changes, drop 0s
  toPlotDF2 <- toPlotDF
  g2 <- ggplot(data=toPlotDF2,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_line(aes(group=PairNR),linetype="longdash",col="darkgray") + geom_point()
  #ggsave(g2,filename = paste('p_results/plot_',c,'_','pairs','.png', sep=''))
  
  if (!is.null(res)) {
    res2 <- data.frame("N"=c, "taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
    res <- rbind.data.frame(res2,res)
  } else {
    res <- data.frame("N"=c,"taxon"=taxName,"means.delta"=mean(tt1)-mean(tt4),"p-value"=wc$p.value)
  }
}
# p-value adjust
res$FDR <- p.adjust(res$p.value,method = "fdr")

write.table(res,file = 'VitB2_Results_Class.csv',row.names = F,sep=",")
``` 



Pathways 
------------- 

Setting working directory 
```
setwd("~/Documents/IBD Weersma/Vitamin B2/Working directory")
```

Importing microbiome taxonomy data 
```
db_VitB2_Pathw <- read.csv(file = '_p_pathways_abundances.tsv',sep='\t',stringsAsFactors = F)
db_VitB2_Pathw = as.data.frame(db_VitB2_Pathw)
```
Now, I want to remove all the different levels of pathways.
```
rownames(db_VitB2_Pathw) = db_VitB2_Pathw$Pathway; db_VitB2_Pathw$Pathway = NULL
```
Remove redundant rows
```db_VitB2_Pathw <- db_VitB2_Pathw[grep("__",rownames(db_VitB2_Pathw),invert = T),]
db_VitB2_Pathw <- db_VitB2_Pathw[grep("unclassified",rownames(db_VitB2_Pathw),invert = T),]

db_VitB2_Pathw = t(db_VitB2_Pathw)
db_VitB2_Pathw = as.data.frame(db_VitB2_Pathw)
```

Make it relative abundances 
```
for (rn in c(1:nrow(db_VitB2_Pathw))) {
  db_VitB2_Pathw[rn,] <- db_VitB2_Pathw[rn,]/sum(db_VitB2_Pathw[rn,])
}
```

Below, is code of Ranko Gacesa, in which he filters the pathways for presence, relative abundance, mean and median.
```
filterHumannDF <- function(inDF,presPerc = 0.1,minMRelAb = 0.01,minMedRelAb=0.0,minSum=90.0, rescale=T,verbose=T) {
  tCols = grep('PWY',colnames(inDF)) # colums with pathways
  tColsNMG = grep('PWY',colnames(inDF),invert = T) # colums without pathways
  # replaces NAs with 0s
  for (c in tCols) {
    inDF[,c][is.na(inDF[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    for (r in seq(1,nrow(inDF))) {
      if (verbose) {
        print(r)
      }
      if (sum(inDF[r,tCols]) == 0) {
        print(paste('r=',r,'sum pathways=0'))
        inDF[r,tCols] <- 0.0
      }
      else {
        inDF[r,tCols] <- inDF[r,tCols]/sum(inDF[r,tCols])*100.0
      }
    }
  }
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    nrnZ = as.numeric(sum(inDF[,c]!=0.0))
    #print(nrnZ)
    if (nrnZ/as.numeric(nrow(inDF)) < presPerc) {
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
  tCols = grep('PWY',colnames(inDF)) # colums with pathways
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'pathways!, ',length(tCols),'pathways left!')); }
  
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
  tCols = grep('PWY',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'pathways!, ',length(tCols),'pathways left!')); }
  
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
  tCols = grep('PWY',colnames(inDF)) # colums with pwy
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'pathways!, ',length(tCols),'pathways left!')); }  
  # filter for final sum
  nrRemoved = 0
  toRemove = c()
  #a = c()
  for (r in seq(1,nrow(inDF))) {
    s <- sum(inDF[r,tCols])
    #a = c(a,s)
    #print(s)
    if (s < minSum) {
      print(paste('r=',r,'sum pathways=',s))
      toRemove <- c(toRemove,r)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[-toRemove,]
  }
  if (verbose) {print (paste(' > sum filter: Removed',length(toRemove),'rows!, ',nrow(inDF),'rows left!')); }
  #  print(sorted(a)[1:100])
  inDF
}
```
```
db_VitB2_Pathw = t(db_VitB2_Pathw)
db_VitB2_Pathw = as.data.frame(db_VitB2_Pathw)
```

Perform arc sine square root transformation (identical to what MaAsLin does)
```
for (c in c(1:nrow(db_VitB2_Pathw))) {
  db_VitB2_Pathw[c,] <- asin(sqrt(db_VitB2_Pathw[c,]))
}
db_VitB2_Pathw[is.na.data.frame(db_VitB2_Pathw)] <- 0.0
```

Only keep samples that have T1 and T4. 
```
t1 <- db_VitB2_Pathw[,c(grep("M1",colnames(db_VitB2_Pathw)))]
colnames(t1) <- gsub("_M1","",colnames(t1))
t1 = as.data.frame(t1)
t2 <- db_VitB2_Pathw[,c(grep("M4",colnames(db_VitB2_Pathw)))]
colnames(t2) <- gsub("_M4","",colnames(t2))
t4 = as.data.frame(t2) 

inCommonT1T2 <- intersect(colnames(t1),colnames(t2))
```

**Now, we will perform pairedd Wilcoxon tests on T1 and T4. Furthermore, we will let R create images right away**. 
```
t1c <- t1[,inCommonT1T2]
t2c <- t2[,inCommonT1T2]

# iterate over entries!
res <- NULL
for (c in c(1:nrow(t1c)) ) {
  taxName <- rownames(t1c)[c]
  print (paste(' -> testing',taxName,' NR = ',c))
  tt1df <- as.data.frame(t(t1c[c,]))
  tt1df$Time <- "T1"
  tt1df$PairNR <- c(1:nrow(tt1df))
  colnames(tt1df) <- c("Value","Time","PairNR")
  tt2df <- as.data.frame(t(t2c[c,]))
  tt2df$Time <- "T2"
  tt2df$PairNR <- c(1:nrow(tt2df))
  colnames(tt2df) <- c("Value","Time","PairNR")
  # get rid of all-zero columns <for plotting>
  print (' --> removing all-zero cases')
  getrid <- c()
  for (r in c(1:nrow(tt2df))) {
    if (tt1df$Value[r] == 0 & tt2df$Value[r] == 0) {
      getrid <- c(getrid,r)
    }
  }
  if (!is.null(getrid)) {
    tt2df <- tt2df[-getrid,]
    tt1df <- tt1df[-getrid,]
  }
  print (nrow(tt1df))
  if (nrow(tt1df) == 0) {
    print ('  WARNING: no cases after removal of 0s')
    next
  }
  print ('  --> doing tests')
  # do tests
  wcT1T2 <- wilcox.test(tt1df$Value,tt2df$Value,alternative = "two.sided",paired = T)
  #wcT1T2 <- t.test(tt1df$Value,tt2df$Value,alternative = "two.sided",paired = T)
  
  print ('     --> tests done')
  # save result
  if (!is.null(res)) {
    res2 <- data.frame("N"=c, "taxon"=taxName,pValue=wcT1T2$p.value)
    #if (is.nan(wc$p.value)) {res2$p.value = 1.0}
    res <- rbind.data.frame(res2,res)
  } else {
    res <- data.frame("N"=c, "taxon"=taxName,pValue=wcT1T2$p.value)
  }
}
res[is.nan.data.frame(res)] <- 1.0
resPA <- res
resPA$FDR <- p.adjust(resPA$pValue,method = "fdr")
#resPA$N <- NULL
resPA <- resPA[order(resPA$pValue),]
#resPA$N <- c(1:nrow(resPA))
write.table(resPA,file = paste('PWY_results_taxa.csv',sep=''),row.names = F,sep=",")

# now plot, sort by FDR
plotNR <- 0
for (c in resPA$N)
{
  plotNR <- plotNR + 1
  taxName <- rownames(t1c)[c]
  print (paste(' -> plotting',taxName,' NR = ',c))
  tt1df <- as.data.frame(t(t1c[c,]))
  tt1df$Time <- "T1"
  tt1df$PairNR <- c(1:nrow(tt1df))
  colnames(tt1df) <- c("Value","Time","PairNR")
  tt2df <- as.data.frame(t(t2c[c,]))
  tt2df$Time <- "T2"
  tt2df$PairNR <- c(1:nrow(tt2df))
  colnames(tt2df) <- c("Value","Time","PairNR")
  
  # get rid of all-zero columns <for plotting>
  print (' --> removing all-zero cases')
  getrid <- c()
  for (r in c(1:nrow(tt1df))) {
    if (tt1df$Value[r] == 0 & tt2df$Value[r] == 0) {
      getrid <- c(getrid,r)
    }
  }
  if (!is.null(getrid)) {
    tt2df <- tt2df[-getrid,]
    tt1df <- tt1df[-getrid,]
  }
  print (nrow(tt1df))
  if (nrow(tt1df) == 0) {
    print ('  WARNING: no cases after removal of 0s')
    next
  }
  toPlotDF <- rbind.data.frame(tt1df,tt2df)
  ypos <- max(toPlotDF$Value)*(1.05)
  pv <- resPA$pValue[plotNR]
  fdr <- resPA$FDR[plotNR]
  lblT1T2 <- paste("Pv=",formatC(pv, format = "g", digits = 2),"; FDR=",formatC(fdr, format = "g", digits = 2),sep="")
  if (as.numeric(fdr) <= 0.1) {lblT1T2 <- paste("*",lblT1T2)}
  if (as.numeric(fdr) <= 1.0e-5) {lblT1T2 <- paste("**",lblT1T2,sep="")}
  
  print ('  --> plotting')
  g <- ggplot(data=toPlotDF,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + geom_boxplot(outlier.alpha = 0.0,width=0.1,alpha=0.5) + 
    geom_jitter(width = 0.25,height = 0.01) + ylab(paste("Rel. Abundance of",taxName)) 
  g <- g + annotate("text", x=1.4,y=ypos, label=lblT1T2,vjust=0 )
  g <- g + annotate("segment", x = 1+0.01, y = ypos*1*0.98, xend = 2-0.01,lineend = "round",
                    yend = ypos*1*0.98,size=0.5,colour="black")
  ggsave(g,filename = paste('Vitb2_pathways',plotNR,'.png', sep=''))
  
  g2 <- ggplot(data=toPlotDF,aes(y=Value,x=Time,col=Time)) + geom_violin(alpha=0.5) + 
    geom_line(aes(group=PairNR),linetype="longdash",col="darkgray") + geom_point() + ylab(paste("Rel. Abundance of",taxName)) 
  g2 <- g2 + annotate("text", x=1.4,y=ypos, label=lblT1T2,vjust=0 )
  g2 <- g2 + annotate("segment", x = 1+0.01, y = ypos*1*0.98, xend = 2-0.01,lineend = "round",
                      yend = ypos*1*0.98,size=0.5,colour="black")
  ggsave(g2,filename = paste('VitB2_pathways',plotNR,'_','pairs','.png', sep=''))
  print ('    --> plotting done')
}
```
