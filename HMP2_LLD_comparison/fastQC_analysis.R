# -----------------------------------------------------
# Analysis of HMP2 data and comparison to our DAG3 Novogene Pilot samples
# -----------------------------------------------------

# load libs & setwd
# note: uses gdata for excel files, fastqcr for fastQC files
setwd("/home/rgacesa/UMCG/Projects/HMP2_IBD/Ranalysis/")
library(fastqcr)
library(gdata)
library(ggplot2)
library(tidyr)
library(stringi)
library(grid)
library(vegan)
library(gplots)
require(VennDiagram)
# =================================================
# -> HELPER FUNCTIONS
# =================================================
nzmean <- function(a){
  mean(a[a!=0])
}

##Function to calculate nº of 0
zsum <- function(a){
  sum (a==0)
}

##Function to calculate nº of non-0
nsum <- function(a){
  sum (a!=0)
}
# ==================================================

# function that finds overlapping species in two 
# metaphlan files
metaPhlanGetSpecies <- function(inFileA, Ttresh=0.001,taxa="S") {
  # -- process 1st input (A)
  # ----------------------------
  taxaA=read.table(inFileA, sep = "\t", row.names = 1, header = T, check.names = F,blank.lines.skip = T)
  taxaA=as.data.frame(t(taxaA))
  # remove strains
  s_speciesA=taxaA[,grep("t__",colnames(taxaA),invert = TRUE)]
  # remove genera
  s_generaA = s_speciesA[,grep('s__',colnames(s_speciesA),invert = TRUE)]
  # remove families
  s_familyA = s_generaA[,grep('g__',colnames(s_generaA),invert = TRUE)]
  # remove orders
  s_orderA = s_familyA[,grep('f__',colnames(s_familyA),invert = TRUE)]
  # remove classes
  s_classA = s_orderA[,grep('o__',colnames(s_orderA),invert = TRUE)]
  #remove phyla
  s_phylaA = s_classA[,grep('c__',colnames(s_classA),invert = TRUE)]
  # grab species and drop 0s
  s_speciesA=s_speciesA[,grep("s__",colnames(s_speciesA))]
  s_speciesAFN=colnames(s_speciesA[,apply(s_speciesA,2,mean) > Ttresh])
  # grab genera and drop 0s
  s_generaA=s_generaA[,grep("g__",colnames(s_generaA))]
  s_generaAFN=colnames(s_generaA[,apply(s_generaA,2,mean) > Ttresh])
  # grab families and drop 0s
  s_familyA=s_familyA[,grep("f__",colnames(s_familyA))]
  s_familyAFN=colnames(s_familyA[,apply(s_familyA,2,mean) > Ttresh])
  # grab orders and drop 0s
  s_orderA=s_orderA[,grep("o__",colnames(s_orderA))]
  # grab classes and drop 0s
  s_classA=s_classA[,grep("c__",colnames(s_classA))]
  # grab phyla and drop 0s
  s_phylaA=s_phylaA[,grep("p__",colnames(s_phylaA))]
  
  s_speciesAFN
  
}


metaPhlanGetOverlapSp <- function(inFileA,inFileB,Ttresh = 0.001,taxa="S") {
  # -- process 1st input (A)
  # ----------------------------
  taxaA=read.table(inFileA, sep = "\t", row.names = 1, header = T, check.names = F,blank.lines.skip = T)
  taxaA=as.data.frame(t(taxaA))
  # remove strains
  s_speciesA=taxaA[,grep("t__",colnames(taxaA),invert = TRUE)]
  # remove genera
  s_generaA = s_speciesA[,grep('s__',colnames(s_speciesA),invert = TRUE)]
  # remove families
  s_familyA = s_generaA[,grep('g__',colnames(s_generaA),invert = TRUE)]
  # remove orders
  s_orderA = s_familyA[,grep('f__',colnames(s_familyA),invert = TRUE)]
  # remove classes
  s_classA = s_orderA[,grep('o__',colnames(s_orderA),invert = TRUE)]
  #remove phyla
  s_phylaA = s_classA[,grep('c__',colnames(s_classA),invert = TRUE)]
  # grab species and drop 0s
  s_speciesA=s_speciesA[,grep("s__",colnames(s_speciesA))]
  s_speciesAFN=colnames(s_speciesA[,apply(s_speciesA,2,mean) > Ttresh])
  # grab genera and drop 0s
  s_generaA=s_generaA[,grep("g__",colnames(s_generaA))]
  s_generaAFN=colnames(s_generaA[,apply(s_generaA,2,mean) > Ttresh])
  # grab families and drop 0s
  s_familyA=s_familyA[,grep("f__",colnames(s_familyA))]
  s_familyAFN=colnames(s_familyA[,apply(s_familyA,2,mean) > Ttresh])
  # grab orders and drop 0s
  s_orderA=s_orderA[,grep("o__",colnames(s_orderA))]
  # grab classes and drop 0s
  s_classA=s_classA[,grep("c__",colnames(s_classA))]
  # grab phyla and drop 0s
  s_phylaA=s_phylaA[,grep("p__",colnames(s_phylaA))]
  
  # -- process 2nd input 
  # ----------------------------
  taxaB=read.table(inFileB, sep = "\t", row.names = 1, header = T, check.names = F,blank.lines.skip = T)
  taxaB=as.data.frame(t(taxaB))
  # remove strains
  s_speciesB=taxaB[,grep("t__",colnames(taxaB),invert = TRUE)]
  # remove genera
  s_generaB = s_speciesB[,grep('s__',colnames(s_speciesB),invert = TRUE)]
  # remove families
  s_familyB = s_generaB[,grep('g__',colnames(s_generaB),invert = TRUE)]
  # remove orders
  s_orderB = s_familyB[,grep('f__',colnames(s_familyB),invert = TRUE)]
  # remove classes
  s_classB = s_orderB[,grep('o__',colnames(s_orderB),invert = TRUE)]
  #remove phyla
  s_phylaB = s_classB[,grep('c__',colnames(s_classB),invert = TRUE)]
  # grab species and drop 0s
  s_speciesB=s_speciesB[,grep("s__",colnames(s_speciesB))]
  s_speciesBFN=colnames(s_speciesB[,apply(s_speciesB,2,mean) > Ttresh])
  # grab genera and drop 0s
  s_generaB=s_generaB[,grep("g__",colnames(s_generaB))]
  s_generaBFN=colnames(s_generaB[,apply(s_generaB,2,mean) > Ttresh])
  # grab families and drop 0s
  s_familyB=s_familyB[,grep("f__",colnames(s_familyB))]
  s_familyBFN=colnames(s_familyB[,apply(s_familyB,2,mean) > Ttresh])
  # grab orders and drop 0s
  s_orderB=s_orderB[,grep("o__",colnames(s_orderB))]
  # grab classes and drop 0s
  s_classB=s_classB[,grep("c__",colnames(s_classB))]
  # grab phyla and drop 0s
  s_phylaB=s_phylaB[,grep("p__",colnames(s_phylaB))]  
 
  venn(list(inFileA=s_speciesAFN,inFileB=s_speciesBFN),names = c(inFileA,inFileB))
}

  



# function to prep taxonomy data from metaphlan file
# - grab MetaPhlan .txt file, returns dataframe with
# numbers of taxonomical units and diversity measures
# (shannon, Bray-Curtis PC 1-5)
metaPhlanGetStats <- function(inFile) {
  taxa_pilot=read.table(inFile, sep = "\t", row.names = 1, header = T, check.names = F,blank.lines.skip = T)
  taxa_pilot=as.data.frame(t(taxa_pilot))
    # remove strains
  s_species=taxa_pilot[,grep("t__",colnames(taxa_pilot),invert = TRUE)]
  # remove genera
  s_genera = s_species[,grep('s__',colnames(s_species),invert = TRUE)]
  # remove families
  s_family = s_genera[,grep('g__',colnames(s_genera),invert = TRUE)]
  # remove orders
  s_order = s_family[,grep('f__',colnames(s_family),invert = TRUE)]
  # remove classes
  s_class = s_order[,grep('o__',colnames(s_order),invert = TRUE)]
  #remove phyla
  s_phyla = s_class[,grep('c__',colnames(s_class),invert = TRUE)]
  
  # grab species
  s_species=s_species[,grep("s__",colnames(s_species))]
  # grab genera
  s_genera=s_genera[,grep("g__",colnames(s_genera))]
  # grab families
  s_family=s_family[,grep("f__",colnames(s_family))]
  # grab orders
  s_order=s_order[,grep("o__",colnames(s_order))]
  # grab classes
  s_class=s_class[,grep("c__",colnames(s_class))]
  # grab phyla
  s_phyla=s_phyla[,grep("p__",colnames(s_phyla))]
  # count taxonomy stuff per sample
  # ---------------------------------------------------
  s_taxnr <- data.frame(cbind(row.names(s_species), apply(X = s_species,MARGIN = 1,FUN = nsum), 
                              apply(X = s_genera,MARGIN = 1,FUN = nsum),
                              apply(X = s_family,MARGIN = 1,FUN = nsum),
                              apply(X = s_order,MARGIN = 1,FUN = nsum),
                              apply(X = s_class,MARGIN = 1,FUN = nsum),
                              apply(X = s_phyla,MARGIN = 1,FUN = nsum)))
  
  # get rid of rows with no species detected
  s_taxnr <- s_taxnr[apply(s_species,MARGIN = 1,sum) != 0,]
  s_species <- s_species[apply(s_species,MARGIN = 1,sum) != 0,]
  
  # add diversity measures
  s_taxnr <- cbind(s_taxnr,diversity(s_species,index="shannon"))
  colnames(s_taxnr) <- c("sample","NrSpecies","NrGenera","NrFamilies","NrOrders","NrClasses","NrPhyla","Shannon")
  s_taxnr$NrSpecies <- as.numeric(as.character(s_taxnr$NrSpecies))
  s_taxnr$NrGenera <- as.numeric(as.character(s_taxnr$NrGenera))
  s_taxnr$NrFamilies <- as.numeric(as.character(s_taxnr$NrFamilies))
  s_taxnr$NrClasses <- as.numeric(as.character(s_taxnr$NrClasses))
  s_taxnr$NrOrders <- as.numeric(as.character(s_taxnr$NrOrders))
  s_taxnr$NrPhyla <- as.numeric(as.character(s_taxnr$NrPhyla))
  s_taxnr$sample <- rownames(s_taxnr)
  # -> get B-C distance
  beta <- vegdist((s_species), method="bray")
  # -> do PCA
  my_pcoa <- as.data.frame(cmdscale(beta, k = 5))
  # add it to bit dataset
  my_pcoa$sample <- row.names(my_pcoa)
  colnames (my_pcoa) <- c("Bray_Curtis_PC1","Bray_Curtis_PC2","Bray_Curtis_PC3","Bray_Curtis_PC4","Bray_Curtis_PC5","sample")
  #my_pcoa$sample <- stri_sub(my_pcoa$sample,2,-1)
  ret <- merge(my_pcoa,s_taxnr)
  ret
}

# ======================================================================
# prepare taxonomy data (DAG3)
# ======================================================================
dagTax <-metaPhlanGetStats('./Pilot_metaphlan.txt')
dagTax$cohort <- "DAG3_Pilot"

# ======================================================================
# prepare taxonomy data (HMP2)
# ======================================================================
hmp1508Tax <- metaPhlanGetStats('./HMP2_IBD_1508_taxprofile.tsv')
hmp1508Tax$cohort <- "HMP1508"
# 316 samples

hmp1729Tax <- metaPhlanGetStats('./HMP2_IBD_1729_taxprofile.tsv')
hmp1729Tax$cohort <- "HMP1729"
# 1338 samples

# ======================================================================
# prepare taxonomy data (LLD)
# ======================================================================
lldTax <- metaPhlanGetStats('./LLD_taxonomy_metaphlan2_092017.txt')
lldTax$cohort <- "LLD"


# join datasets
alltax <- rbind(lldTax,dagTax,hmp1508Tax,hmp1729Tax)
hmpalltax <- rbind(hmp1508Tax,hmp1729Tax)

# box plots: NR species & other parameters VS novogene diagnosis
ggplot(data=alltax,aes(x=cohort,y=NrSpecies,color=cohort)) + ggtitle("") + geom_jitter(width = 0.2,height = 0.2) + geom_boxplot(alpha=0.4,outlier.alpha = 0.0)
ggplot(data=alltax,aes(x=cohort,y=NrGenera,color=cohort)) + ggtitle("") + geom_jitter(width = 0.2,height = 0.2) + geom_boxplot(alpha=0.4,outlier.alpha = 0.0)
ggplot(data=alltax,aes(x=cohort,y=NrFamilies,color=cohort)) + ggtitle("") + geom_jitter(width = 0.2,height = 0.2) + geom_boxplot(alpha=0.4,outlier.alpha = 0.0)
ggplot(data=alltax,aes(x=cohort,y=NrClasses,color=cohort)) + ggtitle("") + geom_jitter(width = 0.2,height = 0.2) + geom_boxplot(alpha=0.4,outlier.alpha = 0.0)
ggplot(data=alltax,aes(x=cohort,y=NrOrders,color=cohort)) + ggtitle("") + geom_jitter(width = 0.2,height = 0.2) + geom_boxplot(alpha=0.4,outlier.alpha = 0.0)
ggplot(data=alltax,aes(x=cohort,y=NrPhyla,color=cohort)) + ggtitle("") + geom_jitter(width = 0.2,height = 0.2) + geom_boxplot(alpha=0.4,outlier.alpha = 0.0)
ggplot(data=alltax,aes(x=cohort,y=Shannon,color=cohort)) + ggtitle("") + geom_jitter(width = 0.2,height = 0.2) + geom_boxplot(alpha=0.4,outlier.alpha = 0.0)

# --------------- BRAY CURTIS!!!
ggplot(data=alltax,aes(x=Bray_Curtis_PC1,y=Bray_Curtis_PC2,color=cohort)) + geom_point(size=2) + ggtitle("Bray-Curtis PCA analysis (PC2 vs PC1)")
ggplot(data=hmpalltax,aes(x=Bray_Curtis_PC1,y=Bray_Curtis_PC2,color=cohort)) + geom_point(size=2) + ggtitle("Bray-Curtis PCA analysis (PC2 vs PC1)")
hmpalltaxjoined <- hmpalltax
hmpalltaxjoined$cohort <- "HMP2"
daghmp2 <- rbind(hmpalltaxjoined,dagTax)
ggplot(data=daghmp2,aes(x=Bray_Curtis_PC1,y=Bray_Curtis_PC2,color=cohort)) + geom_point(size=2) + ggtitle("Bray-Curtis PCA analysis (PC2 vs PC1)")
daghmp1508 <- rbind(hmp1508Tax,dagTax)
ggplot(data=daghmp1508,aes(x=Bray_Curtis_PC1,y=Bray_Curtis_PC2,color=cohort)) + geom_point(size=2) + ggtitle("Bray-Curtis PCA analysis (PC2 vs PC1)")

hmplld <- rbind(hmpalltaxjoined,lldTax)
ggplot(data=hmplld,aes(x=Bray_Curtis_PC1,y=Bray_Curtis_PC2,color=cohort)) + geom_point(size=2) + ggtitle("Bray-Curtis PCA analysis (PC2 vs PC1)")
daglld <- rbind(dagTax,lldTax)
ggplot(data=daglld,aes(x=Bray_Curtis_PC1,y=Bray_Curtis_PC2,color=cohort)) + geom_point(size=2) + ggtitle("Bray-Curtis PCA analysis (PC2 vs PC1)")


c = 0
s2 <- hmp1729Tax$sample
for (i in strsplit(s2,'_')) {
#  print (i[1])
  if (i[1] %in% hmp1508Tax$sample) {c =c +1}
}
print (paste("hmp1729Tax has",c," samples in common with hmp1508"))

# check for overlap
dagTaxS = metaPhlanGetSpecies("Pilot_metaphlan.txt",Ttresh = 0.005)
lldTaxS = metaPhlanGetSpecies("LLD_taxonomy_metaphlan2_092017.txt",Ttresh = 0.005)
hmp1508TaxS = metaPhlanGetSpecies("HMP2_IBD_1508_taxprofile.tsv",Ttresh = 0.005)
hmp1729TaxS = metaPhlanGetSpecies("HMP2_IBD_1729_taxprofile.tsv",Ttresh = 0.005)
hmpall = c(hmp1508TaxS,hmp1729TaxS)

while (dev.cur()>1) dev.off()
vd <- venn.diagram(list(DAG3=dagTaxS,HMP1508=hmp1508TaxS,HMP1729=hmp1729TaxS),filename = NULL,fill = 2:4 , alpha = 0.3)
grid.draw(vd)

while (dev.cur()>1) dev.off()
grid.draw(vd <- venn.diagram(list(HMP1508=hmp1508TaxS,HMP1729=hmp1729TaxS),filename = NULL, fill=2:3,alpha=0.3))

while (dev.cur()>1) dev.off()
grid.draw(venn.diagram(list(DAG3=dagTaxS,HMPall=hmpall),filename = NULL,fill = 2:3 , alpha = 0.3))

while (dev.cur()>1) dev.off()
grid.draw(venn.diagram(list(LLD=lldTaxS,HMPall=hmpall,DAG3=dagTaxS),filename = NULL,fill = 2:4 , alpha = 0.3))

while (dev.cur()>1) dev.off()
grid.draw(venn.diagram(list(LLD=lldTaxS,HMP1508=hmp1508TaxS,HMP1729=hmp1729TaxS,DAG3=dagTaxS),filename = NULL,fill = 2:5 , alpha = 0.3))
