Diet-Microbiome Project 
-------------
 
Creator: Laura Bolte

Year: 2018 


 
 
 1.Phenotypic Metadata  
 -------------

**Subset to relevant columns**  
```
library(readxl)  
Food=read_excel("~/Desktop/Data/Metadata/Final_Metadata_IBD_V46.xlsx",sheet = "Completemetadata_drugs.txt")
Food=Food[-c(1,2,3),] 
Food=Food[!is.na(Food$UMCGIBDDNAID),]                                                                         
Food[Food=="NA"]<-NA  
colnames(Food)
Food=Food[,c(10,14,82,84,85,748,11,757,758,754,4,48,33,512:587,589:683,685:696,773:777,779:785)] #keep these
```
**Remove samples with a read depth below 10 M**  
```
Food=Food[!grepl("yes",Food$ReadDepthBelow10M),]                                                                      
Food=Food[,-10] #remove column ReadDepth
```
**Convert variables to numeric or factors**                                                                            
```
str(Food)              
cols1=c(11:12,32:207)     #food groups, age, PFR to numeric
Food[cols]=lapply(Food[cols1],as.numeric)  
cols2=c(10,13:31)         #dietary practices (yes/no/..), gender to factors
Food[cols2]=lapply(Food[cols2],as.factor)                                                                                                                                                                           
```
**Split file into Food groups and Dietary ways** 
```
Foodgroups=Food[,c(1:12,32:207)]                                                        
write.table(Foodgroups,'../Metadata/Subsetted Files/Food/Foodgroups.txt',sep = '\t')
Dietaryways=Food[,c(1:31)]                                                                          
write.table(Dietaryways,'../Metadata/Subsetted Files/Food/Dietaryways.txt',sep = '\t')
```
**Correct food groups for caloric intake**                                                                                  
*Divide all food groups by caloric intake*                                                                        
```
colnames(Foodgroups)                                           
Foodgroups_kcal=Foodgroups[,c(13:163,165:188)]/Foodgroups$SUMOFKCAL                             
MissingCols=Foodgroups[,c(1:12,164)]                    #col 1-12 no food, col 164 SUMOFKCAL 
Foodgroups_corr=cbind(MissingCols, Foodgroups_kcal)     #add non-food columns                                                                  
```
**Subset to 4 groups based on ID's**                                                                                             
*1. IBD-Food:*                                                        
```
IBD=Foodgroups_corr[!grepl("LLDeep",Foodgroups_corr$UMCGIBDResearchIDorLLDeepID),] #remove LLD cohort
IBD=IBD[!grepl("no", IBD$Sequenced),]                   #remove those not sequenced
colnames(IBD)
IBD=IBD[,c(2:5, 10:188)]                                #keep cols ID, sex, age, PFR, diagnosis, cd, uc, food    
```
*1.1 CD-Food:*                                                                          
```
CD=IBD[!grepl("UC",IBD$UCVersusGeneralPopulation),]     #remove UC patients 
CD=CD[!is.na(CD$CDVersusGeneralPopulation),]            #remove NA in CD column  
colnames(CD)
CD=CD[,c(1,5:183)]                                      #keep all cols except diagnosis, cd, uc
CD_Food=as.data.frame(CD)                               #otherwise warning when setting rownames
rownames(CD_Food)=CD_Food$UMCGIBDDNAID
CD_Food$UMCGIBDDNAID=NULL 
write.table(CD_Food,'../Metadata/Subsetted Files/Food/CD_Food.txt',sep = '\t')
```
*1.2 UC-Food:*                                                                       
```
UC=IBD[!grepl("CD",IBD$CDVersusGeneralPopulation),]    #remove CD patients  
UC=UC[!is.na(UC$UCVersusGeneralPopulation),]           #remove NA in UC column
colnames(UC)
UC=UC[,c(1,5:183)]                                     #keep all cols except diagnosis, cd, uc
UC_Food=as.data.frame(UC)
rownames(UC_Food)=UC_Food$UMCGIBDDNAID                 
UC_Food$UMCGIBDDNAID=NULL                                                                               
write.table(UC_Food,'../Metadata/Subsetted Files/Food/UC_Food.txt',sep = '\t')
```
*2. LifeLines-Food*
```
LLD=Foodgroups_corr[!grepl("UMCGIBD",Foodgroups_corrkcal$UMCGIBDResearchIDorLLDeepID),] #remove IBD cohort
LLD=LLD[!grepl("no",LLD$Sequenced),]                   #remove those not sequenced
colnames(LLD)
LLD=LLD[,c(1,6,10:188)]                                #keep cols ID, sex, age, PFR, IBS, food 
LLD_Food=as.data.frame(LLD)                             
rownames(LLD_Food)=LLD_Food$UMCGIBDResearchIDorLLDeepID
LLD_Food$UMCGIBDResearchIDorLLDeepID=NULL 
```
*2.1 IBS-Food*
```
IBS=LLD_Food[grepl("yes",LLD_Food$Irritable_bowel_syndrome),] #keep only IBS yes  
IBS_Food=IBS[,-1]                                             #remove column IBS                      
write.table(IBS_Food,'../Metadata/Subsetted files/Food/IBS_Food.txt',sep = '\t')
```
*2.2 HC-Food*                                                
```
HC=LLD_Food[!grepl("yes",LLD_Food$Irritable_bowel_syndrome),] #delete IBS yes, keep IBS NA and NO  
HC_Food=HC[,-1]                                               #remove column IBS  
write.table(HC_Food,'../Metadata/Subsetted files/Food/HC_Food.txt',sep = '\t')
```
 
 2.Microbial Abundance Metadata   
 -------------

**2.1 Taxonomy - All Levels**                                                                                                 
*N.B.: Taxonomy files I am using have been filtered (step 5 Medication project) and normalized (step 6 Medication project)*  

**Subset and merge**  
*1. IBD-Tax:*
```
tax_IBD=read.table('../Metadata/IBD_filtered_taxonomy_pheno.txt', header=T, sep='\t')
colnames(tax_IBD)
tax_IBD=tax_IBD[,c(1,49:304)] #keep ID, Taxa
rownames(tax_IBD)=tax_IBD$SID
tax_IBD$SID=NULL 
rowSums(tax_IBD)   #check if abundances are relative. N.B.: here normalized data!    
```
*2. LLD-Tax:*                                                                                   
```
tax_LLD=read.table('../Metadata/LLD_filtered_taxonomy_pheno.txt', header=T, sep='\t')
colnames(tax_LLD)
tax_LLD=tax_LLD[,c(1, 49:279)] #keep ID, Taxa
rownames(tax_LLD)=tax_LLD$SID
tax_LLD$SID=NULL 
tax_LLD2=as.data.frame(t(tax_LLD))  
tax_LLD3=tax_LLD2[,-grep("IBS", colnames(tax_LLD2))] #delete Maastricht cohort                    
tax_LLD4=as.data.frame(t(tax_LLD3))  
tax_LLD=tax_LLD4
rowSums(tax_LLD)
```
*3. Merge Taxonomy IBD and LLD:*                                                                                             
N.B.: Metaanalysis will need equal taxa in every cohort. IBD file has more taxa than LLDeep. These 2 files need to be merged by common columns i.e. taxa:   
```
tax_I=as.data.frame(t(tax_IBD))  #set taxa as row.names 
tax_L=as.data.frame(t(tax_LLD))  #set taxa as row.names
tax=merge(tax_I, tax_L, by = "row.names", all=T)
rownames(tax)=tax$Row.names
tax$Row.names=NULL
tax=as.data.frame(t(tax))
which(is.na(tax))                #NA's are created when merging two df's by columns of different lenght and setting all=TRUE
tax[is.na(tax)] <- 0             #change created NA's to zero 
write.table(tax,'../Metadata/Subsetted Files/Taxonomy/Taxonomy_LLD_&_IBD.txt', sep='\t')
```

**2.2 Species**  

**Species including Strains (n=189)**
```
tax2=as.data.frame(t(tax)) #All Taxa (n=287)
Species_strains=tax2[grep('s__', row.names(tax2)),]              
colSums(Species_strains)
Species_strains=as.data.frame(t(Species_strains))
write.table(Species_strains,'../Metadata/Subsetted Files/Species/Species_LLD_&_IBD.txt', sep='\t')
```
**Species without Strains (n=48)** 
```
#Species_no_strains=Species_Strains[!grepl('t__', row.names(Species_Strains)),]  
```

**2.3 Pathways** 

**Subset and merge**

*1. IBD-Path:*  
```
path_IBD=read.table("../Metadata/IBD_filtered_path_pheno.txt", header=T, sep='\t')
colnames(path_IBD)
path_IBD=path_IBD[,c(1,49:382)] #keep these columns
rownames(path_IBD)=path_IBD$SID
path_IBD$SID=NULL
```
*2. LLD-Path:*
```
path_LLD=read.table("../Metadata/LLD_filtered_path_pheno.txt", header=T, sep='\t')
colnames(path_LLD)
path_LLD=path_LLD[,c(1,49:378)] #keep these columns
rownames(path_LLD)=path_LLD$SID
path_LLD$SID=NULL
path_LLD2=as.data.frame(t(path_LLD))  
path_LLD3=path_LLD2[,-grep("_IBS", colnames(path_LLD2))] #delete Maastricht 
path_LLD=as.data.frame(t(path_LLD3))  
```
*3. Merge Pathway files IBD and LLD:*                                                                                           
*N.B.: Metaanalysis will need equal taxa in every cohort. The IBD file has more pathways than the LLDeep file. These 2 files need to be merged by common columns i.e. pathways*                                                                              
```
path_I=as.data.frame(t(path_IBD))  #make pathways the row.names 
path_L=as.data.frame(t(path_LLD))  #make pathways the row.names
path=merge(path_I, path_L, by = "row.names", all=T)
rownames(path)=path$Row.names
path$Row.names=NULL
path=as.data.frame(t(path))
which(is.na(path))                #NA's are created when merging two dataframes by columns of different lenght and setting all=TRUE
path[is.na(path)] <- 0            #change created NA's to zero 
rowSums(path)
colSums(path)
write.table(path,'../Metadata/Subsetted Files/Pathways/Pathways_LLD_&_IBD_not_rel.txt', sep='\t')
```
**Make pathways relative to sum** 
```
path_rel=t(t(path))/colSums(t(path)) #Or instead of colSums(t(Pathways)) -> rowSums(Pathways)
rowSums(path)
rowSums(path_rel)
#check:for Aerobactinsyn.PWY: 
#IBDFEC0115 in path = 5.93. rowSums(path) = 2099.47
#IBDFEC0115 in path_rel = 0.0028 -> 5.93/2099.47 = 0.0028
path_rel=as.data.frame(path_rel) #safe as data frame
write.table(path_rel,'../Metadata/Subsetted files/Pathways/Pathways_LLD_&_IBD_rel_abundance.txt',sep = '\t')
View(path_rel)
```
**2.4 Growth Rates**                                                                                                                    
```
setwd("~/Desktop/Data/Metadata")
growth=read.table("../Metadata/selected_samples_growthrates_2.txt", header=T, sep='\t')
rownames(growth)=growth$X
growth$X=NULL
growth=as.data.frame(t(growth))
growth=growth[,-grep("_IBS", colnames(growth))] #delete Maastricht 
growth=as.data.frame(t(growth))
write.table(growth,'../Metadata/Subsetted Files/Growth Rates/Growth_Rates_LLD_&_IBD_not_rel.txt', sep='\t')
```

 3.Merge Phenotypic and Microbial Abundance Files   
 -------------
 
 **Maaslin input files require:**
 - tsv format
 - saved in a folder which will be the working directory during the Maaslin runs
 - column names: ID's first, than phenotypes, followed by microbial abundances.                                                   
   In this case: col 1-4 (ID, age, gender, read depth), col 5-180 (foods), col 180-467 (taxa) 
 
*1. CD-Food-Tax:*
```
setwd("~/Desktop/Data/Maaslin Files")#CD_Food_Taxonomy
CD_foodtax = merge (CD_Food, tax, by = "row.names", all = FALSE)  
#all = FALSE implies that only rownames that are common in both dataframes will be kept
write.table(CD_foodtax, "CD_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)
```
*2. UC-Food-Tax:*
```
UC_foodtax = merge (UC_Food, tax, by = "row.names", all = FALSE)
write.table(UC_foodtax, "UC_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)
```
*3. IBS_Food_Tax:* 
```
IBS_foodtax = merge (IBS_Food, tax, by = "row.names", all = FALSE)
write.table(IBS_foodtax, "IBS_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)
```
*4. HC-Food-Tax:* 
```
HC_foodtax = merge (HC_Food, tax, by = "row.names", all = FALSE)
write.table(HC_foodtax, "HC_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)
```
**Merge Food with Species, Pathways and Growth rates, accordingly!**


 
 4.MaAsLin  
 -------------


**Config file**
*Text file to be saved as input.read.config in working directory*
```
Matrix: Metadata
Delimiter: TAB
Read_TSV_Columns:Sex-how_often_tea__1

Matrix: Abundance 
Delimiter: TAB 
Read_TSV_Columns:k__Archaea.p__Euryarchaeota-
```

- For species change Abundance to: *k__Archaea.p__Euryarchaeota.c__Methanobacteria.o__Methanobacteriales.f__Methanobacteriaceae.g__Methanobrevibacter.s__Methanobrevibacter_smithii.t__Methanobrevibacter_smithii_unclassified-*
- For growth rates change Abundance to: *butyrate.producing_bacterium_SSC.2-*
- For pathways change Abundance to: *AEROBACTINSYN.PWY-*

**Remove patients with only NAs in Food (not filled FFQ)**

**Univariate Maaslin Runs**

- Packages: activate all the 3 gamlss packages 
```
setwd("/Users/laurabolte/Desktop/Data/Maaslin Files/")  #All tsv files from step 3 are stored here
library(Maaslin) 
Maaslin('CD_Food_Tax.tsv','CD_Food_Tax_output',strInputConfig ='input.read.config',strForcedPredictors = c('Sex','AgeAtFecalSampling','PFReads'), dSignificanceLevel = 1, dMinAbd = 0, dMinSamp = 0, strModelSelection = "none", fAllvAll=TRUE, fZeroInflated=TRUE, strTransform = "none")
```
- Thereafter for UC, IBS, HC 
- For the purpose of a later Metaanalysis I have set dsign to 1. So there is no restriction on the significance and Maaslin gives coefficients and p-values for all univariate food-tax runs. 
- Growth rates are run without fZeroInflated=TRUE*

 
 5.Import all Maaslin output files into R 
 -------------
 
*Import files from folder, put them into a list, and merge them into one dataframe*

- Make sure only the txt.files with individual diet factors agains all taxonomy are in Maaslin output folders
- Transfer all other files such as PDFs and QC-files into a different folder

*get a list of files in a directory* 
```
setwd("~/Desktop/Data/Maaslin Files/Maaslin_Food_Tax/") #Set working dir to that containing all files that need to be merged
file_list <- list.files()  #makes list of files in the directory 
```
**Merge the files into a single dataframe**
*iterate through list of files in the current working directory and put them together to form a dataframe* 
```
filestomerge <- c()
for (fldr in file_list) {
  #print(paste('FOLDER:',fldr))
  #print(paste('./',fldr))
  targetfolder <- paste('./',fldr,sep='')
  fldrfiles <- list.files(targetfolder)
  #print(fldrfiles)
  for (ifile in fldrfiles) {
    t <- paste(targetfolder,'/',ifile,sep='')
    print(t)
    filestomerge <- c(filestomerge,t)
  }
}

mergeddata <- read.csv2(filestomerge[1],sep = '\t',header = TRUE)
fn <- strsplit(filestomerge[1],'/')[[1]][2]
fn <- strsplit(fn,'_')[[1]][1]
print(fn)
mergeddata$filename <- fn
#for (i in c(2,100,200,300,400,500,600,700)) { #length(filestomerge)) {
for (i in 2:length(filestomerge)) {
  tmpdata <- read.csv2(filestomerge[i],sep = '\t',header = TRUE)
  print (filestomerge[i])
  fn <- strsplit(filestomerge[i],'/')[[1]][2]
  fn <- strsplit(fn,'_')[[1]][1]
  tmpdata$filename <- fn
  print(fn)
  mergeddata <- rbind(mergeddata,tmpdata)
}

View(mergeddata)
mergeddata$N <- NULL
mergeddata$N.not.0 <- NULL
mergeddata$Variable <- NULL
mergeddata$Q.value <- NULL

tCD <- mergeddata[mergeddata$filename=="CD",]
tCD$filename <- NULL
View(tCD)
colnames(tCD) <- c("Tax","Diet","CD_Coef","CD_p")

tUC <- mergeddata[mergeddata$filename=="UC",]
tUC$filename <- NULL
colnames(tUC) <- c("Tax","Diet","UC_Coef","UC_p")

tIBS <- mergeddata[mergeddata$filename=="IBS",]
tIBS$filename <- NULL
colnames(tIBS) <- c("Tax","Diet","IBS_Coef","IBS_p")

tHC <- mergeddata[mergeddata$filename=="HC",]
tHC$filename <- NULL
colnames(tHC) <- c("Tax","Diet","HC_Coef","HC_p")

merge1 <- merge(tCD,tUC,by=c("Tax","Diet"))
merge2 <- merge(merge1,tIBS,by=c("Tax","Diet"))
merge3 <- merge(merge2,tHC,by=c("Tax","Diet"))
View(merge3)
#metatable <- data.frame(row.names = c('DietFactor','Taxa','CD_p','CD_coef','UC_p','UC_coef','IBS_p','IBS_coef','HC_p','HC_coef'))
metatable=merge3
View(metatable)
setwd("~/Desktop/Data/Maaslin Files")
write.csv(metatable, '../Maaslin Files/Maaslin_Food_Tax/MergedLargeTable.csv')
```

 6.Metaanalysis   
 -------------

- Random effect meta-analysis 
- Take into account samples sizes (weights) and coefficients (directions) per group
- Adjust for multiple testing 
- Heterogeneity estimation (Cochran's Q)
       
**Taxonomy** 
```
setwd("/Users/laurabolte/Desktop/Data/Maaslin Files/Maaslin_Food_Tax/")
Taxonomy=read.csv('../Maaslin_Food_Tax/MergedLargeTable.csv', header=TRUE, sep=',')
Taxonomy=Taxonomy[-c(1)]
```

**Invert p-values (1-p) if the coefficient of CD, UC, or IBS is different from the coefficient in hc (largest cohort)**  

```
my_results=matrix(nrow=nrow(Taxonomy),ncol=ncol(Taxonomy))

for (x in 1:nrow(Taxonomy)){  
  if (Taxonomy[x,3]>0 & Taxonomy[x,5]>0 & Taxonomy[x,7]>0 & Taxonomy[x,9]>0 | Taxonomy[x,3]<0 & Taxonomy[x,5]<0 & Taxonomy[x,7]<0 & Taxonomy[x,9]<0) {
    my_results[x,3]=Taxonomy[x,3] 
    my_results[x,4]=Taxonomy[x,4]
    my_results[x,5]=Taxonomy[x,5]
    my_results[x,6]=Taxonomy[x,6]
    my_results[x,7]=Taxonomy[x,7]
    my_results[x,8]=Taxonomy[x,8]
    my_results[x,9]=Taxonomy[x,9]
    my_results[x,10]=Taxonomy[x,10]
  } else { 
    my_results[x,3]=Taxonomy[x,3]
    my_results[x,4]=Taxonomy[x,4]
    my_results[x,5]=Taxonomy[x,5]
    my_results[x,6]=Taxonomy[x,6]
    my_results[x,7]=Taxonomy[x,7]
    my_results[x,8]=Taxonomy[x,8]
    my_results[x,9]=Taxonomy[x,9]
    my_results[x,10]=Taxonomy[x,10]
    if (sign(Taxonomy[x,9]) != sign(Taxonomy[x,3])){  #if coef HC (col9) has a different sign (+ or -) than coef CD (col3)
      my_results[x,4]=1-Taxonomy[x,4]                 #than p CD (col4) will be inverted to 1-p 
      my_results[x,3]=Taxonomy[x,3]                   #and coef CD (col3) stays coef CD
    } 
    if (sign(Taxonomy[x,9]) != sign(Taxonomy[x,5])){
      my_results[x,6]=1-Taxonomy[x,6]
      my_results[x,5]=Taxonomy[x,5]
    } 
    if (sign(Taxonomy[x,9]) != sign(Taxonomy[x,7])) {
      my_results[x,8]=1-Taxonomy[x,8]
      my_results[x,7]=Taxonomy[x,7]
    } 
  }
  
} 

colnames(my_results)=colnames(Taxonomy)
my_results=as.data.frame(my_results)
my_results$Tax=Taxonomy$Tax
my_results$Diet=Taxonomy$Diet 
View(my_results)
write.table(my_results,'../Subsetted Files/Taxonomy_inverted_p.txt',sep = '\t') 
```
**Meta-analyse** 
```
Tax=read.csv('../Metaanalysis/Subsetted Files/Taxonomy_inverted_p.txt', header=TRUE, sep='\t')
my_results_meta=as.data.frame(Tax)
my_results_meta$meta_z=0           #c11
my_results_meta$meta_p=0           #c12
W=c(205,124,223,872)               #Weights Taxonomy, Pathways: CD, UC, IBS, HC

for (x in 1:nrow(my_results_meta)) { 
  Wi=sqrt(W)  
  P=c(Tax[x,4],Tax[x,6],Tax[x,8],Tax[x,10])
  Zi=qnorm(1-(P/2))                #Convert p-values to Z-scores                                
  Z=(sum(Zi*Wi)/sqrt(sum(Wi^2)))   #Meta-zscore                    
  MetaP=2*pnorm(-abs(Z))           #Convert Z-score to p-value
  my_results_meta[x,11]=Z
  my_results_meta[x,12]=MetaP 
}

Test=my_results_meta[922,] #same as METAL 
write.table(my_results_meta,'../Metaanalysis/Results/Tax_meta_unadjusted.txt',sep = '\t')
```

 7.Correct Meta-P-values for Multiple Testing   
 -------------

*Obtain FDR-Adjusted p-values* 

**P-adjust per food group**
*1. Split the dataframe*
*2. Apply the p.adjust function to all row-pairs (taxa) in each food-subset*

```
my_adj=as.data.frame(my_results_meta)
my_adj$meta_padj=0

diets <- unique(my_adj$Diet)        #Select all associations of one food group, to adjust for the number of tests

for(diet in diets){
  my_adj$meta_padj[my_adj$Diet == diet] <- p.adjust(my_adj$meta_p[my_adj$Diet == diet], "fdr") 
}

write.table(my_adj,"../Metaanalysis_Heterogeneity/Results/Tax_meta_adj_all.txt",sep = '\t')
```

 8.Cochran's Q-Test   
 -------------

*Measuring the inconsistency (heterogeneity) of studies’ results. Heterogeneity in meta-analysis refers to the variation in study outcomes between studies. Cochran’s Q is calculated as the weighted sum of squared differences between individual study effects and the pooled effect across studies with the weights being those used in the pooling method. Q is distributed as a chi-square statistic with k (number of studies) minus 1 degrees of freedom*

**Filter for significant results**
*Filter results that are significant after FDR correction. On these results the Cochran's Q-test will be performed*
```
Tax_sign=my_adj[my_adj$meta_padj<=0.05,]        
write.table(Tax_sign,"../Metaanalysis/Results/Tax_meta_adj.txt",sep = '\t')
```

**Cochran's Q-Test**
```
my_results_Q=as.data.frame(Tax_sign)   
my_results_Q$weighted_z=0    #c14
my_results_Q$het_chisq=0     #c15
my_results_Q$het_p=0         #c16

W=c(205, 124, 223, 872)      #weights CD, UC, IBS, HC
totalSample=sum(W)

for (x in 1:nrow(my_results_Q)) { 
  P=c(Tax_sign[x,4],Tax_sign[x,6],Tax_sign[x,8],Tax_sign[x,10])
  Zi=qnorm(1-(P/2))                                 #Convert p-values to Z-scores                                
  weightedZ= sum(sqrt(W)*Zi)                         
  expectedZ= sqrt(W)*weightedZ/totalSample          #Calculate expected Z 
  hetchisq= sum((Zi-expectedZ)*(Zi-expectedZ))      #Q-value 
  hetpval=pchisq(hetchisq, lower.tail = F, df=3)    #P-value of Q
  my_results_Q[x,14]=weightedZ  
  my_results_Q[x,15]=hetchisq 
  my_results_Q[x,16]=hetpval
}

#Test2=my_results_Q[, ] #check if same as METAL in row O_Actinomycetales - fruit
write.table(my_results_Q,"../Metaanalysis/Results/Tax_Meta_Chochrans_unadjusted.txt",sep = '\t')

```
The tool METAL (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922887/) can be used to check results 

 9.Correct Heterogeneity-P-values for Multiple Testing    
 -------------
**Null hypothesis: Is there homogeneity?** 
*If not (=significant het_p), reject these because they are heterogenous*
```
Tax_Q_adj=as.data.frame(my_results_Q)
Tax_Q_adj$het_padj=0

diets <- unique(Tax_Q_adj$Diet)

for(diet in diets){
  Tax_Q_adj$het_padj[Tax_Q_adj$Diet == diet] <- p.adjust(Tax_Q_adj$het_p[Tax_Q_adj$Diet == diet], "bonferroni") 
}
```
**Take out results with significant heterogeneity, keep only homogenous results**
```
Tax_Q_adj_insign=Tax_Q_adj[Tax_Q_adj$het_padj > 0.1,]  #above 0.1 = no significant heterogeneity  
write.table(Tax_Q_adj_insign,"../Metaanalysis/Results/Tax_Meta_Homogeneous.txt",sep = '\t')
```
 
 10.Heatmaps    
 -------------
 
**Modify output files of Meta-analysis**                                                                                      
*CD/UC/IBS/HC:*
- color: their coefficients
- colordepth: their original, non-inverted p-vals (go back to original Maaslin table)                                          
*Meta-p:*
- color: coefficient of HC
- colordepth: meta-p 

*STEP 1: Subset and merge metaanalysis table for heatmaps:*
```
library(reshape2)
library(dplyr)
library(ggplot2)
library(gplots)
library(gridExtra)
library(tidyr)

#1.1: Plot only these results, take meta-p and coefficients from here: 
setwd("~/Desktop/Data/Heatmaps")
Species_results=read.table("../Heatmaps/Spec_Meta_Cochrans_adjusted_homogeneous.txt", sep='\t', header=T)
Species_m_c = select(Species_results,-4,-6,-8,-10,-12:-15) #meta-p values, only significant and homogeneous results 
Species_m_c$NewCol = do.call(paste, c(Species_m_c[c("Tax", "Diet")], sep = "__|")) #join strings from Tax and Diet column into one new column
row.names(Species_m_c)=Species_m_c$NewCol
Species_m_c=Species_m_c[,-8]

#1.2: Take individual, non-inverted p-values from here: 
Species_p=read.csv("~/Desktop/Data/Heatmaps/MergedLargeTable_Species.csv")
Species_p=Species_p[,-1]                 
Species_p=select(Species_p,-3,-5,-7,-9) #non-inverted original p-values in cohorts
Species_p$NewCol = do.call(paste, c(Species_p[c("Tax", "Diet")], sep = "__|")) #join strings from Tax and Diet column into one new column
row.names(Species_p)=Species_p$NewCol
Species_p=Species_p[,-7]

#1.3: Merge by NewCol and remove that column afterwards: 
Species=merge(Species_m_c, Species_p, by = "row.names")  #all=T
Species=select(Species,-1,-9,-10)

#1.4: Shorten column names in excel and import back in R  
write.table(Species,"../Heatmaps/Species_map.txt", sep='\t')
```

**Loop to Plot heatmaps** 
```
Speciesmap=read.table("../Heatmaps/Species_strains_map.txt", sep='\t', header=T)
Speciesmap=Speciesmap[,-12]

Speciesmap$CD_color="grey"
Speciesmap$UC_color="grey"
Speciesmap$IBS_color="grey"
Speciesmap$HC_color="grey"
Speciesmap$Metap_color="grey"

Speciesmap_temp=Speciesmap
diets <- unique(Speciesmap_temp$Diet)

for(diet in diets[1:3]){   #takes only 1 to 3 food groups, for all remove [1:3] 
  Speciesmap=Speciesmap_temp[Speciesmap_temp$Diet == diet,] 
  print(diet)
  print(nrow(Speciesmap))

#CD
if (sum(Speciesmap$CD_p > 0.05 & Speciesmap$CD_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$CD_p > 0.05 & Speciesmap$CD_Coef > 0),]$CD_color<-"1"
}
if (sum(Speciesmap$CD_p <= 0.05 & Speciesmap$CD_p > 0.00005 & Speciesmap$CD_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$CD_p <= 0.05 & Speciesmap$CD_p > 0.00005 & Speciesmap$CD_Coef > 0),]$CD_color<-"2"
}
if (sum(Speciesmap$CD_p <= 0.00005 & Speciesmap$CD_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$CD_p <= 0.00005 & Speciesmap$CD_Coef > 0),]$CD_color<-"3"
}
if (sum(Speciesmap$CD_p > 0.05 & Speciesmap$CD_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$CD_p > 0.05 & Speciesmap$CD_Coef < 0),]$CD_color<-"-1"
}
if (sum(Speciesmap$CD_p <= 0.05 & Speciesmap$CD_p > 0.00005 & Speciesmap$CD_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$CD_p <= 0.05 & Speciesmap$CD_p > 0.00005 & Speciesmap$CD_Coef < 0),]$CD_color <-"-2"
  }
if (sum(Speciesmap$CD_p <= 0.00005 & Speciesmap$CD_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$CD_p <= 0.00005 & Speciesmap$CD_Coef < 0),]$CD_color <-"-3"
}

#UC
if (sum(Speciesmap$UC_p > 0.05 & Speciesmap$UC_Coef > 0) > 0) {
    Speciesmap[(Speciesmap$UC_p > 0.05 & Speciesmap$UC_Coef > 0),]$UC_color<-"1"
}
if (sum(Speciesmap$UC_p <= 0.05 & Speciesmap$UC_p > 0.00005 & Speciesmap$UC_Coef > 0) > 0) {
    Speciesmap[(Speciesmap$UC_p <= 0.05 & Speciesmap$UC_p > 0.00005 & Speciesmap$UC_Coef > 0),]$UC_color<-"2"
}
if (sum(Speciesmap$UC_p <= 0.00005 & Speciesmap$UC_Coef > 0) > 0) {
Speciesmap[(Speciesmap$UC_p <= 0.00005 & Speciesmap$UC_Coef > 0),]$UC_color<-"3"
}
if (sum(Speciesmap$UC_p > 0.05 & Speciesmap$UC_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$UC_p > 0.05 & Speciesmap$UC_Coef < 0),]$UC_color<-"-1"
}
if (sum(Speciesmap$UC_p <= 0.05 & Speciesmap$UC_p > 0.00005 & Speciesmap$UC_Coef < 0) > 0) {
    Speciesmap[(Speciesmap$UC_p <= 0.05 & Speciesmap$UC_p > 0.00005 & Speciesmap$UC_Coef < 0),]$UC_color <-"-2"
}
if (sum(Speciesmap$UC_p <= 0.00005 & Speciesmap$UC_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$UC_p <= 0.00005 & Speciesmap$UC_Coef < 0),]$UC_color <-"-3"
}

#IBS
if (sum(Speciesmap$IBS_p > 0.05 & Speciesmap$IBS_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$IBS_p > 0.05 & Speciesmap$IBS_Coef > 0),]$IBS_color<-"1"
}
if (sum(Speciesmap$IBS_p <= 0.05 & Speciesmap$IBS_p > 0.00005 & Speciesmap$IBS_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$IBS_p <= 0.05 & Speciesmap$IBS_p > 0.00005 & Speciesmap$IBS_Coef > 0),]$IBS_color<-"2"
}
if (sum(Speciesmap$IBS_p <= 0.00005 & Speciesmap$IBS_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$IBS_p <= 0.00005 & Speciesmap$IBS_Coef > 0),]$IBS_color<-"3"
}
if (sum(Speciesmap$IBS_p > 0.05 & Speciesmap$IBS_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$IBS_p > 0.05 & Speciesmap$IBS_Coef < 0),]$IBS_color<-"-1"
}
if (sum(Speciesmap$IBS_p <= 0.05 & Speciesmap$IBS_p > 0.00005 & Speciesmap$IBS_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$IBS_p <= 0.05 & Speciesmap$IBS_p > 0.00005 & Speciesmap$IBS_Coef < 0),]$IBS_color <-"-2"
}
if (sum(Speciesmap$IBS_p <= 0.00005 & Speciesmap$IBS_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$IBS_p <= 0.00005 & Speciesmap$IBS_Coef < 0),]$IBS_color <-"-3"
}

#HC
if (sum(Speciesmap$HC_p > 0.05 & Speciesmap$HC_Coef > 0)> 0) {
  Speciesmap[(Speciesmap$HC_p > 0.05 & Speciesmap$HC_Coef > 0),]$HC_color<-"1"
}
if (sum(Speciesmap$HC_p <= 0.05 & Speciesmap$HC_p > 0.00005 & Speciesmap$HC_Coef > 0)> 0) {
  Speciesmap[(Speciesmap$HC_p <= 0.05 & Speciesmap$HC_p > 0.00005 & Speciesmap$HC_Coef > 0),]$HC_color<-"2"
}
if (sum(Speciesmap$HC_p <= 0.00005 & Speciesmap$HC_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$HC_p <= 0.00005 & Speciesmap$HC_Coef > 0),]$HC_color<-"3"
}
if (sum(Speciesmap$HC_p > 0.05 & Speciesmap$HC_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$HC_p > 0.05 & Speciesmap$HC_Coef < 0),]$HC_color<-"-1"
}
if (sum(Speciesmap$HC_p <= 0.05 & Speciesmap$HC_p > 0.00005 & Speciesmap$HC_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$HC_p <= 0.05 & Speciesmap$HC_p > 0.00005 & Speciesmap$HC_Coef < 0),]$HC_color <-"-2"
}
if (sum(Speciesmap$HC_p <= 0.00005 & Speciesmap$HC_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$HC_p <= 0.00005 & Speciesmap$HC_Coef < 0),]$HC_color <-"-3"
}

#Meta
if (sum(Speciesmap$meta_p > 0.05 & Speciesmap$HC_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$meta_p > 0.05 & Speciesmap$HC_Coef > 0),]$Metap_color<-"1"
}
if (sum( (Speciesmap$meta_p <= 0.05 & Speciesmap$meta_p > 0.00005) & Speciesmap$HC_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$meta_p <= 0.05 & Speciesmap$meta_p > 0.00005 & Speciesmap$HC_Coef > 0),]$Metap_color<-"2"
}
if (sum(Speciesmap$meta_p <= 0.00005 & Speciesmap$HC_Coef > 0) > 0) {
  Speciesmap[(Speciesmap$meta_p <= 0.00005 & Speciesmap$HC_Coef > 0),]$Metap_color<-"3"
}
if (sum((Speciesmap$meta_p > 0.05 & Speciesmap$HC_Coef < 0)) > 0) {
  Speciesmap[(Speciesmap$meta_p > 0.05 & Speciesmap$HC_Coef < 0),]$Metap_color<-"-1"
}
if (sum(Speciesmap$meta_p <= 0.05 & Speciesmap$meta_p > 0.00005 & Speciesmap$HC_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$meta_p <= 0.05 & Speciesmap$meta_p > 0.00005 & Speciesmap$HC_Coef < 0),]$Metap_color <-"-2"
}
if (sum(Speciesmap$meta_p <= 0.00005 & Speciesmap$HC_Coef < 0) > 0) {
  Speciesmap[(Speciesmap$meta_p <= 0.00005 & Speciesmap$HC_Coef < 0),]$Metap_color <-"-3"
}

for_plot=Speciesmap[,c(1,2,12,13,14,15,16)]
a_for_plot=gather(for_plot,"Cohort","Color",CD_color:Metap_color)
g <- ggplot (a_for_plot, aes(x=Cohort, y=Species)) + geom_tile(aes(fill=Color), colour="white") + 
scale_fill_manual(breaks=c("-1","-2","-3","1","2","3"), values=c("#deebf7","#9ecae1","#3182bd", "#fee0d2", "#fc9272", "#de2d26"), name="p-value", labels=c("P > 5e-02","P < 5e-02", "P < 5e-05","P > 5e-02","P < 5e-02", "P < 5e-05")) + theme(panel.background=element_blank(), axis.text=element_text(colour="black")) + labs(x="Dataset", y="Species") + scale_x_discrete (labels=c("CD (205)","UC (154)","IBS (223)", "HC (871)","Meta (1453)")) + ggtitle(paste("Diet:",diet))
ggsave(g,filename = paste('plot_species_diet_',diet,'.png',sep=''),width = 8,height=nrow(a_for_plot)*0.045)

}

Speciesmap=Species_temp
```

 11.Hierarchial Clustering
 -------------
 
**11.1 Subset**
Input files 
rows: matched IDs, columns: foods

**11.2 CLUSTERING DIET** 

1. Load Input Files 
```
setwd("~/Desktop/Data")
All_Food=read.table('../Data/Clustering/Input files/Food_all_cohorts_Final.txt', header=T, sep='\t')
CD_Food=read.table('../Data/Clustering/Input files/CD_Food.txt', header=T, sep='\t')
UC_Food=read.table('../Data/Clustering/Input files/UC_Food.txt', header=T, sep='\t')
IBS_Food=read.table('../Data/Clustering/Input files/IBS_Food.txt', header=T, sep='\t')
HC_Food=read.table('../Data/Clustering/Input files/HC_Food.txt', header=T, sep='\t')
```

2. Load Cluster Function

Optional: Save the function as seperate script and load it from file using the source() function
```
source("Clustering/ClusterFun_Diet_version1.R") 
```

```
do.clustering <- function(diet_data,
                          dist.method="euclidean",        #Ordinary straight line distance between 2 points in the euclidean space. With this distance, euclidean space becomes a metric space. root of (x1-x2)^2 + (y1-y2)^2
                          cluster.method = "complete",    
                          h= 0.5){                        #Clustering height at which the tree is cut into groups. The higher this number, the more stringent and the less groups.
  if (dist.method == "euclidean"){         
    data.dist = dist(t(scale(diet_data)))                 #Scale variables to make them comparable. E.g. SumKcal has values of 1000-2000, group bread e.g. 150 grams per day, how_often_potatoes 5 per day etc.)
                                                          #t: dist() calculating the distance matrix expects vectors to be horizontal. Transform data so that cols=id, rows=foods 
                                                          #Compute distance between each sample 
  } else {
    data.dist = as.dist(1-cor(diet_data)^2)               #else: correlation based distance, considers 2 ojects similar if their features are highly correlated, even though the observed values mau be far apart in euclidean distance. The lower the distance e.g. 0 the more correlated. 
  }
  hclust.object=hclust(data.dist,method = cluster.method) #Cluster similar foods into groups
  plot(as.dendrogram(hclust.object), cex=0.5, hang=-1, width = 9)  #horiz=T, type="triangle", hang=-1, error: "hang" is not a graphical parameter
  cutree_returned = cutree(hclust.object,h=h)                      #Cut the tree resulting from hclust into several group by specifying h 
  cutree_returned
}
```

*3. Apply Cluster Function*
```
All_Food=as.matrix(All_Food)
dev.off()
library(genefilter)
rv <- rowVars(All_Food)
idx <- order(-rv)
#pick distance method, cluster method, and cuttree height. All_Foodault: euclidean, complete 
do.clustering(All_Food[idx,])->s1       #idx=index 
do.clustering(All_Food[idx,],h=53)->s1  #h = the higher the less clusters: 0.2 = many clusters, 60 = 27 clusters 
s1                                      #enter to show clusters

library(RColorBrewer)
d=dist(t(scale(All_Food)))
hmcol=brewer.pal(11,"RdYl")
x=as.matrix(d)
heatmap(d[cutree(d,k=2)==1,], col=hmcol)
```

*4. Cluster Matrix - Obtain values in a table*
```
s2 = data.frame(names = names(s1),cluster=s1)  #s1 = dietary clusters -> as table s2
s2[,1] = as.character(s2[,1])
```

*5. Load Function to Calculate Cluster Representatives or Centroids and Obtain Matrices per cluster stored in "t1"* 
```
select_representatives = function(All_Food,clusters,type = "centroid"){  
  final.output = list()
  library(foreach)                                                     
  clusters = data.frame(names = names(clusters),cluster = clusters)      #clusters equals s2
  clusters[,1] = as.character(clusters[,1])
  
  output = foreach(i = 1:max(clusters[,2]),.combine = cbind) %do% {
    
    cluster.matrix = All_Food[,s2[s2$cluster==i,1]]                      #cluster.matrix = rows: all sample-IDs, cols: all foods belonging to one cluster
    centroid = rowMeans(scale(cluster.matrix))                           #opt.A) centroid of the cluster, per row = the mean (average) of all foods within a cluster, for each individual (per row)
    if (type =="centroid"){
      out = data.frame(value = centroid)
      colnames(out) = paste0("centroid.",i)
      out
    } else if (type == "repr") {                                         #opt.B) representative by euclidean distance between cluster members (foods) 
      cluster.matrix.centroid = cbind(centroid,scale(cluster.matrix))
      dist2centroid=as.matrix(dist(t(cluster.matrix.centroid)))[2:ncol(cluster.matrix.centroid),1,drop=F]
      repr.name = rownames(dist2centroid)[which.min(dist2centroid)]      #takes name (repr.name) of the food with the minimal distance between samples
      out = All_Food[,repr.name,drop=F]                                  #... and couples it to All_Food 
      out                                                                #in case of cluster1: fruit (35.46542) = minimal distance in dist2centroid 
    } else {stop("incorrect type!")}
    
  }
  final.output[["repr.table"]] = output
  final.output[["clusters"]] = clusters
  final.output[["items.matrix"]] = All_Food[,rownames(clusters)]
  colnames(final.output[["items.matrix"]]) = paste0(clusters[,2],".",colnames(final.output[["items.matrix"]]))
  final.output[["items.matrix"]] = final.output[["items.matrix"]][,sort(colnames(final.output[["items.matrix"]]))]
  final.output
}
```

*6. Apply function*
```
t1 = select_representatives(All_Food,s1,type = "centroid") #Or: type = "repr"
Centroids_All_Food=as.data.frame(t1$repr.table)            #Or: Repr_All_Food=as.data.frame(t1$repr.table)
write.table(Centroids_All_Food,'../Data/Clustering/Centroids_All_Foods_H53.txt', sep='\t') 
Clusters_All_Food=as.data.frame(t1$items.matrix)
write.table(Clusters_All_Food,'../Data/Clustering/Clusters_All_Foods_H53.txt', sep='\t')
```

**11.3 CLUSTERING SPECIES** 

*1. Load Input File* 
setwd("~/Desktop/Data")
```
CD_Spec=read.table('../Data/Clustering/Input files/CD_Specstr_final.txt', header=T, sep='\t', stringsAsFactors = F)
UC_Spec=read.table('../Data/Clustering/Input files/UC_Specstr_final.txt', header=T, sep='\t', stringsAsFactors = F)
IBS_Spec=read.table('../Data/Clustering/Input files/IBS_Specstr_final.txt', header=T, sep='\t', stringsAsFactors = F)
HC_Spec=read.table('../Data/Clustering/Input files/HC_Specstr_final.txt', header=T, sep='\t', stringsAsFactors = F)
All_Spec=read.table('../Data/Clustering/Input files/All_Spec_final.txt', header=T, sep='\t')
```

*2. Load Cluster Function Microbiome Clusters*
```
library(vegan)

do.clustering <- function(data,
                          dist.method="euclidean",   
                          cluster.method = "complete", 
                          h= 0.5){              
  if (dist.method == "euclidean"){          
    data.dist = dist(t(data))                    

  } else if (dist.method=="bray-curtis") {
    data.dist = vegdist(t(data))
  } else if (dist == "correlation") {
    data.dist = as.dist(1-cor(data)^2)         
  } else {stop("not the correct distance")}
  hclust.object=hclust(data.dist,method = cluster.method) 
  
  plot(as.dendrogram(hclust.object), cex=0.5, horiz=T)  
  cutree_returned = cutree(hclust.object,h=h)  
  cutree_returned
}
```

*3. Apply Cluster Function*
```
dev.off()
library(genefilter)
rv <- rowVars(All_Spec)
idx <- order(-rv)
#pick distance method, cluster method, and cuttree height. 
#do.clustering(All_Spec[idx,], dist.method = "bray-curtis")->s1         #OR: dist.method = "euclidean", h=4 
do.clustering(All_Spec[idx,], dist.method = "bray-curtis", h=0.8)->s1   #yields 29 clusters, the higher "h", the less clusters 
s1
```

*4. Cluster Matrix - Obtain values in a table*
```
s2 = data.frame(names = names(s1),cluster=s1)  #s1 = Species clusters -> as table s2
s2[,1] = as.character(s2[,1])
```

*5. Load Function to Calculate Cluster Representatives or Centroids*
```
select_representatives = function(All_Spec,clusters,type = "repr"){  
  final.output = list()
  library(foreach)                                                     
  clusters = data.frame(names = names(clusters),cluster = clusters)      #clusters equals s2
  clusters[,1] = as.character(clusters[,1])
  
  output = foreach(i = 1:max(clusters[,2]),.combine = cbind) %do% {
    
    cluster.matrix = All_Spec[,s2[s2$cluster==i,1]]                      #cluster.matrix = rows: all sample-IDs, cols: all species of one cluster
    centroid = rowMeans(scale(cluster.matrix))                           #option A) centroid of the cluster, per row = the mean (average) of all foods within a cluster, for each individual (per row)
    if (type =="centroid"){
      out = data.frame(value = centroid)
      colnames(out) = paste0("centroid.",i)
      out
    } else if (type=="repr") {                                           #option B) representative   
      if (sum(s2$cluster==i)==1) {
        repr.name=rownames(s2[s2$cluster==i,])                           #if the cluster contains only 1 species, take the name of that species as representative
        out = All_Spec[,repr.name,drop=F]                                #... and couple it to All_Food for connection with ID's   
        out                                                              
      } else {
        cluster.matrix.centroid = cbind(centroid,scale(cluster.matrix))  #if the cluster contains more than 1 species
        dist2centroid=as.matrix(dist(t(cluster.matrix.centroid)))[2:ncol(cluster.matrix.centroid),1,drop=F] 
                                                                         #calculate representative, which has minimal euclidean distance to other cluster members  
        repr.name = rownames(dist2centroid)[which.min(dist2centroid)]      
        out = All_Spec[,repr.name,drop=F]                                #... and couple it to All_Food for connection with ID's 
        out
      }                                                                 
    } else {stop("incorrect type!")}
  }
  final.output[["repr.table"]] = output
  final.output[["clusters"]] = clusters
  final.output[["items.matrix"]] = All_Spec[,rownames(clusters)]
  colnames(final.output[["items.matrix"]]) = paste0(clusters[,2],".",colnames(final.output[["items.matrix"]]))
  final.output[["items.matrix"]] = final.output[["items.matrix"]][,sort(colnames(final.output[["items.matrix"]]))]
  final.output
}

t1 = select_representatives(All_Spec,s1,type= "centroid")                   #"OR: type = "repr"
Centroids_All_Spec=as.data.frame(t1$repr)
write.table(Centroids_All_Spec,'../Data/Clustering/Centroids_All_Specs_H_0.8.txt', sep='\t')
Clusters_All_Spec=as.data.frame(t1$items.matrix)
#Clusters_All_Spec=t(Clusters_All_Spec)
write.table(Clusters_All_Spec,'../Data/Clustering/Clusters_All_Specs_H_0.8.txt', sep='\t')
```

**11.4. CORRELATING FOOD CLUSTERS WITH SPECIES CLUSTERS**
```
#1. Load Input Files 
setwd("~/Desktop/Data")
food=read.table('../Data/Clustering/Centr_Corr/Centroids_All_Foods.txt', sep='\t', header=T)
row.names(food)=food$Row.Names
food$Row.Names=NULL

mb=read.table('../Data/Clustering/Centr_Corr/Centroids_All_Specs_h_0.8.txt', sep='\t', header=T)
row.names(mb)=mb$Row.Names
mb$Row.Names=NULL

covariates=read.table("../Data/Clustering/Covariates.txt", sep='\t', header=T) 
```

*2. Check distribution*
```
#plot(covariates$AgeAtFecalSampling) #normal
#hist(covariates$AgeAtFecalSampling)  
#plot(covariates$SUMOFKCAL)          #outliers -> check who and normalize
#hist(covariates$AgeAtFecalSampling)
```

*3. Load libraries*
Consider that in the cluster analysis part, cutting the tree lower (lower h), will create more clusters and stronger correlation

```
library(gdata)
library (corrplot)
library (plyr)
library (psych)
library(RColorBrewer) 
library(foreach) 
library(ppcor)

file=cbind(mb, food, covariates)

#Test: one species eg C1 with one diet cluster eg C1 
test=pcor.test(x=file$Dorea_formicigenerans,y=file$Sweetened_tea_yoghurtdrink,z=file[, c("AgeAtFecalSampling","SUMOFKCAL","Sex2")], method=c("spearman"))
```

*4. Loop through all clusters and perform partial correlations*
```
my_results=matrix(,nrow = 29, ncol=25)   #29 Species clusters and 25 Food clusters
my_results_p=matrix(,nrow = 29, ncol=25)

for (i in 1:29){    #columns containing bacterial clusters
  x=1
  for (a in 30:54){ #columns containing food clusters
    my_temp=pcor.test(x=file[,i],y=file[,a],z=file[, c("AgeAtFecalSampling","SUMOFKCAL","Sex2")], method=c("spearman"))
    my_results[i,x]=my_temp[,1]
    my_results_p[i,x]=my_temp[,2]
    x=x+1
  }
}
colnames(my_results)=colnames(file)[30:54]
rownames(my_results)=colnames(file)[1:29]
colnames(my_results_p)=colnames(file)[30:54]
rownames(my_results_p)=colnames(file)[1:29]
write.table(my_results_p, "~/Desktop/Data/Clustering/Centr_Corr/Corr_P.txt", sep='\t')
write.table(my_results, "~/Desktop/Data/Clustering/Centr_Corr/Corr_R.txt", sep='\t')
```

*6. Adjust for multiple testing*
```
my_results_fdr=matrix(p.adjust(as.vector(as.matrix(my_results_p)), method='fdr',n=725),ncol=25) #29 species clusters times 25 food clusters = 725 tests
colnames(my_results_fdr)=colnames(my_results_p)
rownames(my_results_fdr)=rownames(my_results_p)
View(my_results_fdr)
write.table(my_results_fdr, "~/Desktop/Data/Clustering/Centr_Corr/Corr_FDR.txt", sep='\t')
```

*7. Show number of significant correlations (fdr <0.05)*
```
sum(my_results_fdr <= 0.05)  #53 significant results  #compared to representative method which yields 23 results

#Correlation plot 
dev.off()
#Color
display.brewer.all() #display all colour schemes
col1 = rev(brewer.pal(n=6, name="RdYlBu")) #reverse colors so that positive association = red, negative = blue)
color = c(col1,col1,col1) #repeat col to extend the range
#cl.pos="b"               #plots legenda below plot instead of right side
```

*8. Corrplot* 
```
#Leave insignificant results blanc by insig='blank'
CorrClust_sign=corrplot(my_results, p.mat=my_results_fdr, sig.level=0.05, insig='blank',cl.lim=c(-0.3, 0.3), method="color",tl.cex = 0.6, cl.cex=0.5,cl.align.text = "l",col = color,mar=c(0, 0, 0, 0),na.label= "square", na.label.col = "white")
#mar = c(2,0,1,0)
#title='Diet-Species Clusters'
CorrClust_all=corrplot(my_results, cl.lim=c(-0.3, 0.3), method="color", tl.cex = 0.6, cl.cex=0.5, cl.align.text = "l", na.label= "square", na.label.col = "white", col = color,mar = c(2,0,1,0))

#non-square matrix
CorrClusters=corrplot(my_results, p.mat=my_results_fdr, sig.level=0.05, insig = 'blank', cl.lim=c(-0.3, 0.3), method="color", tl.cex = 0.4, cl.cex=0.7, cl.align.text = "l", na.label= "square", na.label.col = "white", order = "hclust", col = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(200),mar = c(2,0,1,0))
##make colors darker  
CorrClusters=corrplot(cor1, p.mat=pvalue1, sig.level=0.05, insig = 'blank', cl.lim=colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(200),mar = c(2,0,1,0))
col3=brewer.pal(n=8, name="RdYlBu")
color = c(col3,col3,col3) ##repeat col3(-0.3, 0.3), method="color", tl.cex = 0.7, cl.cex=0.7, cl.align.text = "l", na.label= "square", na.label.col = "white", order = "hclust", col = c to extend the range
```
