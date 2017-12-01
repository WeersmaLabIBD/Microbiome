Diet-Microbiome Project 
-------------
 
Creator: Laura Bolte

Year: 2017 


 
 
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
Food=Food[!grepl("yes",Food$ReadDepthBelow10M),]                                                                      
Food=Food[,-10] #remove column ReadDepth

**Convert food groups, age, and PFR to numeric**                                                                            
str(Food)              
cols=c(11:12,32:207)  
Food[cols]=lapply(Food[cols],as.numeric)  
*Save as Foodgroups*                                                                                    
Foodgroups=Food[,c(1:12,32:207)]                                                                   
write.table(Foodgroups,'../Metadata/Subsetted Files/Food/Foodgroups.txt',sep = '\t')
**Convert dietary practices, and gender to factors**                                                               
cols2=c(10,13:31)     
Food[cols2]=lapply(Food[cols2],as.factor)                                                                                      
*Save as Dietaryways*                                                                                     
Dietaryways=Food[,c(1:31)]                                                                           
write.table(Dietaryways,'../Metadata/Subsetted Files/Food/Dietaryways.txt',sep = '\t')

**Correct food groups for caloric intake**                                                                                  
*Divide all food groups by caloric intake*                                                                        
colnames(Foodgroups)                                           
Foodgroups_kcal=Foodgroups[,c(13:163,165:188)]/Foodgroups$SUMOFKCAL                             
MissingCols=Foodgroups[,c(1:12,164)]                    #col 1-12 no food, col 164 SUMOFKCAL 
Foodgroups_corr=cbind(MissingCols, Foodgroups_kcal)     #add non-food columns                                                                  

**Subset to 4 groups based on ID's**                                                                                             
*1. IBD-Food:*                                                        
IBD=Foodgroups_corr[!grepl("LLDeep",Foodgroups_corr$UMCGIBDResearchIDorLLDeepID),]                                          
                                                        #remove LLD cohort
IBD=IBD[!grepl("no", IBD$Sequenced),]                   #remove those not sequenced
colnames(IBD)
IBD=IBD[,c(2:5, 10:188)]                                #keep cols ID, sex, age, PFR, diagnosis, cd, uc, food    

*1.1 CD-Food:*                                                                          
CD=IBD[!grepl("UC",IBD$UCVersusGeneralPopulation),]     #remove UC patients 
CD=CD[!is.na(CD$CDVersusGeneralPopulation),]            #remove NA in CD column  
colnames(CD)
CD=CD[,c(1,5:183)]                                      #keep all cols except diagnosis, cd, uc
CD_Food=as.data.frame(CD)                               #otherwise warning when setting rownames
rownames(CD_Food)=CD_Food$UMCGIBDDNAID
CD_Food$UMCGIBDDNAID=NULL 
write.table(CD_Food,'../Metadata/Subsetted Files/Food/CD_Food.txt',sep = '\t')

*1.2 UC-Food:*                                                                       
UC=IBD[!grepl("CD",IBD$CDVersusGeneralPopulation),]    #remove CD patients  
UC=UC[!is.na(UC$UCVersusGeneralPopulation),]           #remove NA in UC column
colnames(UC)
UC=UC[,c(1,5:183)]                                     #keep all cols except diagnosis, cd, uc
UC_Food=as.data.frame(UC)
rownames(UC_Food)=UC_Food$UMCGIBDDNAID                 
UC_Food$UMCGIBDDNAID=NULL                                                                               
write.table(UC_Food,'../Metadata/Subsetted Files/Food/UC_Food.txt',sep = '\t')

*2. LifeLines-Food*
LLD=Foodgroups_corr[!grepl("UMCGIBD",Foodgroups_corrkcal$UMCGIBDResearchIDorLLDeepID),]                                 
                                                       #remove IBD cohort
LLD=LLD[!grepl("no",LLD$Sequenced),]                   #remove those not sequenced
colnames(LLD)
LLD=LLD[,c(1,6,10:188)]                                #keep cols ID, sex, age, PFR, IBS, food 
LLD_Food=as.data.frame(LLD)                             
rownames(LLD_Food)=LLD_Food$UMCGIBDResearchIDorLLDeepID
LLD_Food$UMCGIBDResearchIDorLLDeepID=NULL 

*2.1 IBS-Food*
IBS=LLD_Food[grepl("yes",LLD_Food$Irritable_bowel_syndrome),] #keep only IBS yes  
IBS_Food=IBS[,-1]                                             #remove column IBS                      
write.table(IBS_Food,'../Metadata/Subsetted files/Food/IBS_Food.txt',sep = '\t')

*2.2 HC-Food*                                                
HC=LLD_Food[!grepl("yes",LLD_Food$Irritable_bowel_syndrome),] #delete IBS yes, keep IBS NA and NO  
HC_Food=HC[,-1]                                               #remove column IBS  
write.table(HC_Food,'../Metadata/Subsetted files/Food/HC_Food.txt',sep = '\t')

 
 2.Microbial Abundance Metadata   
 -------------

**2.1 Taxonomy - All Levels**  
*N.B.: Taxonomy files I am using have been filtered (step 5 Medication project) and normalized (step 6 Medication project)  

**Subset to relevant columns**  
*1. IBD-Tax:*
tax_IBD=read.table('../Metadata/IBD_filtered_taxonomy_pheno.txt', header=T, sep='\t')
colnames(tax_IBD)
tax_IBD=tax_IBD[,c(1,49:304)] #keep ID, Taxa
rownames(tax_IBD)=tax_IBD$SID
tax_IBD$SID=NULL 
rowSums(tax_IBD)   #check if abundances are relative. N.B.: here normalized data!    

*2. LLD-Tax:*                                                                                   
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

*3. Merge Taxonomy IBD and LLD:* 
*N.B.: Metaanalysis will need equal taxa in every cohort. IBD file has more taxa than LLDeep. These 2 files need to be merged by common columns i.e. taxa*   
tax_I=as.data.frame(t(tax_IBD))  #set taxa as row.names 
tax_L=as.data.frame(t(tax_LLD))  #set taxa as row.names
tax=merge(tax_I, tax_L, by = "row.names", all=T)
rownames(tax)=tax$Row.names
tax$Row.names=NULL
tax=as.data.frame(t(tax))
which(is.na(tax))                #NA's are created when merging two df's by columns of different lenght and setting all=TRUE
tax[is.na(tax)] <- 0             #change created NA's to zero 
write.table(tax,'../Metadata/Subsetted Files/Taxonomy/Taxonomy_LLD_&_IBD.txt', sep='\t')

**Subset to relevant columns**  


**2.2 Species**  
*N.B.: Taxonomy files I am using have been filtered (step 5 Medication project) and normalized (step 6 Medication project)  

*All Taxa (n=287)* 
tax2=as.data.frame(t(tax)) 

*1. Species including Strains (n=189):*
Species_strains=tax2[grep('s__', row.names(tax2)),]              
colSums(Species_strains)
Species_strains=as.data.frame(t(Species_strains))
write.table(Species_strains,'../Metadata/Subsetted Files/Species/Species_LLD_&_IBD.txt', sep='\t')

*2. Species without Strains (n=48):* 
#Species_no_strains=Species_Strains[!grepl('t__', row.names(Species_Strains)),]  
#48 Species without strains 

**2.3 Pathways**  
*1. IBD-Path:*  
path_IBD=read.table("../Metadata/IBD_filtered_path_pheno.txt", header=T, sep='\t')
colnames(path_IBD)
path_IBD=path_IBD[,c(1,49:382)] #keep these columns
rownames(path_IBD)=path_IBD$SID
path_IBD$SID=NULL

*2. LLD-Path:*
path_LLD=read.table("../Metadata/LLD_filtered_path_pheno.txt", header=T, sep='\t')
colnames(path_LLD)
path_LLD=path_LLD[,c(1,49:378)] #keep these columns
rownames(path_LLD)=path_LLD$SID
path_LLD$SID=NULL
path_LLD2=as.data.frame(t(path_LLD))  
path_LLD3=path_LLD2[,-grep("_IBS", colnames(path_LLD2))] #delete Maastricht 
path_LLD=as.data.frame(t(path_LLD3))  

*3. Merge Pathway files IBD and LLD:* 
*N.B.: Metaanalysis will need equal taxa in every cohort. The IBD file has more pathways than the LLDeep file. These 2 files need to be merged by common columns i.e. pathways* 

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

*4. Covert into relative abundances i.e. in relation to sum* 
path_rel=t(t(path))/colSums(t(path)) #Or instead of colSums(t(Pathways)) -> rowSums(Pathways)
rowSums(path)
rowSums(path_rel)
#check:for Aerobactinsyn.PWY: 
#IBDFEC0115 in path = 5.93. rowSums(path) = 2099.47
#IBDFEC0115 in path_rel = 0.0028 -> 5.93/2099.47 = 0.0028
path_rel=as.data.frame(path_rel) #safe as data frame
write.table(path_rel,'../Metadata/Subsetted files/Pathways/Pathways_LLD_&_IBD_rel_abundance.txt',sep = '\t')
View(path_rel)

**Growth Rates**
*1. Growth rates IBD and LLD* 
setwd("~/Desktop/Data/Metadata")
growth=read.table("../Metadata/selected_samples_growthrates_2.txt", header=T, sep='\t')
rownames(growth)=growth$X
growth$X=NULL
growth=as.data.frame(t(growth))
growth=growth[,-grep("_IBS", colnames(growth))] #delete Maastricht 
growth=as.data.frame(t(growth))
write.table(growth,'../Metadata/Subsetted Files/Growth Rates/Growth_Rates_LLD_&_IBD_not_rel.txt', sep='\t')


 3.Merge Phenotypic and Microbial Abundance Files   
 -------------
 
 *Maaslin input files require: 
 - saving as tsv in a folder which will be the working directory during the Maaslin runs
 - column names: ID's first, than phenotypes, followed by microbial abundances.                                                   
   In this case: col 5-180 = food, col 180-467 = taxa* 
 
setwd("~/Desktop/Data/Maaslin Files")#CD_Food_Taxonomy

*1. CD-Food-Tax:*
CD_foodtax = merge (CD_Food, tax, by = "row.names", all = FALSE)  
#all = FALSE implies that only rownames that are common in both dataframes will be kept
write.table(CD_foodtax, "CD_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)

*2. UC-Food-Tax:*
UC_foodtax = merge (UC_Food, tax, by = "row.names", all = FALSE)
write.table(UC_foodtax, "UC_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)

*3. IBS_Food_Tax:* 
IBS_foodtax = merge (IBS_Food, tax, by = "row.names", all = FALSE)
write.table(IBS_foodtax, "IBS_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)

*4. HC-Food-Tax:* 
HC_foodtax = merge (HC_Food, tax, by = "row.names", all = FALSE)
write.table(HC_foodtax, "HC_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)

**Merge Food with Species, Pathways and Growth rates, accordingly!**
