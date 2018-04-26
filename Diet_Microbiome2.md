Diet-Microbiome Project (2)
-------------
 
Linear Models for Association Studies in R                                                                                  
Inverse variance based Meta-anlysis                                                                                         
Cluster Analyses

Creator: Laura Bolte

Year: 2018

Last updated: 26 April 2018 


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
 
 **Input files for association analysis and meta-analysis:**
 - save in a folder which will be the working directory  
 - column names: ID's first, than phenotypes, followed by microbial abundances                                                   
   In this case: col 1-4 (ID, age, sex, read depth), col 5-180 (foods), col 180-467 (taxa) 
 
*1. CD-Food-Tax:*
```
setwd("~/Desktop/Data/Input Files") 
CD_foodtax = merge (CD_Food, tax, by = "row.names", all = FALSE)  
#all = FALSE implies that only rownames that are common in both dataframes will be kept
write.table(CD_foodtax, "CD_Food_Tax.tsv", sep= "\t", quote = F, row.names=F)
```
*The same for UC, IBS and HC files*

**Merge Food with Species, Pathways and Growth rates, accordingly!**


4.Association Studies   
 -------------

**Load libraries and files**

*Input files to be saved in working directory, change file to imported (CD/UC/IBS/HC) per test!*
```
library(gamlss)
library(outliers)
setwd("~/Desktop/Data/Association Analyses/Inputfiles")

#Taxonomy
IBD=read.table("../Inputfiles/CD_Food_Tax.txt", sep = "\t", header = T, row.names = 1)
IBD=IBD[,-4]     #data is corrected for kcal, so column kcal can go out 

#Pathways
IBD=read.table("../Inputfiles/CD_Food_Path.txt", sep = "\t", header = T, row.names = 1)
IBD=IBD[,-4]   
colnames(IBD) <- gsub('\\.','_', colnames(IBD))
colnames(IBD) <- gsub('_$', '', colnames(IBD))  #$ = end of line #'^x
```

**QC, iterate over taxa**    
```
setwd("~/Desktop/Data/Association Analyses/")
myOutliers=TRUE
IBD2=IBD
IBD3=IBD2
colnames(IBD)
#Change 179:ncol for the columns containing taxa or pw's 
results_IBD=matrix(,ncol = 9, nrow = length(colnames(IBD2)[179:ncol(IBD2)])) 
x=0                                                                         

#Loop detecting outliers
for (i in 179:ncol(IBD2)){
  x=x+1
  #Get initial statistics of 0 and non-0 per taxa
  results_IBD[x,1]=length(IBD2[,i])
  results_IBD[x,2]=sum(IBD2[,i]!=0)
  #Test outliers, but if all are marked as outliers keep the original values... If you want to perform strict analysis you may want to remove them or pre-filter the table
  while (any(myOutliers==T)){
    if (sum(is.na(IBD2[,i]))==length(IBD2[,i])){
      IBD2[,i]=IBD3[,i]
      break
    }
    outliers=grubbs.test(IBD2[,i])
    myOutliers = outlier( IBD2[,i], logical = TRUE )
    # Threshold for outlier identification
    if (outliers$p.value < 0.05){
      #Transform relative abundances to NA's in case of outliers
      IBD2[,i][myOutliers] <- NA
    } else {
      break
    }
  }
}
```

**Linear Models using the lm function in R** 
-------------
```
#Change 4:178 for the columns containing the phenotypes to test (here 175 foods)  
#Loop per food 
for (a in 4:178){                  
  if (is.na(mean(IBD2[,a]))==T){
    IBD2[,a][is.na(IBD2[,a])] =median(IBD2[,a], na.rm=TRUE)
  }
  results_IBD[,3]=mean(IBD2[,a])    
  results_IBD[,4]=sd(IBD2[,a])     
  food=colnames(IBD2[a])
  results_IBD[,9]=food
  z=0
  
for (b in 179:ncol(IBD2)){
  z=z+1
  if (sum(IBD2[,i]!=0, na.rm = T)!=0){       #sum of non-zeros > 0, means taxon is present in cohort) 
    temp=IBD2[,c(b,1:3,a)] 
    df2<-temp[complete.cases(temp),]
    y=lm(df2[,1]~df2[,5]+Sex+AgeAtFecalSampling+PFReads,data = df2)
    yy=summary(y)
    results_IBD[z,5]=yy$coefficients[2,1]
    results_IBD[z,6]=yy$coefficients[2,2]
    results_IBD[z,7]=yy$coefficients[2,4]
    }
  }
  results_IBD[,8]=p.adjust( results_IBD[,7], method = "fdr")      #adjust for multiple testing (per food for all taxa)
  colnames(results_IBD)=c("N","non-zeros", "Mean", "SD", "Coef", "StdError","pvalue", "qvalue","food")
  rownames(results_IBD)=colnames(IBD2)[179:ncol(IBD2)]
  ## CHANGE PREFIX OF THE RESULT FILE 
  write.table(results_IBD, file=paste('CD',food, sep = '_'), sep="\t", quote=F)
}
```

5.Metaanalysis - Inverse variance based 
 -------------

**Load libraries and input files**
```
library (meta)
library (dplyr) 
#REMOVE how_often_tea_1 and how_often_soda_1 from output folder Association Analyses (double categories) 
setwd("~/Desktop/Data/Association Analyses/Taxonomy/")
setwd("./")
```
**Create Food list** 
*Step 1*
- In Excel: Select column names -> press format -> cells -> custom -> type: \"@\"
- Paste colmun names in R (evt. replace ” by ") 

*Step 2* 
Using the food names as before will give an Error: "invalid 'description' argument.
This is because bread and group_bread OR yoghurt_lf and yoghurt_lf_fruit will match with each other. 
Change file names in folder to e.g. bread_x and yoghurt_lf_x and so forth

*Step 3*
Adapt those food names which needed an X, accordingly in the food list
```
food_list=c("often_breakfast", "often_bread",	"often_hot_meal",	"how_often_alcohol",	"how_often_chocomilk_sweetened_milk_drinks",	"how_often_coffee",	"how_often_fruits",	"how_often_milk_or_buttermilk",	"how_often_muesli",	"how_often_nuts",	"how_often_pasta",	"how_often_pulses",	"how_often_rice",	"how_often_soda",	"how_often_tea",	"how_often_vegetables",	"how_often_yoghurt_milk_based_puddings",	"crackers_x",	"rolls",	"bread_x",	"cheese_20",	"cheese_40",	"cheese_48",	"cheese_other",	"meats_fat",	"meats_other",	"peanutbutter",	"sandwichspread",	"chocolatespreads",	"spreads_sweet",	"egg_cooked",	"egg_baked",	"breakfast_drink",	"cereals_x",	"milk_whole",	"milk_semiskimmed",	"milk_skimmed",	"buttermilk_x",	"chocolatemilk",	"yoghurt_drink_added_sugar",	"yoghurt_drink_other",	"custard_ff",	"yoghurt_ff",	"yoghurt_lf_x",	"yoghurt_lf_fruits",	"fromage_frais_fruits",	"porridge",	"other_custard_yoghurt_fromagefrais",	"icecream_dairy",	"whipped_cream",	"sugar_yoghurt",	"coffee_x",	"sugar_coffee",	"coffeecreamer_hf",	"coffeecreamer_p",	"coffeecreamer_ff",	"milk_coffee",	"tea_x",	"sugar_tea",	"soup_legumes",	"soup_x",	"ready_meals_chinese_indian",	"meals_fast_food",	"ready_meals_other",	"pizza_x",	"pasta_x",	"rice_x",	"legumes_x",	"potato_cooked_mashed",	"potato_baked_fries",	"vegetables_cooked_nobutter",	"herring_salted",	"fish_white_fried",	"fish_lean",	"fish_fatty",	"fish_other",	"meat_x",	"sausage_smoked",	"beef_lean",	"beef_fat",	"pork_lean",	"pork_fat",	"pork_processed",	"chicken",	"other_meat_poultry",	"gravy",	"sauce_hot",	"mayonaise_x",	"non_red_sauces_x",	"saladdressing_f",	"saladdressing_w",	"nut_d",	"cheese_d",	"fruit_x",	"applesauce",	"biscuits_s",	"cake_x",	"pastry_x",	"spiced_cake",	"candybars",	"chocolate_x",	"sweets_x",	"snack_savoury_hot",	"mayonaise_snack",	"non_red_sauces_snack",	"snack_nut",	"crisps_x",	"snack_cheese",	"snack_meats",	"salad_toast",	"softdrink_sugar",	"softdrink_no_sugar",	"fruitjuice", "beer_x",	"beer_af",	"wine_red",	"wine_white",	"wine_fort",	"spirits",	"other_alc_drinks",	"butter_b",	"margarine_lfb",	"butter_ob",	"vegetables_stirfried",	"fish_prepared_fat",	"vegetables_cooked_butter",	"group_alcohol",	"group_breads",	"group_cereals",	"group_cheese",	"group_coffee",	"group_dairy",	"group_eggs",	"group_fish",	"group_fruits",	"group_legumes",	"group_meat",	"group_nonalc_drinks",	"group_nuts",	"group_pasta",	"group_pastry",	"group_potatoes",	"group_prepared_meal",	"group_rice",	"group_sauces",	"group_savoury_snacks",	"group_soup",	"group_spreads",	"group_sugar_sweets",	"group_tea",	"group_vegetables",	"SUMOFEIWITTOT",	"SUMOFEIWITPLANT",	"SUMOFEIWITDIER",	"SUMOFVETTOT",	"SUMOFKHTOT",	"SUMOFALCOHOL",	"Prot_en",	"P_plant_en",	"P_animal_en",	"Fat_en",	"Carb_en",	"Alc_en",	"how_often_boiled_potatos",	"how_often_crisps_savory_crackers",	"how_often_fish",	"how_often_juice",	"how_often_meat",	"how_often_baked_fried_potatoes",	"how_often_chocolate",	"how_often_eggs",	"how_often_icecream",	"how_often_pizza")
```

**Loop to meta-analyse**
```
path="./"
flag=1
for (a in food_list){                     #all in one folder
  file.names = dir(path,pattern=a)
  my_CD=grep("CD",file.names, value = T )
  my_UC=grep("UC",file.names, value = T )
  my_IBS=grep("IBS",file.names, value = T )
  my_HC=grep("HC",file.names, value = T )
  CD=read.table(my_CD, sep="\t", header = T, row.names = 1)
  UC=read.table(my_UC, sep="\t", header = T, row.names = 1)
  IBS=read.table(my_IBS, sep="\t", header = T, row.names = 1)
  HC=read.table(my_HC, sep="\t", header = T, row.names = 1)
  
  CD_UC=merge(CD, UC, by = "row.names")
  rownames(CD_UC)=CD_UC$Row.names
  CD_UC$Row.names=NULL
  
  IBS_HC=merge(IBS, HC, by = "row.names")
  rownames(IBS_HC)=IBS_HC$Row.names
  IBS_HC$Row.names=NULL
  
  all=merge(CD_UC, IBS_HC, by = "row.names")
  
  list_coef=all[,c(1,6,15,24,33)] 
  list_coef[list_coef == 0] <- NA
  list_coef2=na.omit(list_coef) #coef zero means bacteria are not present in cohort -> remove rows if coef is zero 
  
  #selection=all[list_coef2$Row.names]
  selection = all %>% semi_join(list_coef2, by = "Row.names")
  #selection=all
  colnames(selection)=c("Taxa", 
                        "N.CD", "non-zeros.CD", "Mean.CD", "SD.CD", "Coef.CD", "SE.CD", "Pval.CD", "Qval.CD", "food.CD",
                        "N.UC", "non-zeros.UC", "Mean.UC", "SD.UC", "Coef.UC", "SE.UC", "Pval.UC", "Qval.UC", "food.UC",
                        "N.IBS","non-zeros.IBS","Mean.IBS","SD.IBS","Coef.IBS","SE.IBS","Pval.IBS","Qval.IBS","food.IBS", 
                        "N.HC", "non-zeros.HC", "Mean.HC", "SD.HC", "Coef.HC", "SE.HC", "Pval.HC", "Qval.HC", "food.HC")
  #9*4 columns = 36 + 1 (taxa column) = 37  
  
  #Calculate Inverse variance
  #38
  selection$inverse_var.cd=1/selection$SE.CD^2     #METAL Paper 
  #39
  selection$inverse_var.uc=1/selection$SE.UC^2
  #40
  selection$inverse_var.ibs=1/selection$SE.IBS^2   
  #41
  selection$inverse_var.hc=1/selection$SE.HC^2 
  
  #Calculate SE  #42
  selection$se=sqrt(1/
                      (selection$inverse_var.cd+
                         selection$inverse_var.uc+
                         selection$inverse_var.ibs+
                         selection$inverse_var.hc))
  #Calculate Beta #43
  selection$beta=(selection$inverse_var.cd*selection$Coef.CD+
                    selection$inverse_var.uc*selection$Coef.UC+
                    selection$inverse_var.ibs*selection$Coef.IBS+
                    selection$inverse_var.hc*selection$Coef.HC)/
    (selection$inverse_var.cd+
       selection$inverse_var.uc+
       selection$inverse_var.ibs+
       selection$inverse_var.hc)
  #selection$beta=   (selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.mibs*selection$coeff.correcting.all.MIBS+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld)
  
  #Calculate Z-score #44
  selection$Z=selection$beta/selection$se
  
  #Calculate meta p-value #45
  selection$P=2*pnorm(-abs(selection$Z))
  
  #Adjust pvalue with FDR #46
  selection$FDR=p.adjust(selection$P,method = "fdr")   #we are looping for each file, so its per food group for all taxa 
  
  #Create empty columns
  selection$Het.Q="NS" #47
  selection$Het.I2="NS" #48
  selection$Het.Pval="NS" #49
  ```
  
 6.Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
 -------------
  ```
  for (i in 1:length(rownames(selection))){
    if (selection$FDR[i]<0.1){
      #Select coefficients
      TE=c( selection[i,6], selection[i,15], selection[i,24], selection[i,33])  #select cols of coefficients e.g. coef CD= 5
      #Select Standart error
      SE=c( selection[i,7], selection[i,16], selection[i,25], selection[i,34] ) #select cols of standard error e.g. coef CD= 6
      het=metagen(TE,SE)
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      selection[i,48]=het$I2   
      #Match Het.Q column number
      selection[i,47]=het$Q 
      #Calculate p-value from Q calculation // Het.Pval
      selection[i,49]=pchisq(het$Q,df=2,lower.tail=F)  #het q (Het.Pval)
    }
  }
  write.table(selection, file=paste(a, "_meta.txt", sep=""), sep="\t", quote=F)
}
```

7.Merge all results in BASH
 ------------- 
 
**All results incl. non-significant meta-results**

*Results that were significant in one cohort but not in the metaanalysis, are also important to check
in order to spot differences between cohorts*

```
cd /Users/laurabolte/Desktop/Data/Association\ Analyses/Taxonomy/meta 
mv *meta.txt ./meta/
less group_breads_meta.txt | head -1 >> ../tax_all_results.txt
for i in *.txt ; do  less $i >> ../tax_all_results.txt ; done
#Remove repeated row headers:
In Excel: shift first header to right, remove first column, save, import into R.
Tax=read.table("Data/Association Analyses/Taxonomy/tax_all_results.txt", sep='\t', header=T)    
Tax1=Tax[!grepl("N.CD",Tax$Taxa),]                                                            #172 header rows removed
write.table(Tax1,"../Desktop/Data/Association Analyses/tax_all_results.txt", sep='\t')
```

**Only significant results** 
```
cd /Users/laurabolte/Desktop/Data/Association\ Analyses/Taxonomy/meta 
mv *meta.txt ./meta/
ls
less group_breads_meta.txt | head -1 >> ../tax_results_sign.txt 
for i in *.txt ; do  less $i | awk -F "\t" '{if ($47<0.1){print $0}}' >> ../tax_results_sign.txt; done
```

8.Cluster Analyses
 -------------
 
**8.1 Subset**

As previously described. Input files: rows: matched IDs, columns: foods

**8.2 CLUSTERING DIET** 

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

**8.3 CLUSTERING SPECIES** 

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

**8.4. CORRELATING FOOD CLUSTERS WITH SPECIES CLUSTERS**

*1. Load Input Files* 
```
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
Consider that in the cluster analysis part, cutting the tree lower (lower h=), will create more clusters and stronger correlation

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
sum(my_results_fdr <= 0.05)    
```

*8. Corrplot* 
```
dev.off()
#Color
display.brewer.all() #display all colour schemes
col1 = rev(brewer.pal(n=6, name="RdYlBu")) #reverse colors so that positive association = red, negative = blue)
color = c(col1,col1,col1) #repeat col to extend the range
#cl.pos="b"               #plots legenda below plot instead of right side

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
