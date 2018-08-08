
#######################################################################################################
############################# metagenomics data correction for three covariates ####################
#######################################################################################################

########################### IBD #######################################################################

ibd_id_change=read.table("IBD_allID_change.txt",header = T,stringsAsFactors = F,sep = "\t")

################## covariate ####################

ibd_covariate=read.table("IBD_allCovariates_Arnau.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
ibd_covariate=ibd_covariate[ibd_covariate$UMCGIBDResearchIDorLLDeepIDorMIBSID %in% ibd_id_change$ID,]
ibd_covariate=ibd_covariate[,c(2,1,4:11)]
colnames(ibd_covariate)[1]="ID"
ibd_covariate=merge(ibd_covariate,ibd_id_change,all=T,by="ID")[,c(1:10,13)]
ibd_covariate[ibd_covariate$Sample_ID=="209-6495",c(2,10)]=c("female",25)
ibd_covariate[ibd_covariate$Sample_ID=="214-1341",c(2,10)]=c("female",48)
ibd_covariate[ibd_covariate$Sample_ID=="214-1709",c(2,10)]=c("female",37)
rownames(ibd_covariate)=ibd_covariate[,11]
ibd_covariate=ibd_covariate[,c(-1,-11)]
ibd_covariate=ibd_covariate[order(rownames(ibd_covariate)),]
ibd_covariate$AgeAtFecalSampling=as.numeric(ibd_covariate$AgeAtFecalSampling)

#### impute missing NA values
## Solution 1 From Alex If the proportion of NAs is not large, say less than 5%, you can replace NAs by variable mean or median. We usually use median
## Solution 2 From lmchen using imputePCA (not very clear about the principle)
ibd_covariate_imputated=ibd_covariate
ibd_covariate_imputated$Sex=factor(ibd_covariate_imputated$Sex, levels=c("male","female"), labels=c(0,1))
ibd_covariate_imputated$Sex= as.integer(as.character(ibd_covariate_imputated$Sex))

ibd_covariate_imputated$DiseaseLocation=factor(ibd_covariate_imputated$DiseaseLocation, levels=c("colon","ileum","both"), labels=c(0,1,2))
ibd_covariate_imputated$DiseaseLocation= as.integer(as.character(ibd_covariate_imputated$DiseaseLocation))

ibd_covariate_imputated$MedicationAntibiotics=factor(ibd_covariate_imputated$MedicationAntibiotics, levels=c("no","yes"), labels=c(0,1))
ibd_covariate_imputated$MedicationAntibiotics= as.integer(as.character(ibd_covariate_imputated$MedicationAntibiotics))

ibd_covariate_imputated$MedicationPPI=factor(ibd_covariate_imputated$MedicationPPI, levels=c("no","yes"), labels=c(0,1))
ibd_covariate_imputated$MedicationPPI= as.integer(as.character(ibd_covariate_imputated$MedicationPPI))

ibd_covariate_imputated$laxatives=factor(ibd_covariate_imputated$laxatives, levels=c("no","yes"), labels=c(0,1))
ibd_covariate_imputated$laxatives= as.integer(as.character(ibd_covariate_imputated$laxatives))

ibd_covariate_imputated$SmokeAreYouAFormerSmoker=factor(ibd_covariate_imputated$SmokeAreYouAFormerSmoker, levels=c("no","yes"), labels=c(0,1))
ibd_covariate_imputated$SmokeAreYouAFormerSmoker= as.integer(as.character(ibd_covariate_imputated$SmokeAreYouAFormerSmoker))

ibd_covariate_imputated=apply(ibd_covariate_imputated,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
ibd_covariate_imputated=as.data.frame(ibd_covariate_imputated)
ibd_covariate_imputated=ibd_covariate_imputated[,c(1,2,9)]
write.table(ibd_covariate_imputated,file = "IBD_imputated_covariate.txt",sep = "\t",quote = F,row.names = T)

################## pathway data ####################

#### extract samples based on ID
library(stringr)
ibd_pathway=read.table("IBD_humann2_pathways_uniref90_082017.txt",header = T,stringsAsFactors = F,check.names = F,sep = "\t",comment.char = "",row.names = 1)
ibd_path=ibd_pathway[grep("__",rownames(ibd_pathway),invert=T),]
rownames(ibd_path)=str_replace_all(rownames(ibd_path)," ","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),"/","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),":","")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),",","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),"\\|","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),"\\(","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),"\\)","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),"&","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),";","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),"\\[","_")
rownames(ibd_path)=str_replace_all(rownames(ibd_path),"]","_")
ibd_path_ID=as.data.frame(unlist(colnames(ibd_path)),stringsAsFactors = F)
colnames(ibd_path_ID)="path_ID"
ibd_id_change$path_ID = sapply(strsplit(as.character(ibd_id_change$Original),'_'), "[", 1)
ibd_id_change$path_ID=paste(ibd_id_change$path_ID,"_kneaddata_merged_Abundance",sep = "")
ibd_pathway=ibd_path[,ibd_path_ID$path_ID %in% ibd_id_change$path_ID]

#### caculate relative abundance
ibd_pathway = ibd_pathway[rowSums(ibd_pathway > 0) > ncol(ibd_pathway)*0.25,]
ibd_pathway = sweep(ibd_pathway,2,colSums(ibd_pathway),"/")
ibd_pathway[ibd_pathway == 0] = NA
ibd_pathway = log(ibd_pathway)
ibd_pathway=t(ibd_pathway)
ibd_pathway=as.data.frame(ibd_pathway,stringsAsFactors = F)
rownames(ibd_pathway)=ibd_id_change$Sample_ID
ibd_pathway=ibd_pathway[order(rownames(ibd_pathway)),]

#### linear model correction
ibd_pathways_correct = apply(ibd_pathway,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = ibd_covariate_imputated[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

################## taxa data ####################

ibd_taxa=read.table("IBD_taxonomy_metaphlan2_082017.txt",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
ibd_taxa=ibd_taxa[,colnames(ibd_taxa) %in% ibd_id_change$Original]
colnames(ibd_taxa)=ibd_id_change$Sample_ID

#### cut off based on present rate among samples and arcsin transformation
cutoff=0.1
ibd_taxa = ibd_taxa[rowSums(ibd_taxa > 0) >= cutoff*ncol(ibd_taxa),]
ibd_taxa[ibd_taxa == 0] = NA
ibd_taxa=ibd_taxa/100
ibd_taxa=sqrt(ibd_taxa)
ibd_taxa=asin(ibd_taxa)

#### clean taxa. If multiple levels have the same abundance, keep the lowest level one.
taxa_clean=list()
for(i in 1:ncol(ibd_taxa)){
  
  x.clean=unlist(rownames(ibd_taxa)[!duplicated(ibd_taxa[,i],fromLast = T)])
  taxa_clean=union(taxa_clean,x.clean)
  
}
ibd_taxa_clean=ibd_taxa[rownames(ibd_taxa) %in% taxa_clean,]
ibd_taxa_clean=as.data.frame(t(ibd_taxa_clean))
ibd_taxa_clean=ibd_taxa_clean[order(rownames(ibd_taxa_clean)),]

#### linear model correction
ibd_taxa_correct = apply(ibd_taxa_clean,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = ibd_covariate_imputated[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  if(ncol(covariate.subset)==ncol(covariate.subset.matrix)){
    covariate.subset = covariate.subset[,apply(covariate.subset.matrix,2,sd) !=0,drop = FALSE]
  }
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

corrected_data=cbind(ibd_taxa_correct,ibd_pathways_correct)
corrected_data = as.data.frame(t(corrected_data))
corrected_data2 = cbind(rownames(corrected_data),corrected_data)
colnames(corrected_data2)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_data2),
                   Gene = rownames(corrected_data2),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_data2, file = "IBD_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "IBD_numeric.txt.annot",sep="\t",row.names=F,quote = F)
write.table(ibd_id_change[,c(4,4)],file="IBD_coupling_file.txt",row.names = F,quote = F,sep = "\t",col.names = F)




########################### LLD #######################################################################

lld_id_change=read.table("LLD_allID_change.txt",header = T,stringsAsFactors = F,sep = "\t")

################## covariate ####################

lld_covariate=read.table("LLD_allCovariates_Arnau.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
lld_covariate=lld_covariate[lld_covariate$SID %in% lld_id_change$exome,]
rownames(lld_covariate)=lld_covariate[,1]
lld_covariate=lld_covariate[,-1]
lld_covariate=lld_covariate[order(rownames(lld_covariate)),]

#### impute missing NA values
## Solution 1 From Alex If the proportion of NAs is not large, say less than 5%, you can replace NAs by variable mean or median. We usually use median
## Solution 2 From lmchen using imputePCA (not very clear about the principle)
lld_covariate_imputated=lld_covariate
lld_covariate_imputated$Sex=factor(lld_covariate_imputated$Sex, levels=c("2","1"), labels=c(0,1))
lld_covariate_imputated$Sex= as.integer(as.character(lld_covariate_imputated$Sex))

lld_covariate_imputated$DiagnosisCurrent=factor(lld_covariate_imputated$DiagnosisCurrent, levels=c("1_HC","2_IBS"), labels=c(0,1))
lld_covariate_imputated$DiagnosisCurrent= as.integer(as.character(lld_covariate_imputated$DiagnosisCurrent))

lld_covariate_imputated=apply(lld_covariate_imputated,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
lld_covariate_imputated=as.data.frame(lld_covariate_imputated)
lld_covariate_imputated=lld_covariate_imputated[,c(2,3,4)]
write.table(lld_covariate_imputated,file = "LLD_imputated_covariate.txt",sep = "\t",quote = F,row.names = T)

################## pathway data ####################

#### extract samples based on ID
library(stringr)
lld_pathway=read.table("LLD_humann2_Uniref90_092017.txt",header = T,stringsAsFactors = F,check.names = F,sep = "\t",comment.char = "",row.names = 1)
lld_path=lld_pathway[grep("__",rownames(lld_pathway),invert=T),]
rownames(lld_path)=str_replace_all(rownames(lld_path)," ","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),"/","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),":","")
rownames(lld_path)=str_replace_all(rownames(lld_path),",","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),"\\|","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),"\\(","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),"\\)","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),"&","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),";","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),"\\[","_")
rownames(lld_path)=str_replace_all(rownames(lld_path),"]","_")
lld_id_change$path_ID = sapply(strsplit(as.character(lld_id_change$ID),'_'), "[", 1)
lld_id_change$path_ID=paste(lld_id_change$path_ID,"_kneaddata_merged_Abundance",sep = "")
lld_pathway=lld_path[,colnames(lld_path) %in% lld_id_change$path_ID]
colnames(lld_pathway)=lld_id_change$exome

#### caculate relative abundance
lld_pathway = lld_pathway[rowSums(lld_pathway > 0) > ncol(lld_pathway)*0.25,]
lld_pathway = sweep(lld_pathway,2,colSums(lld_pathway),"/")
lld_pathway[lld_pathway == 0] = NA
lld_pathway = log(lld_pathway)
lld_pathway=t(lld_pathway)
lld_pathway=as.data.frame(lld_pathway,stringsAsFactors = F)
lld_pathway=lld_pathway[order(rownames(lld_pathway)),]

#### linear model correction
lld_pathways_correct = apply(lld_pathway,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = lld_covariate_imputated[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

################## taxa data ####################

lld_taxa=read.table("LLD_taxonomy_metaphlan2_092017.txt",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
lld_taxa=lld_taxa[,colnames(lld_taxa) %in% lld_id_change$metaphlan]
colnames(lld_taxa)=lld_id_change$exome

#### cut off based on present rate among samples and arcsin transformation
cutoff=0.1
lld_taxa = lld_taxa[rowSums(lld_taxa > 0) >= cutoff*ncol(lld_taxa),]
lld_taxa[lld_taxa == 0] = NA
lld_taxa=lld_taxa/100
lld_taxa=sqrt(lld_taxa)
lld_taxa=asin(lld_taxa)

#### clean taxa. If multiple levels have the same abundance, keep the lowest level one.
taxa_clean=list()
for(i in 1:ncol(lld_taxa)){
  
  x.clean=unlist(rownames(lld_taxa)[!duplicated(lld_taxa[,i],fromLast = T)])
  taxa_clean=union(taxa_clean,x.clean)
  
}
lld_taxa_clean=lld_taxa[rownames(lld_taxa) %in% taxa_clean,]
lld_taxa_clean=as.data.frame(t(lld_taxa_clean))
lld_taxa_clean=lld_taxa_clean[order(rownames(lld_taxa_clean)),]

#### linear model correction
lld_taxa_correct = apply(lld_taxa_clean,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = lld_covariate_imputated[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  if(ncol(covariate.subset)==ncol(covariate.subset.matrix)){
    covariate.subset = covariate.subset[,apply(covariate.subset.matrix,2,sd) !=0,drop = FALSE]
  }
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

corrected_data=cbind(lld_taxa_correct,lld_pathways_correct)
corrected_data = as.data.frame(t(corrected_data))
corrected_data2 = cbind(rownames(corrected_data),corrected_data)
colnames(corrected_data2)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_data2),
                   Gene = rownames(corrected_data2),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_data2, file = "LLD_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "LLD_numeric.txt.annot",sep="\t",row.names=F,quote = F)
write.table(lld_id_change[,c(2,2)],file="LLD_coupling_file.txt",row.names = F,quote = F,sep = "\t",col.names = F)
