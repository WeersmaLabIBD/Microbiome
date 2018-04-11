####################################################################################################################
#------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------cookbook metagenomics-----------------------------------------#
#------------------------------------------------------------------------------------------------------------------#
####################################################################################################################

# this is the most convincible script
# based on Alex's pipeline with a little change
# sample numbers and names are double-checked, so no worries any more

###---------------------------IBD----------------------------------------------------------------------------------

# 1454 taxa and 544 samples
ibd_meta=read.table(file = "IBD_taxonomy_metaphlan2_082017.txt",sep = "\t",header = T,check.names = F)

ibd_meta_ID=as.data.frame(unlist(colnames(ibd_meta)[-1]),stringsAsFactors = F)
colnames(ibd_meta_ID)[1]="Original"

# check ID
# Dianne file: from G****_metaphlan to UMCG****
# But in Dianne file, 12 samples have only 6 UMCG**** names(each 2 have the same UMCG names), which means now the number of samples is 538
old=read.table("IBD_Dienna.txt",sep="\t",header = T)
old=old[!duplicated(old$ID), ]
old=old[!duplicated(old$Original), ]

# new.txt: from UMCG**** to 214-2637 (for example)
new=read.table("IBD_new.txt",sep = "\t",header = T)
new=new[!duplicated(new$ID), ]

# 12 samples with metagenomics don't have exome data, so now the number of samples is 526
# duplicated names in UMCG****
id_change=merge(old,new,by="ID")
id_change=merge(id_change,ibd_meta_ID,by="Original")

# remove ancestry outliers
outliers=read.table("IBD_ancenstry_outliers.txt",header = T,sep = "\t")
colnames(outliers)[1]="Sample_ID"
id_change=id_change[!id_change$Sample_ID %in% outliers$Sample_ID,]

# extract final 504 samples with both exome and metagenomics data
ibd_sub_meta=ibd_meta[,colnames(ibd_meta) %in% id_change$Original]

rownames(ibd_sub_meta)=ibd_meta$ID
id_change=id_change[id_change$Original %in% colnames(ibd_sub_meta),]
id_change=id_change[!duplicated(id_change$ID), ]

# create pheno file , remove any special markers
pheno=read.table("eqtl_IBD_pheno.txt",sep="\t",header = T)
pheno=pheno[-1,]
colnames(pheno)[1:2]=c("ID","Sample_ID")
pheno=pheno[!duplicated(pheno$Sample_ID), ]
pheno_sub=merge(id_change,pheno,by="ID")
pheno_sub=pheno_sub[!duplicated(pheno_sub$ID), ]

# add depth info
depth=read.table("IBD_reads_summary.txt",sep = "\t",header = T)
depth$Original=paste(depth$Sample_ID,"metaphlan",sep = "_")
pheno_sub=merge(pheno_sub,depth,by="Original")

covariate=pheno_sub[,c(4,6,8,13)]
colnames(covariate)[1]="ID"

covariate=covariate[!duplicated(covariate$ID), ]
rownames(covariate)=covariate[,1]
covariate=covariate[,-1]

write.table(covariate,sep = "\t",file = "eqtl_IBD_covariate.txt",quote = F)

# cut off based on presence among samples

cutoff=0.1

taxa=ibd_sub_meta

mm = taxa[rowSums(taxa > 0) >= cutoff*ncol(taxa),]
mm[mm == 0] = NA
mm = log(mm)

write.table(t(mm),file = "IBD_filtered_logTrans.txt",sep = "\t",quote = F)

bb=mm
bb[is.na(bb)]=0

# generate coupling

coupling=id_change[,c(4,4)]

write.table(coupling,file = "coupling_file.txt",row.names = F,col.names = F,quote = F,sep = "\t")


# generate linkage used in correcting step

coupling=id_change[,c(4,1)]
write.table(coupling,file = "eqtl_IBD_linkage.txt",row.names = F,col.names = F,quote = F,sep = "\t")

# correct and split

tax=t(mm)
coupling = read.table("eqtl_IBD_linkage.txt",colClasses = "character",sep = "\t",check.names = F)

has_both = (coupling[,1] %in% rownames(covariate)) & (coupling[,2] %in% rownames(tax))
coupling= coupling[has_both,]
tax = tax[coupling[,2],]
covariate = covariate[coupling[,1],,drop = FALSE]
rownames(tax) = rownames(covariate)
covariate = covariate[rownames(tax),,drop = FALSE]

covariate$Gender=as.character(covariate$Gender)
covariate$Gender[which(covariate$Gender=="Female")]=0
covariate$Gender[which(covariate$Gender=="Male")]=1

corrected_data = apply(tax,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  if(ncol(covariate.subset)==ncol(covariate.subset.matrix)){
    covariate.subset = covariate.subset[,apply(covariate.subset.matrix,2,sd) !=0,drop = FALSE]
  }
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

#correct quantitative data for covariates
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

#generate presence/absence tables
binary_data = apply(tax,2,function(x){as.integer(!is.na(x))})
dimnames(binary_data) = dimnames(tax)
binary_data = binary_data[,colSums(binary_data)>0.1*nrow(binary_data)&colSums(binary_data)<0.9*nrow(binary_data)]
binary_data = as.data.frame(t(binary_data))+100
binary_data = cbind(paste0(rownames(binary_data),".binary"),binary_data)

colnames(binary_data)[1] = "-"
binary_annot = data.frame(platform = "RDP",
                          HT12v3.ArrayAddress = binary_data[,1],
                          Gene = binary_data[,1],
                          Chr = 4,
                          ChrStart = 1000,
                          ChrEnd = 1000)

write.table(corrected_data2, file = "tax_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "tax_numeric.txt.annot",sep="\t",row.names=F,quote = F)

write.table(binary_data, file = "tax_binary.txt",sep="\t",row.names = F,quote = F)
write.table(binary_annot,file = "tax_binary.txt.annot",sep="\t",row.names=F,quote = F)


###--------------------------------------------------------LLD-----------------------------------------------------


lld_meta=read.table(file = "LLD_taxonomy_metaphlan2_092017.txt",sep = "\t",header = T,check.names = F)

lld_pheno=read.table(file = "eqtl_LLD_pheno.txt",as.is = T,header = T,check.names = F)
colnames(lld_pheno)[1]="exome"

# check ID
new=read.table("LLD_ID_list.txt",sep="\t",header = F)
colnames(new)="exome"
old=read.table("LLD_metaphlan_list.txt",sep = "\t",header = F,check.names = F)
change=read.table("rename_LLD.txt",sep = "\t",header = T,check.names = F)

depth=read.table("depth_LLD.txt",sep = "\t",header = T)
colnames(depth)[1]="exome"

change$metaphlan=paste(change$ID,"metaphlan",sep = "_")
colnames(old)="metaphlan"

compare=merge(change,old,by="metaphlan")
colnames(compare)[3]="exome"
intersect=merge(compare,new,by="exome")
intersect=merge(intersect,depth,by="exome")

outliers=read.table("LLD_ancenstry_outliers.txt",sep = "\t",header = T)
intersect=intersect[!intersect$exome %in% outliers$IID,]

pheno=merge(intersect,lld_pheno,by="exome")[,c(1,5,7,9)]
rownames(pheno)=pheno[,1]

covariate=merge(pheno,intersect,by="exome")
covariate=covariate[,1:4]
rownames(covariate)=covariate$exome
covariate=covariate[,-1]

write.table(covariate,file = "eqtl_lld_covariate.txt" ,sep = "\t",quote = F)

# cut off based on presence among samples

lld_sub_meta=lld_meta[,colnames(lld_meta) %in% intersect$metaphlan]

rownames(lld_sub_meta)=lld_meta$ID

cutoff=0.1

taxa=lld_sub_meta

mm = taxa[rowSums(taxa > 0) >= cutoff*ncol(taxa),]
mm[mm == 0] = NA
mm = log(mm)

write.table(t(mm),file = "LLD_filtered_logTrans.txt",sep = "\t",quote = F)

# generate coupling

coupling=intersect[,c(1,1)]

write.table(coupling,file = "coupling_file.txt",row.names = F,col.names = F,sep = "\t",quote = F)
write.table(intersect[,1:2],"eqtl_LLD_linkage.txt",row.names = F,col.names = F,sep = "\t",quote = F)

# correct and split

tax=t(mm)

coupling = read.table("eqtl_LLD_linkage.txt",colClasses = "character",sep = "\t")

has_both = (coupling[,1] %in% rownames(covariate)) & (coupling[,2] %in% rownames(tax))
coupling= coupling[has_both,]
tax = tax[coupling[,2],]
covariate = covariate[coupling[,1],,drop = FALSE]
rownames(tax)=rownames(covariate)

covariate$Gender[which(covariate$Gender=="Female")]=0
covariate$Gender[which(covariate$Gender=="Male")]=1

corrected_data = apply(tax,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  if(ncol(covariate.subset)==ncol(covariate.subset.matrix)){
    covariate.subset = covariate.subset[,apply(covariate.subset.matrix,2,sd) !=0,drop = FALSE]
  }
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

#correct quantitative data for covariates
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

#generate presence/absence tables
binary_data = apply(tax,2,function(x){as.integer(!is.na(x))})
dimnames(binary_data) = dimnames(tax)
binary_data = binary_data[,colSums(binary_data)>0.1*nrow(binary_data)&colSums(binary_data)<0.9*nrow(binary_data)]
binary_data = as.data.frame(t(binary_data))+100
binary_data = cbind(paste0(rownames(binary_data),".binary"),binary_data)

colnames(binary_data)[1] = "-"
binary_annot = data.frame(platform = "RDP",
                          HT12v3.ArrayAddress = binary_data[,1],
                          Gene = binary_data[,1],
                          Chr = 4,
                          ChrStart = 1000,
                          ChrEnd = 1000)

write.table(corrected_data2, file = "tax_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "tax_numeric.txt.annot",sep="\t",row.names=F,quote = F)

write.table(binary_data, file = "tax_binary.txt",sep="\t",row.names = F,quote = F)
write.table(binary_annot,file = "tax_binary.txt.annot",sep="\t",row.names=F,quote = F)


