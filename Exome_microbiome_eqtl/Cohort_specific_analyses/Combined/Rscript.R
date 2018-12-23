#====================================================
#               Prepare genotype and probe
#====================================================

ibd_genotype=read.table("cohort_specific.IBD.dosage.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
lld_genotype=read.table("cohort_specific.LLD.dosage.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
ibd_genotype_target=read.table("cohort_specific.IBD.dosage.target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
lld_genotype_target=read.table("cohort_specific.LLD.dosage.target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
rm=read.table("IBD_ID_change_rm_SB.txt",header = F,sep = "\t",stringsAsFactors = F)
snp=ibd_genotype[,1,drop=F]
ibd_genotype=ibd_genotype[,colnames(ibd_genotype) %in% rm$V2]
ibd_genotype=cbind(snp,ibd_genotype)
snp=ibd_genotype_target[,1,drop=F]
ibd_genotype_target=ibd_genotype_target[,colnames(ibd_genotype_target) %in% rm$V2]
ibd_genotype_target=cbind(snp,ibd_genotype_target)

ibd_genotype=rbind(ibd_genotype,ibd_genotype_target)
lld_genotype=rbind(lld_genotype,lld_genotype_target)
ibd_genotype=unique(ibd_genotype)
lld_genotype=unique(lld_genotype)

rownames(ibd_genotype)=ibd_genotype[,1]
ibd_genotype=ibd_genotype[,-1]
ibd_genotype=as.data.frame(t(ibd_genotype))

rownames(lld_genotype)=lld_genotype[,1]
lld_genotype=lld_genotype[,-1]
lld_genotype=as.data.frame(t(lld_genotype))

ibd_geno_new=matrix(ncol = ncol(ibd_genotype),nrow = nrow(ibd_genotype))
for(i in 1:ncol(ibd_genotype)){
  
  sum_0=sum(ibd_genotype[,i]==0)
  sum_1=sum(ibd_genotype[,i]==1)
  sum_2=sum(ibd_genotype[,i]==2)
  maf=min(sum_0,sum_1,sum_2)
  if(sum_0==maf){
    ibd_geno_new[ibd_genotype[,i]==0,i]=2
    ibd_geno_new[ibd_genotype[,i]==2,i]=0
    ibd_geno_new[ibd_genotype[,i]==1,i]=1
  }else{
    ibd_geno_new[,i]=ibd_genotype[,i]
  }
}
ibd_geno_new=as.data.frame(ibd_geno_new,row.names = rownames(ibd_genotype))
colnames(ibd_geno_new)=colnames(ibd_genotype)

lld_geno_new=matrix(ncol = ncol(lld_genotype),nrow = nrow(lld_genotype))
for(i in 1:ncol(lld_genotype)){
  
  sum_0=sum(lld_genotype[,i]==0)
  sum_1=sum(lld_genotype[,i]==1)
  sum_2=sum(lld_genotype[,i]==2)
  maf=min(sum_0,sum_1,sum_2)
  if(sum_0==maf){
    lld_geno_new[lld_genotype[,i]==0,i]=2
    lld_geno_new[lld_genotype[,i]==2,i]=0
    lld_geno_new[lld_genotype[,i]==1,i]=1
  }else{
    lld_geno_new[,i]=lld_genotype[,i]
  }
}
lld_geno_new=as.data.frame(lld_geno_new,row.names = rownames(lld_genotype))
colnames(lld_geno_new)=colnames(lld_genotype)

ibd_probe=read.table("IBD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
ibd_probe=as.data.frame(t(ibd_probe))
ibd_probe[ibd_probe==0]=NA
lld_probe=read.table("LLD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
lld_probe=as.data.frame(t(lld_probe))
lld_probe[lld_probe==0]=NA

ibd_specific_1=read.table("IBD_specific.txt",sep = "\t",header = T,stringsAsFactors = F)
ibd_specific_2=read.table("ibd_specific_target.txt",sep = "\t",header = T,stringsAsFactors = F)
ibd_specific=rbind(ibd_specific_1,ibd_specific_2)
ibd_specific=unique(ibd_specific)
lld_specific_1=read.table("LLD_specific.txt",sep = "\t",header = T,stringsAsFactors = F)
lld_specific_2=read.table("lld_specific_target.txt",sep = "\t",header = T,stringsAsFactors = F)
lld_specific=rbind(lld_specific_1,lld_specific_2)
lld_specific=unique(lld_specific)
specific=rbind(ibd_specific,lld_specific)

remove=specific[grep("k__",specific$ProbeName,invert = T),]
remove=remove[grep("unclassified",remove$ProbeName),]
specific=specific[!specific$ProbeName %in% remove$ProbeName,]

#==========================================================
#                  WES permutation 1000
#==========================================================

specific_wes=rbind(ibd_specific_1,lld_specific_1)
remove=specific_wes[grep("k__",specific_wes$ProbeName,invert = T),]
remove=remove[grep("unclassified",remove$ProbeName),]
specific_wes=specific_wes[!specific_wes$ProbeName %in% remove$ProbeName,]


ibd_geno_new=ibd_geno_new[,colnames(ibd_geno_new) %in% specific_wes$SNP]
ibd_probe=ibd_probe[,colnames(ibd_probe) %in% specific_wes$ProbeName]
ibd_probe=ibd_probe[rownames(ibd_probe) %in% rownames(ibd_geno_new),]
ibd_probe=ibd_probe[order(rownames(ibd_probe)),]
ibd_geno_new=ibd_geno_new[order(rownames(ibd_geno_new)),]
specific_wes=specific_wes[specific_wes$ProbeName %in% colnames(ibd_probe),]

lld_geno_new=lld_geno_new[,colnames(lld_geno_new) %in% specific_wes$SNP]
lld_probe=lld_probe[,colnames(lld_probe) %in% specific_wes$ProbeName]
lld_probe=lld_probe[rownames(lld_probe) %in% rownames(lld_geno_new),]
lld_probe=lld_probe[order(rownames(lld_probe)),]
lld_geno_new=lld_geno_new[order(rownames(lld_geno_new)),]
specific_wes=specific_wes[specific_wes$ProbeName %in% colnames(lld_probe),]

library(foreach)
library(iterators)
library(RVAideMemoire)

probe_wes=rbind(ibd_probe,lld_probe)
genotype_wes=rbind(ibd_geno_new,lld_geno_new)
disease_wes=genotype_wes[,1:2,drop=F]
colnames(disease_wes)[1]="Status"
disease_wes$Status=as.character(disease_wes$Status)
disease_wes$Status[1:435]=1
disease_wes$Status[436:nrow(disease_wes)]=0
disease_wes[,-1]=NULL

# real one
interaction_wes = foreach(i=1:nrow(specific_wes),.combine = rbind) %do%  {
  probe.sub=as.data.frame(probe_wes[,specific_wes$ProbeName[i]],drop=F,row.names = rownames(probe_wes))
  probe.sub=na.omit(probe.sub)
  snp.sub=specific_wes$SNP[i]
  print(specific_wes[i,])
  genotype.sub=genotype_wes[rownames(genotype_wes) %in% rownames(probe.sub),]
  genotype.sub=genotype.sub[,snp.sub,drop=F]
  disease_wes.sub=disease_wes[rownames(disease_wes) %in% rownames(probe.sub),,drop=F]
  lm1 = lm(unlist(probe.sub) ~ as.numeric(disease_wes.sub$Status) * as.numeric(genotype.sub[,1]))
  summary1 = summary(lm1)$coef
  return.string = data.frame(ProbeName = specific_wes$ProbeName[i], SNPName = specific_wes$SNP[i],
                             P_disease = summary1[2,4],P_genotype = summary1[3,4],
                             Beta_disease=summary1[2,1],Beta_genotype=summary1[3,1],P_interaction = summary1[4,4])
  
}
cutoff=0.05/nrow(specific_wes)
interaction_wes_sig=interaction_wes[interaction_wes$P_interaction<cutoff,]
write.table(interaction_wes_sig,file = "Interaction_WES.sig.txt",sep = "\t",row.names = F,quote = F)

# permutation
permutation=matrix(nrow = 1000,ncol = 1)
permutation=as.data.frame(permutation)
permutation[1,1]=nrow(interaction_wes_sig)
rownames(permutation)[1]="Real_test"
colnames(permutation)="Number of Significance"

for(n in 1:999){
  
  print(paste("Permutation_number is",n,sep = " "))
  test_set=disease_wes
  test_set$Status=sample(test_set[,1])
  interaction_permutation = foreach(i=1:nrow(specific_wes),.combine = rbind) %do%  {
    probe.sub=as.data.frame(probe_wes[,specific_wes$ProbeName[i]],drop=F,row.names = rownames(probe_wes))
    probe.sub=na.omit(probe.sub)
    snp.sub=specific_wes$SNP[i]
    genotype.sub=genotype_wes[rownames(genotype_wes) %in% rownames(probe.sub),]
    genotype.sub=genotype.sub[,snp.sub,drop=F]
    disease_wes.sub=test_set[rownames(test_set) %in% rownames(probe.sub),,drop=F]
    lm1 = lm(unlist(probe.sub) ~ as.numeric(disease_wes.sub$Status) * as.numeric(genotype.sub[,1]))
    summary1 = summary(lm1)$coef
    return.string = data.frame(ProbeName = specific_wes$ProbeName[i], SNPName = specific_wes$SNP[i],
                               P_disease = summary1[2,4],P_genotype = summary1[3,4],
                               Beta_disease=summary1[2,1],Beta_genotype=summary1[3,1],P_interaction = summary1[4,4])
  }
  count=nrow(interaction_permutation[interaction_permutation$P_interaction<cutoff,])
  permutation[n+1,1]=unlist(count)
  rownames(permutation)[n+1]=paste("Permutation",n,sep = "_")
  
}




#==========================================================
#                  Targeted permutation 1000
#==========================================================

specific_tar=rbind(ibd_specific_2,lld_specific_2)
remove=specific_tar[grep("k__",specific_tar$ProbeName,invert = T),]
remove=remove[grep("unclassified",remove$ProbeName),]
specific_tar=specific_tar[!specific_tar$ProbeName %in% remove$ProbeName,]

ibd_geno_new=ibd_geno_new[,colnames(ibd_geno_new) %in% specific_tar$SNP]
ibd_probe=ibd_probe[,colnames(ibd_probe) %in% specific_tar$ProbeName]
ibd_probe=ibd_probe[rownames(ibd_probe) %in% rownames(ibd_geno_new),]
ibd_probe=ibd_probe[order(rownames(ibd_probe)),]
ibd_geno_new=ibd_geno_new[order(rownames(ibd_geno_new)),]
specific_tar=specific_tar[specific_tar$ProbeName %in% colnames(ibd_probe),]

lld_geno_new=lld_geno_new[,colnames(lld_geno_new) %in% specific_tar$SNP]
lld_probe=lld_probe[,colnames(lld_probe) %in% specific_tar$ProbeName]
lld_probe=lld_probe[rownames(lld_probe) %in% rownames(lld_geno_new),]
lld_probe=lld_probe[order(rownames(lld_probe)),]
lld_geno_new=lld_geno_new[order(rownames(lld_geno_new)),]
specific_tar=specific_tar[specific_tar$ProbeName %in% colnames(lld_probe),]

library(foreach)
library(iterators)
library(RVAideMemoire)

probe_tar=rbind(ibd_probe,lld_probe)
genotype_tar=rbind(ibd_geno_new,lld_geno_new)
disease_tar=genotype_tar[,1:2,drop=F]
colnames(disease_tar)[1]="Status"
disease_tar$Status=as.character(disease_tar$Status)
disease_tar$Status[1:435]=1
disease_tar$Status[436:nrow(disease_tar)]=0
disease_tar[,-1]=NULL

# real one
interaction_tar = foreach(i=1:nrow(specific_tar),.combine = rbind) %do%  {
  probe.sub=as.data.frame(probe_tar[,specific_tar$ProbeName[i]],drop=F,row.names = rownames(probe_tar))
  probe.sub=na.omit(probe.sub)
  snp.sub=specific_tar$SNP[i]
  print(specific_tar[i,])
  genotype.sub=genotype_tar[rownames(genotype_tar) %in% rownames(probe.sub),]
  genotype.sub=genotype.sub[,snp.sub,drop=F]
  disease_tar.sub=disease_tar[rownames(disease_tar) %in% rownames(probe.sub),,drop=F]
  lm1 = lm(unlist(probe.sub) ~ as.numeric(disease_tar.sub$Status) * as.numeric(genotype.sub[,1]))
  summary1 = summary(lm1)$coef
  return.string = data.frame(ProbeName = specific_tar$ProbeName[i], SNPName = specific_tar$SNP[i],
                             P_disease = summary1[2,4],P_genotype = summary1[3,4],
                             Beta_disease=summary1[2,1],Beta_genotype=summary1[3,1],P_interaction = summary1[4,4])
  
}
cutoff=0.05/nrow(specific_tar)
interaction_tar_sig=interaction_tar[interaction_tar$P_interaction<cutoff,]
write.table(interaction_tar_sig,file = "Interaction_tar.sig.txt",sep = "\t",quote = F,row.names = F)

# permutation
permutation=matrix(nrow = 1000,ncol = 1)
permutation=as.data.frame(permutation)
permutation[1,1]=nrow(interaction_tar_sig)
rownames(permutation)[1]="Real_test"
colnames(permutation)="Number of Significance"

for(n in 1:999){
  
  print(paste("Permutation_number is",n,sep = " "))
  test_set=disease_tar
  test_set$Status=sample(test_set[,1])
  interaction_permutation = foreach(i=1:nrow(specific_tar),.combine = rbind) %do%  {
    probe.sub=as.data.frame(probe_tar[,specific_tar$ProbeName[i]],drop=F,row.names = rownames(probe_tar))
    probe.sub=na.omit(probe.sub)
    snp.sub=specific_tar$SNP[i]
    genotype.sub=genotype_tar[rownames(genotype_tar) %in% rownames(probe.sub),]
    genotype.sub=genotype.sub[,snp.sub,drop=F]
    disease_tar.sub=test_set[rownames(test_set) %in% rownames(probe.sub),,drop=F]
    lm1 = lm(unlist(probe.sub) ~ as.numeric(disease_tar.sub$Status) * as.numeric(genotype.sub[,1]))
    summary1 = summary(lm1)$coef
    return.string = data.frame(ProbeName = specific_tar$ProbeName[i], SNPName = specific_tar$SNP[i],
                               P_disease = summary1[2,4],P_genotype = summary1[3,4],
                               Beta_disease=summary1[2,1],Beta_genotype=summary1[3,1],P_interaction = summary1[4,4])
  }
  count=nrow(interaction_permutation[interaction_permutation$P_interaction<cutoff,])
  permutation[n+1,1]=unlist(count)
  rownames(permutation)[n+1]=paste("Permutation",n,sep = "_")
  
}

#=========================================================
#                   calculate CI 
#=========================================================

ibd = foreach(i=1:nrow(specific),.combine = rbind) %do%  {
  
  probe=specific$ProbeName[i]
  snp=specific$SNP[i]
  print(paste(probe,snp,sep = "  "))
  probe.sub=ibd_probe[,probe,drop=F]-100
  genotype.sub=ibd_geno_new[,snp,drop=F]
  probe.sub=na.omit(probe.sub)
  genotype.sub=genotype.sub[rownames(genotype.sub) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  mm=spearman.ci(probe.sub[,1],genotype.sub[,1])
  upper=mm$conf.int[[1]]
  lower=mm$conf.int[[2]]
  nn=cor.test(probe.sub[,1],genotype.sub[,1],method = "spearman")
  Correlationcoeffcient=nn$estimate[[1]]
  pvalue=nn$p.value
  return.string=data.frame(ibd_Probe=probe,ibd_SNP=snp,ibd_PValue=pvalue,
                           ibd_CorrelationCoeffcient=Correlationcoeffcient,
                           ibd_CI_upper=upper,ibd_CI_lower=lower)
  
}

lld = foreach(i=1:nrow(specific),.combine = rbind) %do%  {
  
  probe=specific$ProbeName[i]
  snp=specific$SNP[i]
  print(paste(probe,snp,sep = "  "))
  probe.sub=lld_probe[,probe,drop=F]-100
  genotype.sub=lld_geno_new[,snp,drop=F]
  probe.sub=na.omit(probe.sub)
  genotype.sub=genotype.sub[rownames(genotype.sub) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  mm=spearman.ci(probe.sub[,1],genotype.sub[,1])
  upper=mm$conf.int[[1]]
  lower=mm$conf.int[[2]]
  nn=cor.test(probe.sub[,1],genotype.sub[,1],method = "spearman")
  Correlationcoeffcient=nn$estimate[[1]]
  pvalue=nn$p.value
  return.string=data.frame(lld_Probe=probe,lld_SNP=snp,lld_PValue=pvalue,
                           lld_CorrelationCoeffcient=Correlationcoeffcient,
                           lld_CI_upper=upper,lld_CI_lower=lower)
  
}

ibd$name=paste(ibd$ibd_Probe,ibd$ibd_SNP)
lld$name=paste(lld$lld_Probe,lld$lld_SNP)

table=merge(ibd,lld,by="name",all = F)
table$data=NA
table$data[table$ibd_SNP %in% ibd_specific$SNP]="IBD"
table$data[table$lld_SNP %in% lld_specific$SNP]="LLD"

interaction_1=read.table("Interaction_WES.sig.txt",header = T,sep = "\t",stringsAsFactors = F)
interaction_2=read.table("Interaction_tar.sig.txt",header = T,sep = "\t",stringsAsFactors = F)

interaction=rbind(interaction_1,interaction_2)
table$sig=NA
table$sig[table$ibd_SNP %in% interaction$SNPName]="Interaction Significant"
table$sig[table$lld_SNP %in% interaction$SNPName]="Interaction Significant"
table=na.omit(table)

library(ggplot2)
png("Cohort_specific_mbQTL_correlation_coeffecient.png",width = 1000,height = 1000)
ggplot(data = table,aes(table$lld_CorrelationCoeffcient,table$ibd_CorrelationCoeffcient)) + 
  geom_point(aes(fill=data,shape=sig),shape=21,size=4) + 
  scale_fill_manual(values = c( "red", "#56B4E9"))+
  geom_errorbar(aes(ymin = table$ibd_CI_lower,ymax = table$ibd_CI_upper),width=0, size=0.3,color="grey") + 
  geom_errorbarh(aes(xmin = table$lld_CI_lower,xmax = table$lld_CI_upper),size=0.3,color="grey") +
  theme_bw()+
  theme_bw()+theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(), 
                   panel.background = element_blank())+guides(fill=guide_legend(title=NULL))+
  geom_vline(xintercept = 0,color="grey")+geom_hline(yintercept = 0,color="grey")+
  #theme(legend.position="rightbottom")+
  xlab("LLD cohort (Correlation Coeffcient with 95% CI)")+
  ylab("IBD cohort (Correlation Coeffcient with 95% CI)")+
  #guides(fill=F)+
  theme(axis.title = element_text(size = 30, vjust = 0.5, hjust = 0.5))+
  theme(axis.text = element_text(size = 30, vjust = 0.5, hjust = 0.5))
dev.off()

