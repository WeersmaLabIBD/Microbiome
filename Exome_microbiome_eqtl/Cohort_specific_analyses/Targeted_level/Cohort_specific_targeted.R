#####################################################################################################################
#####################################  cohort specific  #############################################################
#####################################################################################################################
# attention: genotype does not include all, it is extracted from genotyper_hamonizer outcome

### IBD_specific
ibd_genotype=read.table("IBD_dosage_target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
lld_genotype=read.table("LLD_dosage_target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genotype=merge(ibd_genotype,lld_genotype,all = F,by="SNP")
rownames(genotype)=genotype[,1]
genotype=genotype[,-1]
genotype=as.data.frame(t(genotype))
genotype=genotype[order(rownames(genotype)),]

ibd_probe=read.table("IBD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(ibd_probe)[1]="Probe"
lld_probe=read.table("LLD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(lld_probe)[1]="Probe"
probe=merge(ibd_probe,lld_probe,all = F,by="Probe")
ibd_specific=read.table("IBD_specific.txt",sep = "\t",header = T)
genotype=genotype[,colnames(genotype) %in% ibd_specific$SNP]
probe=probe[probe$Probe %in% ibd_specific$ProbeName,]
rownames(probe)=probe[,1]
probe=probe[,-1]
probe=as.data.frame(t(probe))
probe=probe[rownames(probe) %in% rownames(genotype),]
probe=probe[order(rownames(probe)),]
probe[probe==0]=NA
ibd_specific=ibd_specific[ibd_specific$ProbeName %in% colnames(probe),]

disease_status=genotype[,1:2,drop=F]
colnames(disease_status)[1]="Status"
disease_status$Status=as.character(disease_status$Status)
disease_status$Status[1:500]=1
disease_status$Status[501:nrow(disease_status)]=0
disease_status=disease_status[order(rownames(disease_status)),]

library(foreach)
library(iterators)
library(parallel)
library(doParallel)

ibd = foreach(i=1:nrow(ibd_specific),.combine = rbind) %do%  {
  probe.sub=as.data.frame(probe[,ibd_specific$ProbeName[i]],drop=F,row.names = rownames(probe))
  probe.sub=na.omit(probe.sub)
  genotype.sub=genotype[rownames(genotype) %in% rownames(probe.sub),]
  disease_status.sub=disease_status[rownames(disease_status) %in% rownames(probe.sub),]
  lm1 = lm(unlist(probe.sub) ~ as.numeric(disease_status.sub$Status) * as.numeric(genotype.sub[,ibd_specific$SNP[i]]))
  summary1 = summary(lm1)$coef
  return.string = data.frame(ProbeName = ibd_specific$ProbeName[i], SNPName = ibd_specific$SNP[i],P_disease = summary1[2,4],P_genotype = summary1[3,4],Beta_disease=summary1[2,1],Beta_genotype=summary1[3,1],P_interaction = summary1[4,4])
}

### LLD_specific
ibd_genotype=read.table("IBD_dosage_target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
lld_genotype=read.table("LLD_dosage_target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genotype=merge(ibd_genotype,lld_genotype,all = F,by="SNP")
rownames(genotype)=genotype[,1]
genotype=genotype[,-1]
genotype=as.data.frame(t(genotype))
genotype=genotype[order(rownames(genotype)),]

ibd_probe=read.table("IBD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(ibd_probe)[1]="Probe"
lld_probe=read.table("LLD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(lld_probe)[1]="Probe"
probe=merge(ibd_probe,lld_probe,all = F,by="Probe")
lld_specific=read.table("LLD_specific.txt",sep = "\t",header = T)
genotype=genotype[,colnames(genotype) %in% lld_specific$SNP]
probe=probe[probe$Probe %in% lld_specific$ProbeName,]
rownames(probe)=probe[,1]
probe=probe[,-1]
probe=as.data.frame(t(probe))
probe=probe[rownames(probe) %in% rownames(genotype),]
probe=probe[order(rownames(probe)),]
probe[probe==0]=NA
lld_specific=lld_specific[lld_specific$ProbeName %in% colnames(probe),]

disease_status=genotype[,1:2,drop=F]
colnames(disease_status)[1]="Status"
disease_status$Status=as.character(disease_status$Status)
disease_status$Status[1:500]=1
disease_status$Status[501:nrow(disease_status)]=0
disease_status=disease_status[order(rownames(disease_status)),]

library(foreach)
library(iterators)
library(parallel)
library(doParallel)

lld = foreach(i=1:nrow(lld_specific),.combine = rbind) %do%  {
  probe.sub=as.data.frame(probe[,lld_specific$ProbeName[i]],drop=F,row.names = rownames(probe))
  probe.sub=na.omit(probe.sub)
  genotype.sub=genotype[rownames(genotype) %in% rownames(probe.sub),]
  disease_status.sub=disease_status[rownames(disease_status) %in% rownames(probe.sub),]
  lm1 = lm(unlist(probe.sub) ~ as.numeric(disease_status.sub$Status) * as.numeric(genotype.sub[,lld_specific$SNP[i]]))
  summary1 = summary(lm1)$coef
  return.string = data.frame(ProbeName = lld_specific$ProbeName[i], SNPName = lld_specific$SNP[i],P_disease = summary1[2,4],P_genotype = summary1[3,4],Beta_disease=summary1[2,1],Beta_genotype=summary1[3,1],P_interaction = summary1[4,4])
}

### merge
ibd$Group="IBDcohort_significant"
lld$Group="LLDcohort_significant"
all=rbind(ibd,lld)
write.table(all,file="targeted_cohort_specific_mbQTL_interactions.txt",row.names = FALSE,quote = FALSE,sep = "\t")

significant_P_interaction=read.table("targeted_cohort_specific_mbQTL_interactions.txt",header = T,sep = "\t",stringsAsFactors = F)
cutoff=0.05/nrow(significant_P_interaction)
significant_P_interaction=significant_P_interaction[significant_P_interaction$P_interaction<cutoff,]
ibd_genotype=read.table("IBD_dosage_target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
lld_genotype=read.table("LLD_dosage_target.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genotype=merge(ibd_genotype,lld_genotype,all = F,by="SNP")
rownames(genotype)=genotype[,1]
genotype=genotype[,-1]
genotype=as.data.frame(t(genotype))
genotype=genotype[order(rownames(genotype)),]

ibd_allele=read.table("IBD_genotype_target.txt",header = T,stringsAsFactors = F,check.names = F,sep="\t")
lld_allele=read.table("LLD_genotype_target.txt",header = T,stringsAsFactors = F,check.names = F,sep="\t")
ibd_allele[ibd_allele=="0/0"]=NA
lld_allele[lld_allele=="0/0"]=NA

ibd_probe=read.table("IBD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(ibd_probe)[1]="Probe"
lld_probe=read.table("LLD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(lld_probe)[1]="Probe"
probe=merge(ibd_probe,lld_probe,all = F,by="Probe")

# standered plot
library(ggplot2)
pdf("Targeted_cohort_specific.pdf",width=10, height=7)
for(i in 1:nrow(significant_P_interaction)){
  
  geno.sub=genotype[,significant_P_interaction$SNPName[i],drop=F]
  geno.sub$Group=c(rep("IBD",500),rep("LLD",920))
  ibd.allele.sub=ibd_allele[ibd_allele$SNP==colnames(geno.sub)[1],,drop=F]
  ibd.allele.sub=as.data.frame(t(ibd.allele.sub))
  colnames(ibd.allele.sub)=ibd.allele.sub[1,1]
  ibd.allele.sub=ibd.allele.sub[-1,,drop=F]
  lld.allele.sub=lld_allele[lld_allele$SNP==colnames(geno.sub)[1],,drop=F]
  lld.allele.sub=as.data.frame(t(lld.allele.sub))
  colnames(lld.allele.sub)=lld.allele.sub[1,1]
  lld.allele.sub=lld.allele.sub[-1,,drop=F]
  allele.sub=rbind(ibd.allele.sub,lld.allele.sub)
  allele.sub=allele.sub[rownames(allele.sub) %in% rownames(geno.sub),,drop=F]
  allele.sub=allele.sub[order(rownames(allele.sub)),,drop=F]
  geno.sub=geno.sub[order(rownames(geno.sub)),,drop=F]
  geno.sub=merge(allele.sub,geno.sub,by="row.names",all=T)
  rownames(geno.sub)=geno.sub$Row.names
  geno.sub=geno.sub[,-1]
  
  pro.sub=probe[probe$Probe==significant_P_interaction$ProbeName[i],,drop=F]
  pro.sub=as.data.frame(t(pro.sub),stringsAsFactors = F)
  colnames(pro.sub)=pro.sub[1,1]
  pro.sub=pro.sub[-1,,drop=F]
  pro.sub[pro.sub==0]=NA
  pro.sub=na.omit(pro.sub)
  table=merge(pro.sub,geno.sub,all = F,by="row.names")
  table=na.omit(table)
  colnames(table)=c("sample","probe","genotype","dosage","group")
  mm=as.data.frame(table(unlist(table$genotype)))
  mm[mm==0]=NA
  mm=na.omit(mm)
  nn=rbind(mm[mm$Var1=="A/A",],mm[mm$Var1=="T/T",],mm[mm$Var1=="G/G",],mm[mm$Var1=="C/C",])
  nn=nn[order(nn$Var1),]
  table$genotype=factor(table$genotype,levels = c(as.character(nn$Var1[1]),as.character(mm$Var1[!mm$Var1 %in% nn$Var1]),as.character(nn$Var1[2])))
  table$probe=scale(as.numeric(table$probe))
  table=na.omit(table)
  table$face="IBD+LLD"
  
  p=ggplot(table, aes(genotype,probe)) + 
    geom_boxplot()+
    facet_grid(~group)+
    theme_bw()+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
    stat_boxplot(geom = "errorbar",linetype="dotted")+
    geom_boxplot(outlier.size = 0,size=0,outlier.colour="white")+
    geom_dotplot(binaxis = "y", stackdir = "center",dotsize = 2,alpha=0.5,binwidth=0.03,aes(fill=group))+
    scale_fill_manual(values=c("red","blue"))+
    guides(fill=FALSE)+
    theme(plot.title = element_text(size = 25,vjust = 0.5, hjust = 0.5))+
    ylab(significant_P_interaction$ProbeName[i])+xlab("")+ylim(-2,3)+
    theme(axis.title = element_text(size = 10, vjust = 0.5, hjust = 0.5))+
    theme(axis.text = element_text(size = 10, vjust = 0.5, hjust = 0.5))+
    theme(strip.text = element_text(size=25, face="bold"))
 print(p)
}
dev.off()

# testing plot
#library(ggplot2)
#pdf("hehe.pdf",width=10, height=7)
#for(i in 1:nrow(significant_P_interaction)){
  
#  geno.sub=genotype[,significant_P_interaction$SNPName[i],drop=F]
#  geno.sub$Group=c(rep("IBD",500),rep("LLD",920))
#  ibd.allele.sub=ibd_allele[ibd_allele$SNP==colnames(geno.sub)[1],,drop=F]
#  ibd.allele.sub=as.data.frame(t(ibd.allele.sub))
#  colnames(ibd.allele.sub)=ibd.allele.sub[1,1]
#  ibd.allele.sub=ibd.allele.sub[-1,,drop=F]
#  lld.allele.sub=lld_allele[lld_allele$SNP==colnames(geno.sub)[1],,drop=F]
#  lld.allele.sub=as.data.frame(t(lld.allele.sub))
#  colnames(lld.allele.sub)=lld.allele.sub[1,1]
#  lld.allele.sub=lld.allele.sub[-1,,drop=F]
#  allele.sub=rbind(ibd.allele.sub,lld.allele.sub)
#  allele.sub=allele.sub[rownames(allele.sub) %in% rownames(geno.sub),,drop=F]
#  allele.sub=allele.sub[order(rownames(allele.sub)),,drop=F]
#  geno.sub=geno.sub[order(rownames(geno.sub)),,drop=F]
#  geno.sub=merge(allele.sub,geno.sub,by="row.names",all=T)
#  rownames(geno.sub)=geno.sub$Row.names
#  geno.sub=geno.sub[,-1]
  
#  pro.sub=probe[probe$Probe==significant_P_interaction$ProbeName[i],,drop=F]
#  pro.sub=as.data.frame(t(pro.sub),stringsAsFactors = F)
#  colnames(pro.sub)=pro.sub[1,1]
#  pro.sub=pro.sub[-1,,drop=F]
#  pro.sub[pro.sub==0]=NA
#  pro.sub=na.omit(pro.sub)
#  table=merge(pro.sub,geno.sub,all = F,by="row.names")
#  table=na.omit(table)
#  colnames(table)=c("sample","probe","genotype","dosage","group")
#  mm=as.data.frame(table(unlist(table$genotype)))
#  mm[mm==0]=NA
#  mm=na.omit(mm)
#  nn=rbind(mm[mm$Var1=="A/A",],mm[mm$Var1=="T/T",],mm[mm$Var1=="G/G",],mm[mm$Var1=="C/C",])
#  nn=nn[order(nn$Var1),]
#  table$genotype=factor(table$genotype,levels = c(as.character(nn$Var1[1]),as.character(mm$Var1[!mm$Var1 %in% nn$Var1]),as.character(nn$Var1[2])))
#  table$probe=scale(as.numeric(table$probe))
#  table=na.omit(table)
#  table$face="IBD+LLD"
  
#  p1=ggplot(table, aes(genotype,probe)) + 
#    geom_boxplot()+
#    facet_grid(~group)+
#    theme_bw()+
#    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
#    stat_boxplot(geom = "errorbar",linetype="dotted")+
#    geom_boxplot(outlier.size = 0,size=0,outlier.colour="white")+
#    geom_dotplot(binaxis = "y", stackdir = "center",dotsize = 2,alpha=0.5,binwidth=0.03,aes(fill=group))+
#    scale_fill_manual(values=c("red","blue"))+
#    guides(fill=FALSE)+
#    theme(plot.title = element_text(size = 25,vjust = 0.5, hjust = 0.5))+
#    ylab(significant_P_interaction$ProbeName[i])+xlab("")+ylim(-2,3)+
#    theme(axis.title = element_text(size = 10, vjust = 0.5, hjust = 0.5))+
#    theme(axis.text = element_text(size = 10, vjust = 0.5, hjust = 0.5))+
#    theme(strip.text = element_text(size=25, face="bold"))
#  lm1=summary(lm(table[table$group=="IBD",]$probe~table[table$group=="IBD",]$dosage))$coef
#  lm2=summary(lm(table[table$group=="LLD",]$probe~table[table$group=="LLD",]$dosage))$coef
#  p2=ggplot(table, aes(genotype,probe)) + 
#    theme_bw()+
#    facet_grid(~face)+
#    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
#    geom_dotplot(binaxis = "y", stackdir = "center",dotsize = 2,alpha=0.5,binwidth=0.03,aes(fill=group))+
#    scale_fill_manual(values=c("red","blue"))+
#    guides(fill=FALSE)+
#    theme(plot.title = element_text(size = 25,vjust = 0.5, hjust = 0.5))+
#    ylab("")+xlab("")+ylim(-2,3)+
#    theme(axis.title = element_text(size = 10, vjust = 0.5, hjust = 0.5))+
#    theme(axis.text = element_text(size = 10, vjust = 0.5, hjust = 0.5))+
#    theme(strip.text = element_text(size=25, face="bold"))+
#    geom_abline(intercept=lm1[1,1],slope=lm1[2,1],color="red",size=2,linetype="dashed")+
#    geom_abline(intercept=lm2[1,1],slope=lm2[2,1],color="blue",size=2,linetype="dashed")
#  print(multiplot(p1,p2,cols = 2))
#}
#dev.off()

# double check
#for(i in 1:nrow(significant_P_interaction)){
#  geno.sub=ibd_genotype[ibd_genotype$SNP==significant_P_interaction$SNPName[i],,drop=F]
#  geno.sub=as.data.frame(t(geno.sub),stringsAsFactors = F)
#  colnames(geno.sub)=geno.sub[1,1]
#  geno.sub=geno.sub[-1,,drop=F]
#  pro.sub=ibd_probe[ibd_probe$Probe==significant_P_interaction$ProbeName[i],,drop=F]
#  pro.sub=as.data.frame(t(pro.sub),stringsAsFactors = F)
#  colnames(pro.sub)=pro.sub[1,1]
#  pro.sub=pro.sub[-1,,drop=F]
#  pro.sub[pro.sub==0]=NA
#  pro.sub=na.omit(pro.sub)
#  table=merge(pro.sub,geno.sub,all = F,by="row.names")
#  colnames(table)=c("sampe","probe","genotype")
#  table$probe=scale(as.numeric(table$probe))
#  ggplot(table, aes(genotype,probe)) + 
#    geom_boxplot()+
#    theme_bw()+
#    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
#    stat_boxplot(geom = "errorbar",linetype="dotted")+
#    geom_boxplot(outlier.size = 0,size=0,outlier.colour="white")+
#    geom_dotplot(binaxis = "y", stackdir = "center",dotsize = 2,alpha=0.3,binwidth=0.03)
#}
