# mbQTL mapping for low-frequency and rare mutations

-------------------------------------------

 
 To perform mbQTL mapping for rare mutations, you need the following files, including:
 
 1) IBD/LLD_info.txt: each SNP with gene name, for example,
 
 ```
 SNP1     GENE1
 SNP2     GENE1
 SNP3     GENE1
 SNP4     GENE2
 SNP5     GENE3
 ```
 
 2) IBD/LLD_numeric.txt: corrected bacterial data
 
 3) IBD/LLD_low_PTV.table.dosage.txt: dosage of low frequency and rare PTV variants
 
# Step.1 Read files, load R packages SKAT and MetaSKAT

```
library(SKAT)
library(MetaSKAT)

ibd_info=read.table("IBD_info.txt",header = F,sep = "\t",stringsAsFactors = F,check.names = F)
lld_info=read.table("LLD_info.txt",header = F,sep = "\t",stringsAsFactors = F,check.names = F)
colnames(ibd_info)=c("Name","gene")
colnames(lld_info)=c("Name","gene")
ibd_probe=read.table("IBD_numeric.txt",sep = "\t",header = T, stringsAsFactors= F,check.names = F,row.names = 1)
lld_probe=read.table("LLD_numeric.txt",sep = "\t",header = T,stringsAsFactors= F,check.names = F,row.names = 1)
ibd_probe=as.data.frame(t(ibd_probe))
lld_probe=as.data.frame(t(lld_probe))

ibd_dosage=read.table("IBD_low_PTV.table.dosage.txt",header = T,check.names = F,stringsAsFactors = F,sep = "\t",row.names = 1)
lld_dosage=read.table("LLD_low_PTV.talbe.dosage.txt",header = T,check.names = F,stringsAsFactors = F,sep = "\t",row.names = 1)

ibd_dosage=as.data.frame(t(ibd_dosage))
lld_dosage=as.data.frame(t(lld_dosage))
lld_pheno=as.data.frame(lld_dosage[,1])
ibd_pheno=as.data.frame(ibd_dosage[,1])
lld_pheno[,1]=1
ibd_pheno[,1]=1
ibd_probe=ibd_probe[rownames(ibd_probe) %in% rownames(ibd_dosage),]
lld_probe=lld_probe[rownames(lld_probe) %in% rownames(lld_dosage),]
ibd_dosage[ibd_dosage==-1]=NA
lld_dosage[lld_dosage==-1]=NA

all_info=rbind(ibd_info,lld_info)
all_info=all_info[!duplicated(all_info),]
all_info=all_info[order(all_info$gene),]
ibd_info=ibd_info[order(ibd_info$gene),]
lld_info=lld_info[order(lld_info$gene),]
has_both=colnames(lld_probe)[colnames(lld_probe) %in% colnames(ibd_probe)]
genes=all_info[!duplicated(all_info$gene),]

```

# Step.2 Single variant test

```
ibd_dosage[is.na(ibd_dosage)]="QQ"
for(i in 1:ncol(ibd_dosage)){
  
  ibd_dosage.sub=na.omit(ibd_dosage[,i])
  ibd_dosage.sub=as.numeric(ibd_dosage.sub[ibd_dosage.sub!="QQ"])
  sum_dosage=sum(ibd_dosage.sub)
  if(sum_dosage<800){
    
    ibd_dosage[ibd_dosage[,i]==2,i]="M"
    ibd_dosage[ibd_dosage[,i]==0,i]="m"
    ibd_dosage[ibd_dosage[,i]=="M",i]=0
    ibd_dosage[ibd_dosage[,i]=="m",i]=2
    
  }
}
ibd_dosage[ibd_dosage=="QQ"]=NA

lld_num=as.numeric(nrow(lld_dosage))
lld_dosage[is.na(lld_dosage)]="QQ"
for(i in 1:ncol(lld_dosage)){
  
  lld_dosage.sub=na.omit(lld_dosage[,i])
  lld_dosage.sub=as.numeric(lld_dosage.sub[lld_dosage.sub!="QQ"])
  sum_dosage=sum(lld_dosage.sub)
  if(sum_dosage<800){
    
    lld_dosage[lld_dosage[,i]==2,i]="M"
    lld_dosage[lld_dosage[,i]==0,i]="m"
    lld_dosage[lld_dosage[,i]=="M",i]=0
    lld_dosage[lld_dosage[,i]=="m",i]=2
    
  }
}
lld_dosage[lld_dosage=="QQ"]=NA

singleton=merge(ibd_info,lld_info,by="Name",all=F)
library(foreach)
single_ibd = foreach(i=1:nrow(singleton),.combine = rbind) %do%  {
  
  variant=singleton$Name[i]
  genotype.sub=ibd_dosage[,colnames(ibd_dosage)==variant,drop=F]
  aa=length(genotype.sub[genotype.sub==0,])
  bb=length(genotype.sub[genotype.sub==2,])
  homozygous=min(aa,bb)
  if(homozygous>=2){
    
    return.string=data.frame(SNP=variant,Gene=singleton$gene.x[i])
    
  }
}

single_lld = foreach(i=1:nrow(singleton),.combine = rbind) %do%  {
  
  variant=singleton$Name[i]
  genotype.sub=lld_dosage[,colnames(lld_dosage)==variant,drop=F]
  aa=length(genotype.sub[genotype.sub==0,])
  bb=length(genotype.sub[genotype.sub==2,])
  homozygous=min(aa,bb)
  if(homozygous>=2){
    
    return.string=data.frame(SNP=variant,Gene=singleton$gene.x[i])
    
  }
}  

singleton=merge(single_ibd,single_lld,by="SNP",all = F)
write.table(singleton,"Single_variants.txt",sep = "\t",quote = F,row.names = F)

singleton$SNP=as.character(singleton$SNP)
singleton$Gene.x=as.character(singleton$Gene.x)
ibd_spearman = foreach(i=1:nrow(singleton),.combine = rbind) %do%  {
  ibd.geno.sub=ibd_dosage[,colnames(ibd_dosage)==singleton$SNP[i],drop=F]
  print(paste(singleton$SNP[i],"Start",sep = "-----"))
  tmp=foreach(n=has_both,.combine = rbind)  %do% {
    
    ibd.pheno.sub=ibd_probe[,colnames(ibd_probe)==n,drop=F]
    ibd.pheno.sub[ibd.pheno.sub==0]=NA
    ibd.pheno.sub=na.omit(ibd.pheno.sub)
    ibd.dose.sub=ibd.geno.sub[rownames(ibd.geno.sub) %in% rownames(ibd.pheno.sub),,drop=F]
    ibd.dose.sub=ibd.dose.sub[order(rownames(ibd.dose.sub)),,drop=F]
    ibd.pheno.sub=ibd.pheno.sub[order(rownames(ibd.pheno.sub)),,drop=F]
    number=nrow(ibd.dose.sub)
    mm=cor.test(as.numeric(ibd.pheno.sub[,1]),as.numeric(ibd.dose.sub[,1]),method = "spearman")
    return.string=data.frame(PValue=mm$p.value,CorrelationCoefficient=mm$estimate,SNP=singleton$SNP[i],Gene=singleton$Gene.x[i],Probe=n,
                             Number=number)
  }
  
  return.string=data.frame(PValue=tmp$PValue,CorrelationCoefficient=tmp$CorrelationCoefficient,SNP=tmp$SNP,Gene=tmp$Gene,Probe=tmp$Probe,
                           Number=tmp$Number)
  
}

lld_spearman = foreach(i=1:nrow(singleton),.combine = rbind) %do%  {
  lld.geno.sub=lld_dosage[,colnames(lld_dosage)==singleton$SNP[i],drop=F]
  print(paste(singleton$SNP[i],"Start",sep = "-----"))
  tmp=foreach(n=has_both,.combine = rbind)  %do% {
    
    lld.pheno.sub=lld_probe[,colnames(lld_probe)==n,drop=F]
    lld.pheno.sub[lld.pheno.sub==0]=NA
    lld.pheno.sub=na.omit(lld.pheno.sub)
    lld.dose.sub=lld.geno.sub[rownames(lld.geno.sub) %in% rownames(lld.pheno.sub),,drop=F]
    lld.dose.sub=lld.dose.sub[order(rownames(lld.dose.sub)),,drop=F]
    lld.pheno.sub=lld.pheno.sub[order(rownames(lld.pheno.sub)),,drop=F]
    number=nrow(ibd.dose.sub)
    mm=cor.test(as.numeric(lld.pheno.sub[,1]),as.numeric(lld.dose.sub[,1]),method = "spearman")
    return.string=data.frame(PValue=mm$p.value,CorrelationCoefficient=mm$estimate,SNP=singleton$SNP[i],Gene=singleton$Gene.x[i],Probe=n,
                             Number=number)
  }
  
  return.string=data.frame(PValue=tmp$PValue,CorrelationCoefficient=tmp$CorrelationCoefficient,SNP=tmp$SNP,Gene=tmp$Gene,Probe=tmp$Probe,
                           Number=tmp$Number)
  
}

library(metap)
ibd_spearman$name=paste(ibd_spearman$Gene,ibd_spearman$Probe)
lld_spearman$name=paste(lld_spearman$Gene,lld_spearman$Probe)
ibd_spearman=ibd_spearman[order(ibd_spearman$name),]
lld_spearman=lld_spearman[order(lld_spearman$name),]
ibd_spearman[is.na(ibd_spearman)]=0.99999
lld_spearman[is.na(lld_spearman)]=0.99999
ibd_spearman$PValue[ibd_spearman$PValue==1]=0.99999
lld_spearman$PValue[lld_spearman$PValue==1]=0.99999
meta_result=foreach(i=1:nrow(ibd_spearman),.combine = rbind) %do% {
  print(i)
  p_vector=c(ibd_spearman$PValue[i],lld_spearman$PValue[i])
  w_vector=c(ibd_spearman$Number[i],lld_spearman$Number[i])
  meta_z=sumz(p_vector,weights = w_vector)
  return.string=data.frame(Meta_Z=meta_z$z,meta_P=meta_z$p,IBD_Coefficient=ibd_spearman$CorrelationCoefficient[i],LLD_Coefficient=lld_spearman$CorrelationCoefficient[i],
                           IBD_P=ibd_spearman$PValue[i],LLD_P=lld_spearman$PValue[i],SNP=ibd_spearman$SNP[i],Gene=ibd_spearman$Gene[i],Probe=ibd_spearman$Probe[i])
}

meta_result=meta_result[order(meta_result$meta_P),]
cutoff=0.05/nrow(singleton)
meta_sig=meta_result[meta_result$meta_P<cutoff,]
write.table(meta_sig,file = "Single_variant_sig.txt",row.names = F,quote = F,sep = "\t")
```

# Step.3 Gene-based burden test

Burden test in each cohort
```
ibd_non_singleton=ibd_info[duplicated(ibd_info$gene),]
lld_non_singleton=lld_info[duplicated(lld_info$gene),]
ibd_gene=unique(ibd_non_singleton$gene)
lld_gene=unique(lld_non_singleton$gene)
non_singleton=union(ibd_gene,lld_gene)
all_info=all_info[all_info$gene %in% non_singleton,]
genes=all_info[!duplicated(all_info$gene),]
library(foreach)
burden_ibd = foreach(i=1:length(has_both),.combine = rbind) %do%  {
  
  probe=has_both[i]
  print(paste("IBD",probe,sep = "-----"))
  ibd_probe.sub=ibd_probe[,probe,drop=F]
  ibd_probe.sub[ibd_probe.sub==0]=NA
  ibd_probe.sub=na.omit(ibd_probe.sub)
  
  ibd_dosage.sub=ibd_dosage[rownames(ibd_dosage) %in% rownames(ibd_probe.sub),]
  ibd_probe.sub=ibd_probe.sub[order(rownames(ibd_probe.sub)),,drop=F]
  ibd_dosage.sub=ibd_dosage.sub[order(rownames(ibd_dosage.sub)),,drop=F]
  
  sub.string=foreach(n=genes$gene,.combine = rbind) %do%  {
    
    z_sample=all_info[all_info$gene==n,]
    z_sample=ibd_dosage.sub[,colnames(ibd_dosage.sub) %in% z_sample$Name,drop=F]
    
    if(ncol(z_sample)!=0){
      
      obj=SKAT_Null_Model(as.matrix(ibd_probe.sub) ~ 1, out_type = "C")
      outcome=SKAT(as.matrix(z_sample),obj,r.corr=1)
      return.string=data.frame(PValue=outcome$p.value,NMarker=outcome$param$n.marker,VMarker=outcome$param$n.marker.test,Genes=n)
      
    }else{
      
      return.string=data.frame(PValue=NA,NMarker=NA,VMarker=NA,Genes=n)
    }
  }
  
  return.string = data.frame(Pvalue=sub.string$PValue,NMarker=sub.string$NMarker,VMarker=sub.string$VMarker,Genes=sub.string$Genes,Probe=probe)
  
}
write.table(burden_ibd,file = "IBD_burden_raw.txt",row.names = F,sep = "\t",quote = F)

library(foreach)
burden_lld = foreach(i=1:length(has_both),.combine = rbind) %do%  {
  
  probe=has_both[i]
  print(paste("LLD",probe,sep = "-----"))
  lld_probe.sub=lld_probe[,probe,drop=F]
  lld_probe.sub[lld_probe.sub==0]=NA
  lld_probe.sub=na.omit(lld_probe.sub)
  
  lld_dosage.sub=lld_dosage[rownames(lld_dosage) %in% rownames(lld_probe.sub),]
  lld_probe.sub=lld_probe.sub[order(rownames(lld_probe.sub)),,drop=F]
  lld_dosage.sub=lld_dosage.sub[order(rownames(lld_dosage.sub)),,drop=F]
  
  sub.string=foreach(n=genes$gene,.combine = rbind) %do%  {
    
    z_sample=all_info[all_info$gene==n,]
    z_sample=lld_dosage.sub[,colnames(lld_dosage.sub) %in% z_sample$Name,drop=F]
    
    if(ncol(z_sample)!=0){
      
      obj=SKAT_Null_Model(as.matrix(lld_probe.sub) ~ 1, out_type = "C")
      outcome=SKAT(as.matrix(z_sample),obj,r.corr=1)
      return.string=data.frame(PValue=outcome$p.value,NMarker=outcome$param$n.marker,VMarker=outcome$param$n.marker.test,Genes=n)
      
    }else{
      
      return.string=data.frame(PValue=NA,NMarker=NA,VMarker=NA,Genes=n)
    }
  }
  
  return.string = data.frame(Pvalue=sub.string$PValue,NMarker=sub.string$NMarker,VMarker=sub.string$VMarker,Genes=sub.string$Genes,Probe=probe)
  
}
write.table(burden_lld,file = "LLD_burden_raw.txt",row.names = F,sep = "\t",quote = F)

burden_ibd=read.table("IBD_burden_raw.txt",sep = "\t",header = T,stringsAsFactors = F)
burden_lld=read.table("LLD_burden_raw.txt",sep = "\t",header = T,stringsAsFactors = F)
```

Meta_test_1, discovery in ibd, replication in lld
```
ibd_sig=burden_ibd[burden_ibd$Pvalue<0.001,,drop=F]
lld_sig=burden_lld[burden_lld$Pvalue<0.05,,drop=F]
ibd_sig=na.omit(ibd_sig)
lld_sig=na.omit(lld_sig)

ibd_sig$name=paste(ibd_sig$Genes,ibd_sig$Probe,sep = "QQQ")
lld_sig$name=paste(lld_sig$Genes,lld_sig$Probe,sep = "QQQ")
meta_feature=intersect(ibd_sig$name,lld_sig$name)
nn=strsplit(as.character(meta_feature),"QQQ")

library(SKAT)
library(MetaSKAT)
library(foreach)
burden_meta_test_1 = foreach(i=1:length(nn),.combine = rbind) %do%  {
  
  probe=nn[[i]][2]
  variant=nn[[i]][1]
  print(paste("Meta",paste(variant,probe,sep = "--versus--"),sep = "-----"))
  ibd_probe.sub=ibd_probe[,probe,drop=F]
  ibd_probe.sub[ibd_probe.sub==0]=NA
  ibd_probe.sub=na.omit(ibd_probe.sub)
  lld_probe.sub=lld_probe[,probe,drop=F]
  lld_probe.sub[lld_probe.sub==0]=NA
  lld_probe.sub=na.omit(lld_probe.sub)
  
  ibd_dosage.sub=ibd_dosage[rownames(ibd_dosage) %in% rownames(ibd_probe.sub),]
  lld_dosage.sub=lld_dosage[rownames(lld_dosage) %in% rownames(lld_probe.sub),]
  
  ibd_probe.sub=ibd_probe.sub[order(rownames(ibd_probe.sub)),]
  ibd_dosage.sub=ibd_dosage.sub[order(rownames(ibd_dosage.sub)),,drop=F]
  lld_probe.sub=lld_probe.sub[order(rownames(lld_probe.sub)),]
  lld_dosage.sub=lld_dosage.sub[order(rownames(lld_dosage.sub)),,drop=F]
  y.list=list(as.vector(ibd_probe.sub),as.vector(lld_probe.sub))
  
  ibd_pheno.sub=ibd_pheno[1:nrow(ibd_dosage.sub),,drop=F]
  lld_pheno.sub=lld_pheno[1:nrow(lld_dosage.sub),,drop=F]
  colnames(ibd_pheno.sub)="pheno"
  colnames(lld_pheno.sub)="pheno"
  x.list=list(as.matrix(ibd_pheno.sub),as.matrix(lld_pheno.sub))
  
  z_sample=all_info[all_info$gene==variant,]
  lld_z_sample=lld_dosage.sub[,colnames(lld_dosage.sub) %in% z_sample$Name,drop=F]
  ibd_z_sample=ibd_dosage.sub[,colnames(ibd_dosage.sub) %in% z_sample$Name,drop=F]
  lld_z_sample=as.data.frame(t(lld_z_sample))
  ibd_z_sample=as.data.frame(t(ibd_z_sample))
  z_sample=merge(ibd_z_sample,lld_z_sample,all = T,by="row.names")
  rownames(z_sample)=z_sample$Row.names
  z_sample=z_sample[,-1]
  z_sample=as.data.frame(t(z_sample))
  z_sample[is.na(z_sample)]=2
  z_sample <- data.frame(lapply(z_sample,as.numeric),check.names = F,row.names = rownames(z_sample))
  
  obj<-Meta_Null_Model(y.list,x.list, n.cohort=2, out_type="C")
  meta.sub=MetaSKAT_wZ(as.matrix(z_sample), obj,r.corr=1)
  return.string=data.frame(PValue=meta.sub$p.value,Genes=variant,Feature=probe)
  
}
write.table(burden_meta_test_1,file = "Burden_meta_raw_test_1.txt",quote = F,row.names = F,sep = "\t")
cutoff=0.05/980
meta_sig_test_1=burden_meta_test_1[burden_meta_test_1$PValue<cutoff,]
meta_sig_test_1=na.omit(meta_sig_test_1)
write.table(meta_sig_test_1,file = "Meta_significant_test_1.txt",quote = F,row.names = F,sep = "\t")
```

Meta_test_2, discovery in lld, replication in ibd
```
ibd_sig=burden_ibd[burden_ibd$Pvalue<0.05,,drop=F]
lld_sig=burden_lld[burden_lld$Pvalue<0.001,,drop=F]
ibd_sig=na.omit(ibd_sig)
lld_sig=na.omit(lld_sig)

ibd_sig$name=paste(ibd_sig$Genes,ibd_sig$Probe,sep = "QQQ")
lld_sig$name=paste(lld_sig$Genes,lld_sig$Probe,sep = "QQQ")
meta_feature=intersect(ibd_sig$name,lld_sig$name)
nn=strsplit(as.character(meta_feature),"QQQ")

library(SKAT)
library(MetaSKAT)
library(foreach)
burden_meta_test_2 = foreach(i=1:length(nn),.combine = rbind) %do%  {
  
  probe=nn[[i]][2]
  variant=nn[[i]][1]
  print(paste("Meta",paste(variant,probe,sep = "--versus--"),sep = "-----"))
  ibd_probe.sub=ibd_probe[,probe,drop=F]
  ibd_probe.sub[ibd_probe.sub==0]=NA
  ibd_probe.sub=na.omit(ibd_probe.sub)
  lld_probe.sub=lld_probe[,probe,drop=F]
  lld_probe.sub[lld_probe.sub==0]=NA
  lld_probe.sub=na.omit(lld_probe.sub)
  
  ibd_dosage.sub=ibd_dosage[rownames(ibd_dosage) %in% rownames(ibd_probe.sub),]
  lld_dosage.sub=lld_dosage[rownames(lld_dosage) %in% rownames(lld_probe.sub),]
  
  ibd_probe.sub=ibd_probe.sub[order(rownames(ibd_probe.sub)),]
  ibd_dosage.sub=ibd_dosage.sub[order(rownames(ibd_dosage.sub)),,drop=F]
  lld_probe.sub=lld_probe.sub[order(rownames(lld_probe.sub)),]
  lld_dosage.sub=lld_dosage.sub[order(rownames(lld_dosage.sub)),,drop=F]
  y.list=list(as.vector(ibd_probe.sub),as.vector(lld_probe.sub))
  
  ibd_pheno.sub=ibd_pheno[1:nrow(ibd_dosage.sub),,drop=F]
  lld_pheno.sub=lld_pheno[1:nrow(lld_dosage.sub),,drop=F]
  colnames(ibd_pheno.sub)="pheno"
  colnames(lld_pheno.sub)="pheno"
  x.list=list(as.matrix(ibd_pheno.sub),as.matrix(lld_pheno.sub))
  
  z_sample=all_info[all_info$gene==variant,]
  lld_z_sample=lld_dosage.sub[,colnames(lld_dosage.sub) %in% z_sample$Name,drop=F]
  ibd_z_sample=ibd_dosage.sub[,colnames(ibd_dosage.sub) %in% z_sample$Name,drop=F]
  lld_z_sample=as.data.frame(t(lld_z_sample))
  ibd_z_sample=as.data.frame(t(ibd_z_sample))
  z_sample=merge(ibd_z_sample,lld_z_sample,all = T,by="row.names")
  rownames(z_sample)=z_sample$Row.names
  z_sample=z_sample[,-1]
  z_sample=as.data.frame(t(z_sample))
  z_sample[is.na(z_sample)]=2
  z_sample <- data.frame(lapply(z_sample,as.numeric),check.names = F,row.names = rownames(z_sample))
  
  obj<-Meta_Null_Model(y.list,x.list, n.cohort=2, out_type="C")
  meta.sub=MetaSKAT_wZ(as.matrix(z_sample), obj,r.corr=1)
  return.string=data.frame(PValue=meta.sub$p.value,Genes=variant,Feature=probe)
  
}
write.table(burden_meta_test_2,file = "Burden_meta_raw_test_2.txt",quote = F,row.names = F,sep = "\t")
cutoff=0.05/980
meta_sig_test_2=burden_meta_test_2[burden_meta_test_2$PValue<cutoff,]
meta_sig_test_2=na.omit(meta_sig_test_2)
write.table(meta_sig_test_2,file = "Meta_significant_test_2.txt",quote = F,row.names = F,sep = "\t")
```

# Step.4 Plot figures for checking

```
meta_sig=read.table("Meta_significant.txt",sep = "\t",header = T,stringsAsFactors = F)
# change 0 to 2 in some snp loci
ibd_num=as.numeric(nrow(ibd_dosage))
ibd_dosage[is.na(ibd_dosage)]="QQ"
for(i in 1:ncol(ibd_dosage)){
  
  ibd_dosage.sub=na.omit(ibd_dosage[,i])
  ibd_dosage.sub=as.numeric(ibd_dosage.sub[ibd_dosage.sub!="QQ"])
  sum_dosage=sum(ibd_dosage.sub)
  if(sum_dosage<800){
    
    ibd_dosage[ibd_dosage[,i]==2,i]="M"
    ibd_dosage[ibd_dosage[,i]==0,i]="m"
    ibd_dosage[ibd_dosage[,i]=="M",i]=0
    ibd_dosage[ibd_dosage[,i]=="m",i]=2
    
  }
}
ibd_dosage[ibd_dosage=="QQ"]=NA

lld_num=as.numeric(nrow(lld_dosage))
lld_dosage[is.na(lld_dosage)]="QQ"
for(i in 1:ncol(lld_dosage)){
  
  lld_dosage.sub=na.omit(lld_dosage[,i])
  lld_dosage.sub=as.numeric(lld_dosage.sub[lld_dosage.sub!="QQ"])
  sum_dosage=sum(lld_dosage.sub)
  if(sum_dosage<800){
    
    lld_dosage[lld_dosage[,i]==2,i]="M"
    lld_dosage[lld_dosage[,i]==0,i]="m"
    lld_dosage[lld_dosage[,i]=="M",i]=0
    lld_dosage[lld_dosage[,i]=="m",i]=2
    
  }
}
lld_dosage[lld_dosage=="QQ"]=NA

lld_dosage_num <- data.frame(lapply(lld_dosage,as.numeric),check.names = F)
rownames(lld_dosage_num)=rownames(lld_dosage)
ibd_dosage_num <- data.frame(lapply(ibd_dosage,as.numeric),check.names = F)
rownames(ibd_dosage_num)=rownames(ibd_dosage)

pdf("Burden_9_factors_2_arounds.pdf")
foreach(i=1:nrow(meta_sig),.combine = rbind) %do%  {
  
  gene.sub=meta_sig$Genes[i]
  feature.sub=meta_sig$Feature[i]
  p.sub=meta_sig$PValue[i]
  info.sub=all_info[all_info$gene==gene.sub,]
  for(n in 1:nrow(info.sub)){
    snp.sub=paste("snp",n,sep="_")
    assign(snp.sub,info.sub$Name[n])
    
  }
  print(paste(gene.sub,feature.sub,sep = "---"))
  ibd_geno=ibd_info[ibd_info$gene==gene.sub,,drop=F]
  aa=as.numeric(nrow(ibd_geno))
  ibd_geno=ibd_dosage_num[,ibd_geno$Name,drop=F]
  ibd_geno$dosage=rowSums(ibd_geno)
  ibd_geno$Genotype[ibd_geno$dosage==2*aa]=gene.sub
  ibd_geno$Genotype[ibd_geno$dosage!=2*aa]=paste(gene.sub,"*",sep = "")
  
  ibd_pro=as.data.frame(ibd_probe[,colnames(ibd_probe)==feature.sub],row.names = rownames(ibd_probe))
  ibd_pro[,1][which(ibd_pro[,1]==0)]=NA
  ibd_pro=na.omit(ibd_pro)
  ibd_pro[,1]=scale(ibd_pro[,1])
  ibd_violin=merge(ibd_pro,ibd_geno,all=F,by="row.names")
  ibd_violin=na.omit(ibd_violin)
  factor_1=factor(ibd_violin$Genotype)
  
  p1=ggplot(ibd_violin, aes(factor_1,ibd_violin[,2])) + 
    geom_boxplot()+
    theme_bw()+
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize = 6,alpha=0.4,binwidth=0.01,aes(fill=factor_1))+
    guides(fill=FALSE)+ylab(feature.sub)+xlab("")+
    annotate("text", x=0.8, y=4, label=paste("N=",nrow(total),sep = ""),size=5)+
    annotate("text", x=0.8, y=3.6, label=paste("PValue=",round(P.sub,5),sep = ""),size=5)
  
  lld_geno=lld_info[lld_info$gene==gene.sub,,drop=F]
  bb=as.numeric(nrow(lld_geno))
  lld_geno=lld_dosage_num[,lld_geno$Name,drop=F]
  lld_geno$dosage=rowSums(lld_geno)
  lld_geno$Genotype[lld_geno$dosage==2*bb]=gene.sub
  lld_geno$Genotype[lld_geno$dosage!=2*bb]=paste(gene.sub,"*",sep = "")
  
  lld_pro=as.data.frame(lld_probe[,colnames(lld_probe)==feature.sub],row.names = rownames(lld_probe))
  lld_pro[,1][which(lld_pro[,1]==0)]=NA
  lld_pro=na.omit(lld_pro)
  lld_pro[,1]=scale(lld_pro[,1])
  lld_violin=merge(lld_pro,lld_geno,all=F,by="row.names")
  lld_violin=na.omit(lld_violin)
  factor_2=factor(lld_violin$Genotype)
  
  p2=ggplot(lld_violin, aes(factor_2,lld_violin[,2])) + 
    geom_boxplot()+
    theme_bw()+
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize = 6,alpha=0.4,binwidth=0.01,aes(fill=factor_2))+
    guides(fill=FALSE)+ylab("")+xlab("")
  
  colnames(ibd_violin)[c(1,2)]=c("sample","Feature")
  colnames(lld_violin)[c(1,2)]=c("sample","Feature")
  violin=rbind(ibd_violin[,c(1,2,ncol(ibd_violin))],lld_violin[,c(1,2,ncol(lld_violin))])
  violin=na.omit(violin)
  factor=factor(violin$Genotype)
  p3=ggplot(violin, aes(factor,violin[,2])) + 
    geom_boxplot()+
    theme_bw()+
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize = 6,alpha=0.4,binwidth=0.01,aes(fill=factor))+
    guides(fill=FALSE)+ylab("")+xlab("")
  print(multiplot(p1,p2,p3,cols = 3))
}
dev.off()
```


