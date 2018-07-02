# BURDEN TEST
-------------------------------------------

 This is for simple burden test
 
 This procedure is based on a Nature article "A polygenic burden of rare disruptive mutations in schizophrenia"
 
 Creator: Shixian
 

 Burden for LLD cohort
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
      if(length(which(z_sample!=2))>5){
        
        obj=SKAT_Null_Model(as.matrix(lld_probe.sub) ~ 1, out_type = "C")
        outcome=SKAT(as.matrix(z_sample),obj,r.corr=1)
        return.string=data.frame(PValue=outcome$p.value,NMarker=outcome$param$n.marker,VMarker=outcome$param$n.marker.test,Genes=n)
      }
    }else{
      
      return.string=data.frame(PValue=NA,NMarker=NA,VMarker=NA,Genes=n)
    }
  }
  
  return.string = data.frame(Pvalue=sub.string$PValue,NMarker=sub.string$NMarker,VMarker=sub.string$VMarker,Genes=sub.string$Genes,Probe=probe)
  
}
write.table(burden_lld,file = "LLD_burden_raw.txt",row.names = F,sep = "\t",quote = F)
```

Burden for IBD cohort
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
        if(length(which(z_sample!=2))>5){
    
        obj=SKAT_Null_Model(as.matrix(ibd_probe.sub) ~ 1, out_type = "C")
        outcome=SKAT(as.matrix(z_sample),obj,r.corr=1)
        return.string=data.frame(PValue=outcome$p.value,NMarker=outcome$param$n.marker,VMarker=outcome$param$n.marker.test,Genes=n)
    }
    }else{
    
      return.string=data.frame(PValue=NA,NMarker=NA,VMarker=NA,Genes=n)
    }
  }

    return.string = data.frame(Pvalue=sub.string$PValue,NMarker=sub.string$NMarker,VMarker=sub.string$VMarker,Genes=sub.string$Genes,Probe=probe)

}
write.table(burden_ibd,file = "IBD_burden_raw.txt",row.names = F,sep = "\t",quote = F)
```

Extract variants and features for meta analyses
```
awk '{if ($1<0.01) print }' IBD_burden_raw.txt > IBD_burden_significant.txt
awk '{if ($1<0.01) print }' LLD_burden_raw.txt > LLD_burden_significant.txt
```

Burden for meta
```
ibd_sig=read.table("IBD_burden_significant.txt",sep = "\t",header = F,stringsAsFactors = F)
lld_sig=read.table("LLD_burden_significant.txt",sep = "\t",header = F,stringsAsFactors = F)
ibd_sig$name=paste(ibd_sig$V4,ibd_sig$V5,sep = "--")
lld_sig$name=paste(lld_sig$V4,lld_sig$V5,sep = "--")
meta_feature=intersect(ibd_sig$name,lld_sig$name)
nn=strsplit(as.character(meta_feature),"--")

library(foreach)
burden_meta = foreach(i=1:length(nn),.combine = rbind) %do%  {
  
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
    
    obj<-Meta_Null_Model(y.list,x.list, n.cohort=2, out_type="C")
    meta.sub=MetaSKAT_wZ(as.matrix(z_sample), obj,r.corr=1)
    return.string=data.frame(PValue=meta.sub$p.value,Genes=variant,Feature=probe)
  
}
```


