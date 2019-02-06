```
# =====================================
# conditional analysis on MYRF
# =====================================
```

ibd_genotype=read.table("conditional.IBD.genotype.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
lld_genotype=read.table("conditional.LLD.genotype.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
rm=read.table("IBD_ID_change_rm_SB.txt",header = F,sep = "\t",stringsAsFactors = F)
ibd_genotype=ibd_genotype[,colnames(ibd_genotype) %in% rm$V2]
ibd_genotype=as.data.frame(t(ibd_genotype))
lld_genotype=as.data.frame(t(lld_genotype))

ibd_probe=read.table("IBD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
lld_probe=read.table("LLD_numeric.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)

library(foreach)

ibd_recalculation = foreach(i=1:nrow(ibd_probe),.combine = rbind) %do%  {
  probe=rownames(ibd_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(ibd_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=ibd_genotype[,"rs2238001",drop=F]
  genotype.sub=genotype.sub[rownames(genotype.sub) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  spearman=cor.test(probe.sub[,probe],genotype.sub$rs2238001,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs2238001",
                             CorrelationCoefficient = spearman$estimate,Pvalue = spearman$p.value,
                             Number=nrow(probe.sub))
}


lld_recalculation = foreach(i=1:nrow(lld_probe),.combine = rbind) %do%  {
  probe=rownames(lld_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(lld_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=lld_genotype[,"rs2238001",drop=F]
  genotype.sub=genotype.sub[rownames(genotype.sub) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  spearman=cor.test(probe.sub[,probe],genotype.sub$rs2238001,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs2238001",
                             CorrelationCoefficient = spearman$estimate,Pvalue = spearman$p.value,
                             Number=nrow(probe.sub))
}

library(metap)
feature=intersect(ibd_recalculation$Probe,lld_recalculation$Probe)
meta_recalculation = foreach(i=1:length(feature),.combine = rbind) %do%  {
  probe=feature[i]
  print(probe)
  ibd.sub=ibd_recalculation[ibd_recalculation$Probe==probe,]
  lld.sub=lld_recalculation[lld_recalculation$Probe==probe,]
  p_vector=c(ibd.sub$Pvalue,lld.sub$Pvalue)
  w_vector=c(ibd.sub$Number,lld.sub$Number)
  mm=sumz(p_vector,weights = w_vector)
  return.string = data.frame(Probe = probe, SNPName = "rs2238001",Pvalue=mm$p)
}
meta_recalculation=meta_recalculation[order(meta_recalculation$Pvalue),]

# ======================
#  partial correlation
# ======================
# regress out rs9735635,rs2240287
library(ppcor)
ibd_partial = foreach(i=1:nrow(ibd_probe),.combine = rbind) %do%  {
  probe=rownames(ibd_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(ibd_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=ibd_genotype[rownames(ibd_genotype) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  x.resid = resid(lm(probe.sub[,probe] ~ genotype.sub$rs9735635+genotype.sub$rs2240287))
  mm=cor.test(x.resid,genotype.sub$rs2238001,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs2238001",
                             CorrelationCoefficient = mm$estimate,Pvalue = mm$p.value,
                             Number=nrow(probe.sub),RegressedSNP="rs9735635,rs2240287")
}
lld_partial = foreach(i=1:nrow(lld_probe),.combine = rbind) %do%  {
  probe=rownames(lld_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(lld_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=lld_genotype[rownames(lld_genotype) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  x.resid = resid(lm(probe.sub[,probe] ~ genotype.sub$rs9735635+genotype.sub$rs2240287))
  mm=cor.test(x.resid,genotype.sub$rs2238001,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs2238001",
                             CorrelationCoefficient = mm$estimate,Pvalue = mm$p.value,
                             Number=nrow(probe.sub),RegressedSNP="rs9735635,rs2240287")
}
meta_partial = foreach(i=1:length(feature),.combine = rbind) %do%  {
  probe=feature[i]
  print(probe)
  ibd.sub=ibd_partial[ibd_partial$Probe==probe,]
  lld.sub=lld_partial[lld_partial$Probe==probe,]
  p_vector=c(ibd.sub$Pvalue,lld.sub$Pvalue)
  w_vector=c(ibd.sub$Number,lld.sub$Number)
  mm=sumz(p_vector,weights = w_vector)
  return.string = data.frame(Probe = probe, SNPName = "rs2238001",mm$p)
}
meta_partial=meta_partial[order(meta_partial$mm.p),]

# regress out rs2238001,rs2240287
ibd_partial = foreach(i=1:nrow(ibd_probe),.combine = rbind) %do%  {
  probe=rownames(ibd_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(ibd_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=ibd_genotype[rownames(ibd_genotype) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  x.resid = resid(lm(probe.sub[,probe] ~ genotype.sub$rs2238001+genotype.sub$rs2240287))
  mm=cor.test(x.resid,genotype.sub$rs9735635,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs9735635",
                             CorrelationCoefficient = mm$estimate,Pvalue = mm$p.value,
                             Number=nrow(probe.sub),RegressedSNP="rs2238001,rs2240287")
}
lld_partial = foreach(i=1:nrow(lld_probe),.combine = rbind) %do%  {
  probe=rownames(lld_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(lld_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=lld_genotype[rownames(lld_genotype) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  x.resid = resid(lm(probe.sub[,probe] ~ genotype.sub$rs2238001+genotype.sub$rs2240287))
  mm=cor.test(x.resid,genotype.sub$rs9735635,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs9735635",
                             CorrelationCoefficient = mm$estimate,Pvalue = mm$p.value,
                             Number=nrow(probe.sub),RegressedSNP="rs2238001,rs2240287")
}
feature=c("PWY-5173_superpathway_of_acetyl-CoA_biosynthesis","TCA-GLYOX-BYPASS_superpathway_of_glyoxylate_bypass_and_TCA")
meta_partial = foreach(i=1:length(feature),.combine = rbind) %do%  {
  probe=feature[i]
  print(probe)
  ibd.sub=ibd_partial[ibd_partial$Probe==probe,]
  lld.sub=lld_partial[lld_partial$Probe==probe,]
  p_vector=c(ibd.sub$Pvalue,lld.sub$Pvalue)
  w_vector=c(ibd.sub$Number,lld.sub$Number)
  mm=sumz(p_vector,weights = w_vector)
  return.string = data.frame(Probe = probe, SNPName = "rs9735635",mm$p)
}
meta_partial=meta_partial[order(meta_partial$mm.p),]

# regress out rs2238001,rs9735635
ibd_partial = foreach(i=1:nrow(ibd_probe),.combine = rbind) %do%  {
  probe=rownames(ibd_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(ibd_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=ibd_genotype[rownames(ibd_genotype) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  x.resid = resid(lm(probe.sub[,probe] ~ genotype.sub$rs2238001+genotype.sub$rs9735635))
  mm=cor.test(x.resid,genotype.sub$rs2240287,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs2240287",
                             CorrelationCoefficient = mm$estimate,Pvalue = mm$p.value,
                             Number=nrow(probe.sub),RegressedSNP="rs2238001,rs9735635")
}
lld_partial = foreach(i=1:nrow(lld_probe),.combine = rbind) %do%  {
  probe=rownames(lld_probe)[i]
  print(probe)
  probe.sub=as.data.frame(t(lld_probe[i,,drop=F]))
  probe.sub[which(probe.sub==0),probe]=NA
  probe.sub=na.omit(probe.sub)
  genotype.sub=lld_genotype[rownames(lld_genotype) %in% rownames(probe.sub),,drop=F]
  probe.sub=probe.sub[rownames(probe.sub) %in% rownames(genotype.sub),,drop=F]
  probe.sub=probe.sub[order(rownames(probe.sub)),,drop=F]
  genotype.sub=genotype.sub[order(rownames(genotype.sub)),,drop=F]
  x.resid = resid(lm(probe.sub[,probe] ~ genotype.sub$rs2238001+genotype.sub$rs9735635))
  mm=cor.test(x.resid,genotype.sub$rs2240287,method = "spearman")
  return.string = data.frame(Probe = probe, SNPName = "rs2240287",
                             CorrelationCoefficient = mm$estimate,Pvalue = mm$p.value,
                             Number=nrow(probe.sub),RegressedSNP="rs2238001,rs9735635")
}
feature=c("PWY-5173_superpathway_of_acetyl-CoA_biosynthesis","TCA-GLYOX-BYPASS_superpathway_of_glyoxylate_bypass_and_TCA")
meta_partial = foreach(i=1:length(feature),.combine = rbind) %do%  {
  probe=feature[i]
  print(probe)
  ibd.sub=ibd_partial[ibd_partial$Probe==probe,]
  lld.sub=lld_partial[lld_partial$Probe==probe,]
  p_vector=c(ibd.sub$Pvalue,lld.sub$Pvalue)
  w_vector=c(ibd.sub$Number,lld.sub$Number)
  mm=sumz(p_vector,weights = w_vector)
  return.string = data.frame(Probe = probe, SNPName = "rs2240287",mm$p)
}

