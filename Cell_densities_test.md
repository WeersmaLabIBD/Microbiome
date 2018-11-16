---
title: "Exploring microbial cell quantification"
author: "Arnau Vich"
date: "11/13/2018"
output: html_document
---

# Microbial quantification: total cell density effect 

Load Packages 
--
```{r}
library(vegan)
library (plyr)
library (reshape2)
library (ggplot2)
library (ppcor)
library (ggforce)
library(dplyr)
```

Load data
--

```{r }
#meta_sp=read.table("~/Desktop/QMP/data/metaphlan_sp.txt", sep="\t", header=T, row.names = 1, check.names = F)
#With categorical values
#phenos=read.table("~/Desktop/QMP/phenotypes/db2.txt", sep="\t", heade =T, row.names = 1)
#phenos_num=read.table("~/Desktop/QMP/phenotypes/db2_num.txt", sep="\t", heade =T, row.names = 1)
#species_select=read.table("~/Desktop/QMP/phenotypes/bacteria_of_interest.txt", sep="\t", heade =T, row.names = 1)
#Only numeric values
```

Alternatively load the saved R work space
--

```{r }
load ("~/Desktop/QMP/hopefully_final_results/R_space.RData")
``` 

Calculate Bray-Curtis distances
---

```{r,eval = FALSE}
beta_meta=vegdist(t(meta_sp), method="bray")
#Calculate the first 5 components
my_pcoa_meta=cmdscale(beta_meta, k = 5)
#Transform to dataframe
my_pcoa_meta=as.data.frame(my_pcoa_meta)
#Merge with metadata 
my_pcoa_meta=merge(my_pcoa_meta,phenos, by="row.names")
```
Plot first 2 components, connecting dots per individual (1 ind. = 2 time-points), and colour per individual

```{r, echo=FALSE}
ggplot(my_pcoa_meta,aes(V1,V2, group=Individual_ID))  +  geom_point (size=2, aes(colour=time_point)) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_line( alpha=0.2)
```


Repeat it for other phenotypes, such as location or cell density (FISH counts)

``` {r, echo=FALSE}

ggplot(my_pcoa_meta,aes(V1,V2, group=Individual_ID))  +  geom_point (size=2, aes(colour=as.factor(MontrealL))) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_line( alpha=0.2)

#Remove sample 37 because FISH failed! 
my_pcoa2=my_pcoa_meta[-62,]
my_pcoa2=my_pcoa2[-62,]

ggplot(my_pcoa2,aes(V1,V2, group=Individual_ID))  +  geom_point (size=2, aes(colour=FISH_total)) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_line( alpha=0.2) + scale_color_gradient(low="red", high = "blue")

#Also log transform FISH
ggplot(my_pcoa2,aes(V1,V2, group=Individual_ID))  +  geom_point (size=2, aes(colour=log(FISH_total))) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_line( alpha=0.2) + scale_color_gradient(low="red", high = "blue")
```

Prepare phenotypes for ADONIS test, log transform FISH and Sequencing Depth and remove sample 37
----

```{r,eval = FALSE}
#phenos_num2=phenos_num[-62,]
#phenos_num2=phenos_num2[-62,]
#phenos_num2$FISH_log=log(phenos_num2$FISH_total)
#phenos_num2$RD_log=log(phenos_num2$PF_Reads)

phenos2=phenos[-62,]
phenos2=phenos2[-62,]
phenos2$FISH_log=log(phenos_num2$FISH_total)
phenos2$RD_log=log(phenos_num2$PF_Reads)
```

Prepare microbiome data (species level) for ADONIS test 
----

```{r,eval = FALSE}
#taxa_data=as.data.frame(t(meta_sp))
#Check that phenotypes and microbiome tables are concordant in the number of samples
#taxa_data=taxa_data[rownames(phenos_num2),]
#Calculate (again) Bray Curtis distance matrix
#dist.matrix <- vegdist(taxa_data,method = "bray")

taxa_data=as.data.frame(t(meta_sp))
#Check that phenotypes and microbiome tables are concordant in the number of samples
taxa_data=taxa_data[rownames(phenos2),]
#Calculate (again) Bray Curtis distance matrix
dist.matrix <- vegdist(taxa_data,method = "bray")
```

ADONIS test
---

```{r,eval = FALSE}
#Start an empty matrix

adonis_meta <- matrix(ncol = 7, nrow=ncol(phenos_num2))

#Calculate variance explained per phenotype after 10.000 permutations
for (i in 1:ncol(phenos2)){
  ad<-adonis(dist.matrix ~ phenos2[,i],permutations=10000)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}
adonis_meta= as.data.frame(adonis_meta)
rownames(adonis_meta) = colnames(phenos2)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
```
- First row (Pnumber) is the individual inter-variance, which explains around 80% of the variation

- As expected, ileum resection is in the top 10  

- Cell density (FISH_data) explains ~7%

``` {r, echo=FALSE}
ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + ylab ("Explained variance") + xlab ("Factor")  + theme(text = element_text(size=8))
```

Normalise data (arcsine sqroot for taxa and log for FISH and Read depth)
---

```{r,eval = FALSE}
taxa_data=taxa_data/100
taxa_data= asin(sqrt(taxa_data))
```

Filter taxa prevalent in less than 10% of the samples 
----

```{r,eval = FALSE}

taxa_tmp=as.data.frame(t(taxa_data))
taxa_tmp$prevalence=((rowSums(taxa_tmp!=0))/(ncol(taxa_tmp)))*100
taxa_filt=taxa_tmp[taxa_tmp$prevalence>10,]
taxa_filt$prevalence=NULL
yy=merge(phenos2,as.data.frame(t(taxa_filt)), by="row.names")
```

Test linear models with and without FISH covariate
---- 

```{r,eval = FALSE}
associations=matrix( ncol=4,nrow=174)
associations_fish=matrix( ncol=5,nrow=174)
c=1
for (b in 73:ncol(yy)){

  v=lm(yy[,b]~Age+BMI+Sex+PPI, data=yy)
  a=summary(v)
  #pval Age
  associations[c,1]=a$coefficients[2,4]
  #pval BMI
  associations[c,2]=a$coefficients[3,4]
  #pval Sex
  associations[c,3]=a$coefficients[4,4]
  #pval PPI
  associations[c,4]=a$coefficients[5,4]

  vv=lm(yy[,b]~FISH_log+Age+BMI+Sex+PPI, data=yy)
  aa=summary(vv)

  associations_fish[c,1]=aa$coefficients[2,4]
  #pval Age
  associations_fish[c,2]=aa$coefficients[3,4]
  #pval BMI
  associations_fish[c,3]=aa$coefficients[4,4]
  #pval Sex
  associations_fish[c,4]=aa$coefficients[5,4]
  #pval PPI
  associations_fish[c,5]=aa$coefficients[6,4]
  c=c+1
}

associations=as.data.frame(associations)
rownames(associations)=colnames(yy)[73:245]
colnames(associations)=c("age","bmi","sex","PPI")
associations_fish=as.data.frame(associations_fish)
rownames(associations_fish)=colnames(yy)[73:245]
colnames(associations_fish)=c("fish","age","bmi","sex","PPI")
```

Predict FISH data using Shannon and read depth
----

```{r,eval = FALSE}
#Calculate shannon
my_richness=diversity(taxa_data,index = "shannon")
yy$shannon=my_richness
#Move to the first column
yy2=select(yy,shannon,everything())

lm_pred=lm(FISH_log~shannon+RD_log, data=yy)
yy$FISH_pred=predict(lm_pred)
yy=select(yy,FISH_pred,everything())
```
```{r,echo=FALSE,message=FALSE}

cor.test(yy$FISH_log, yy$FISH_pred, method="spearman")

```

```{r, echo=FALSE}
ggplot(yy, aes(FISH_log, FISH_pred)) + geom_point() + theme_bw() + xlab("Cell density") + ylab("Cell density predicted from Shannon + Read Depth")
```

Predict FISH data using Shannon and DNA concentration
----

```{r,eval = FALSE}
lm_pred2=lm(FISH_total~shannon_meta+Conc, data=zz)
zz$FISH_pred2=predict(lm_pred2)
```

```{r,echo=FALSE, message=FALSE, eval=FALSE}

cor.test(zz$FISH_total, zz$FISH_pred2, method="spearman")

```

```{r, echo=FALSE, eval=FALSE}
ggplot(zz, aes(FISH_total, FISH_pred2)) + geom_point() + theme_bw() + xlab("Cell density") + ylab("Cell density predicted from Shannon + DNA concentration")

```

Correlation bacteria ~ phenotype and partial correlations using FISH data as covariate
----


Subset phenotypes


```{r,eval = FALSE}
##Age, Sex, BMI, Ileal resection, Minimum 2 resectins, PPI,Calprotectin, log_FISH
phenos_subset=phenos_num2[,c(18,35,34,27,26,23,59,67)]
phenos_subset=merge(phenos_subset,as.data.frame(t(taxa_filt)), by="row.names")
row.names(phenos_subset)=phenos_subset$Row.names
phenos_subset$Row.names=NULL
```

Start matrix for results

```{r,eval = FALSE}
correlation_results=matrix( ncol=28,nrow=173)
```

Loop to correlate taxa~pheno and taxa~pheno+FISH

```{r,eval = FALSE}
c=1
for (b in 9:ncol(phenos_subset)){

  #Correlation Age
  x1=cor.test(phenos_subset[,b],phenos_subset[,1], method="spearman")
  correlation_results[c,1]=x1$p.value
  correlation_results[c,2]=x1$estimate
  p1=pcor.test(phenos_subset[,b],phenos_subset[,1],phenos_subset[,8],method="spearman")
  correlation_results[c,3]=p1[1,2]
  correlation_results[c,4]=p1[1,1]

  #Correlation Sex
  x1=cor.test(phenos_subset[,b],phenos_subset[,2], method="spearman")
  correlation_results[c,5]=x1$p.value
  correlation_results[c,6]=x1$estimate
  p1=pcor.test(phenos_subset[,b],phenos_subset[,2],phenos_subset[,8],method="spearman")
  correlation_results[c,7]=p1[1,2]
  correlation_results[c,8]=p1[1,1]

  #Correlation BMI
  x1=cor.test(phenos_subset[,b],phenos_subset[,3], method="spearman")
  correlation_results[c,9]=x1$p.value
  correlation_results[c,10]=x1$estimate
  p1=pcor.test(phenos_subset[,b],phenos_subset[,3],phenos_subset[,8],method="spearman")
  correlation_results[c,11]=p1[1,2]
  correlation_results[c,12]=p1[1,1]

  #Ileal resections(1=colon, 2=ileum)
  x1=cor.test(phenos_subset[,b],phenos_subset[,4], method="spearman")
  correlation_results[c,13]=x1$p.value
  correlation_results[c,14]=x1$estimate
  p1=pcor.test(phenos_subset[,b],phenos_subset[,4],phenos_subset[,8],method="spearman")
  correlation_results[c,15]=p1[1,2]
  correlation_results[c,16]=p1[1,1]

  #Min. 2 resections
  x1=cor.test(phenos_subset[,b],phenos_subset[,5], method="spearman")
  correlation_results[c,17]=x1$p.value
  correlation_results[c,18]=x1$estimate
  p1=pcor.test(phenos_subset[,b],phenos_subset[,5],phenos_subset[,8],method="spearman")
  correlation_results[c,19]=p1[1,2]
  correlation_results[c,20]=p1[1,1]

  #Correlation PPI
  x1=cor.test(phenos_subset[,b],phenos_subset[,6], method="spearman")
  correlation_results[c,21]=x1$p.value
  correlation_results[c,22]=x1$estimate
  p1=pcor.test(phenos_subset[,b],phenos_subset[,6],phenos_subset[,8],method="spearman")
  correlation_results[c,23]=p1[1,2]
  correlation_results[c,24]=p1[1,1]

  #Correlation Calprotectin
  x1=cor.test(phenos_subset[,b],phenos_subset[,7], method="spearman")
  correlation_results[c,25]=x1$p.value
  correlation_results[c,26]=x1$estimate
  p1=pcor.test(phenos_subset[,b],phenos_subset[,7],phenos_subset[,8],method="spearman")
  correlation_results[c,27]=p1[1,2]
  correlation_results[c,28]=p1[1,1]

  c=c+1
}

correlation_results=as.data.frame(correlation_results)
rownames(correlation_results)=colnames(phenos_subset)[9:181]
colnames(correlation_results)=c("age_pval", "age_estimate", "age_fish_pval", "age_fish_estimate", "sex_pval", "sex_estimate", "sex_fish_pval", "sex_fish_estimate","bmi_pval", "bmi_estimate", "bmi_fish_pval", "bmi_fish_estimate","resec_pval", "resec_estimate", "resec_fish_pval", "resec_fish_estimate","min2_pval", "min2_estimate", "min2_fish_pval", "min2_fish_estimate", "PPI_pval", "PPI_estimate", "PPI_fish_pval", "PPI_fish_estimate","Calp_pval", "Calp_estimate", "Calp_fish_pval", "Calp_fish_estimate"  )
```

- PPI effect is large and not influenced by cell density 


```{r ,echo=FALSE, warning=FALSE }
correlation_results$fdr=p.adjust(correlation_results$PPI_fish_pval, method = "BH")
correlation_results$significant="No"
correlation_results[correlation_results$fdr<0.05,]$significant="Yes"
ggplot(correlation_results, aes(PPI_pval,PPI_fish_pval)) + geom_point() + theme_bw() + ylab("pval correlation PPI +FISH") + xlab("pval correlation PPI")  + facet_zoom (x=PPI_pval<0.05, y=PPI_fish_pval<0.05, zoom.size = 0.5, split = F) + scale_color_manual(values=c("black", "red"))
```



- Sex effect is not influenced by cell density 


```{r ,echo=FALSE, warning=FALSE }
correlation_results$fdr=p.adjust(correlation_results$sex_fish_pval, method = "BH")
correlation_results$significant="No"
correlation_results[correlation_results$fdr<0.05,]$significant="Yes"

ggplot(correlation_results, aes(sex_pval,sex_fish_pval, color=significant)) + geom_point() + theme_bw() + ylab("pval correlation Sex +FISH") + xlab("pval correlation Sex (male=1, female=2)")  + facet_zoom (x=sex_pval<0.05, y=sex_fish_pval<0.05, zoom.size = 0.5, split = F)+ scale_color_manual(values=c("black", "red"))
```



- Age show a moderate influence 


```{r ,echo=FALSE, warning=FALSE }
correlation_results$fdr=p.adjust(correlation_results$age_fish_pval, method = "BH")
correlation_results$significant="No"
correlation_results[correlation_results$fdr<0.05,]$significant="Yes"

ggplot(correlation_results, aes(age_pval,age_fish_pval,color=significant)) + geom_point() + theme_bw() + ylab("pval correlation Age +FISH") + xlab("pval correlation Age")  + facet_zoom (x=age_pval<0.05, y=age_fish_pval<0.05, zoom.size = 0.5, split = F) + scale_color_manual(values=c("black", "red"))
```

- Big differences when taking cell density into consideration in patients with Ileal resection! 


```{r ,echo=FALSE, warning=FALSE }
correlation_results$fdr=p.adjust(correlation_results$bmi_fish_pval, method = "BH")
correlation_results$significant="No"
correlation_results[correlation_results$fdr<0.05,]$significant="Yes"

ggplot(correlation_results, aes(bmi_pval,bmi_fish_pval,color=significant)) + geom_point() + theme_bw() + ylab("pval correlation BMI+FISH") + xlab("pval correlation BMI") + facet_zoom (x=bmi_pval<0.05, y=bmi_fish_pval<0.05, zoom.size = 0.5, split = F) + scale_color_manual(values=c("black", "red")) 
```

- Ileal resection show a big influence 


```{r ,echo=FALSE, warning=FALSE }
correlation_results$fdr=p.adjust(correlation_results$resec_fish_pval, method = "BH")
correlation_results$significant="No"
correlation_results[correlation_results$fdr<0.05,]$significant="Yes"

ggplot(correlation_results, aes(resec_pval,resec_fish_pval,color=significant)) + geom_point() + theme_bw() + ylab("pval correlation Ileal resection +FISH") + xlab("pval correlation Ileal resection") + facet_zoom (x=resec_pval<0.05, y=resec_fish_pval<0.05, zoom.size = 0.5, split = F)+ scale_color_manual(values=c("black", "red")) 
```


```{r ,echo=FALSE, warning=FALSE }
correlation_results$fdr=p.adjust(correlation_results$min2_fish_pval, method = "BH")
correlation_results$significant="No"
correlation_results[correlation_results$fdr<0.05,]$significant="Yes"

ggplot(correlation_results, aes(min2_pval,min2_fish_pval,color=significant)) + geom_point() + theme_bw() + ylab("pval correlation More than 1 intestinal resection (yes/no) +FISH") + xlab("pval correlation More than 1 intestinal resection") + facet_zoom (x=min2_pval<0.05, y=min2_fish_pval<0.05, zoom.size = 0.5, split = F)+ scale_color_manual(values=c("black", "red")) 
```


```{r ,echo=FALSE, warning=FALSE, eval=F }
correlation_results$fdr=p.adjust(correlation_results$Calp_fish_pval, method = "BH")
correlation_results$significant="No"
correlation_results[correlation_results$fdr<0.05,]$significant="Yes"
```


```{r ,echo=FALSE, warning=FALSE}
ggplot(correlation_results, aes(Calp_pval,Calp_fish_pval,color=significant)) + geom_point() + theme_bw() + ylab("pval correlation Calprotectin levels +FISH") + xlab("pval correlation Calprotectin levels") + facet_zoom (x=Calp_pval<0.05, y=Calp_fish_pval<0.05, zoom.size = 0.5, split = F)+ scale_color_manual(values=c("black", "red")) 
```

Test each phenotype against FISH
----

```{r,eval = FALSE}
phenos_vs_fish=phenos2
phenos_vs_fish$Individual_ID=NULL
phenos_vs_fish$FISH_total=NULL
phenos_vs_fish$PF_Reads=NULL
phenos_vs_fish$AbsEcoli=NULL
phenos_vs_fish$AbsFprau=NULL
phenos_vs_fish$FISHAveragePreRatioEcoliFprau=NULL
phenos_vs_fish$FISHAveragePreRatioErecEcoli=NULL
phenos_vs_fish$FISHEcoli=NULL
phenos_vs_fish$FISHErec=NULL
phenos_vs_fish$FISHFprau=NULL
phenos_vs_fish$FISH_Rectale=NULL
phenos_vs_fish$Faecalibacterium_prausnitzii_meta=NULL
phenos_vs_fish$Escherichia_coli_meta=NULL
phenos_vs_fish$Eubacterium_rectale_meta=NULL
phenos_vs_fish=select(phenos_vs_fish,FISH_log,everything())
phenos_vs_fish=select(phenos_vs_fish,RD_log,everything())

my_results=matrix(nrow=100, ncol=4)
a=1
for (i in 3:ncol(phenos_vs_fish)){
  if (is.numeric(phenos_vs_fish[,i])){
    my_test=cor.test(phenos_vs_fish[,2], phenos_vs_fish[,i], method="spearman")
    my_results[a,4]=my_test$p.value
    my_results[a,2]=my_test$estimate
    my_results[a,1]=colnames(phenos_vs_fish)[i]
    a=a+1
  } else {
    my_test=pairwise.wilcox.test(phenos_vs_fish[,2], phenos_vs_fish[,i] ,p.adjust.method = "none")
    my_test2=melt(my_test$p.value)
    for (x in 1:nrow(my_test2)){
      my_results[a,4]=as.character(my_test2[x,3])
      my_results[a,3]=as.character(my_test2[x,2])
      my_results[a,2]=as.character(my_test2[x,1])
      my_results[a,1]=colnames(phenos_vs_fish)[i]
      a=a+1  
    }
  }
}
my_results2=as.data.frame(my_results)
my_results2=my_results2[c(1:61),]
my_results2=my_results2[-c(22),]
my_results2=my_results2[-c(18),]
my_results2$V4=as.numeric(as.character(my_results2$V4))
my_results2$FDR=p.adjust(my_results2$V4, method="BH")
```



```{r ,echo=FALSE, warning=FALSE }

ggplot(phenos2, aes(DNA_conc,FISH_total)) + geom_point() + theme_bw()

ggplot(phenos2, aes(Min2Resection,FISH_total)) + geom_boxplot(outlier.size = NA) + theme_bw() + geom_jitter(aes(color=IlealResection)) + xlab ("Patient with at leat 2 resections") + ylab("Cell density") + scale_color_manual(values=c("black", "red")) 

ggplot(phenos2, aes(Min2Resection,FISH_total)) + geom_boxplot(outlier.size = NA) + theme_bw() + geom_jitter(aes(color=Colectomy)) + xlab ("Patient with at leat 2 resections") + ylab("Cell density") + scale_color_manual(values=c("black", "red")) 

ggplot(phenos2, aes(as.factor(liqstool1),FISH_total)) + geom_boxplot(outlier.size = NA) + theme_bw() + xlab("Number liquid stools before sampling") + ylab("Cell density") + geom_jitter()
```


Correlate bacteria to Read Depth and Cell density
-----

```{r,eval = FALSE}
phenos_subset=phenos_num2[,c(18,35,34,27,26,23,59,67,68)]
phenos_subset=merge(phenos_subset,as.data.frame(t(taxa_filt)), by="row.names")
row.names(phenos_subset)=phenos_subset$Row.names
phenos_subset$Row.names=NULL

species_tito=matrix( ncol=4,nrow=173)
 c=1
> for (b in 10:ncol(phenos_subset)){
 
   #Correlation RD
   x1=cor.test(phenos_subset[,b],phenos_subset[,9], method="spearman")
   species_tito[c,1]=x1$p.value
   species_tito[c,2]=x1$estimate
 
   #Correlation FISH
   x1=cor.test(phenos_subset[,b],phenos_subset[,8], method="spearman")
   species_tito[c,3]=x1$p.value
   species_tito[c,4]=x1$estimate
   c=c+1
}
 
species_tito=as.data.frame(species_tito) 
rownames(species_tito)=colnames(phenos_subset)[10:ncol(phenos_subset)]
colnames(species_tito)=c("RD_pvalue", "RD_beta", "FISH_pval", "FISH_beta")
mean_abundance=as.data.frame(colMeans(phenos_subset)[10:ncol(phenos_subset)])
colnames(mean_abundance)=("rel_abund")
species_tito=merge(species_tito,mean_abundance, by="row.names")
```


```{r ,echo=FALSE, warning=FALSE }
ggplot (species_tito, aes(FISH_beta, RD_beta))  + geom_point(aes(size=rel_abund))  + theme_bw() + geom_abline(intercept =0 , slope = 1, color="red") + xlim(c(-0.6,0.6)) + ylim (c(-0.6,0.6)) + xlab("Spearman coefficient bacteria ~ cell density (FISH)") + ylab ("Spearman coefficient bacteria ~ Sequencing depth")
```


Check relation cell counts and selected bacteria
--------


```{r,eval = FALSE}
other_phen=phenos[c(1,6,31)]
species_select=merge(species_select, other_phen, by="row.names")
```

- Fecalibacterium prausnitzii

```{r ,echo=FALSE, warning=FALSE }
ggplot(species_select, aes(FISH_total,FISH_FP_rel, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("F.Prau. relative abundance (FISH)")

ggplot(species_select, aes(FISH_total,Faecalibacterium_prausnitzii_metaphlan, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("F.Prau. relative abundance (MGS)") + ylim(0,20)

ggplot(species_select, aes(FISH_total,FISH_PRAU, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("F.Prau.counts (FISH)")
```

- Eubacterium rectale

```{r ,echo=FALSE, warning=FALSE }
ggplot(species_select, aes(FISH_total,FISH_ER_rel, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("E.Rectale relative abundance (FISH)")

ggplot(species_select, aes(FISH_total,Eubacterium_rectale_metaphlan, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("E.Rectale relative abundance (MGS)") 

ggplot(species_select, aes(FISH_total,FISH_Rectale, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("E.Rectale counts (FISH)")
```


- Escherichia coli

```{r ,echo=FALSE, warning=FALSE }
ggplot(species_select, aes(FISH_total,FISH_EC_rel, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("E.coli relative abundance (FISH)")

ggplot(species_select, aes(FISH_total,Escherichia_coli_metaphlan, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("E.coli relative abundance (MGS)") 

ggplot(species_select, aes(FISH_total,FISH_Coli, group=Individual_ID, color=time_point)) + geom_point() + theme_bw() + geom_line( alpha=0.2, color="black")+ scale_color_manual(values=c("black", "red")) + xlab("Cell density (FISH)") + ylab("E.coli counts (FISH)")
```



Check p-value distribution (change the phenotype you like)
----

```{r,eval = FALSE}
sel_correlation= correlation_results[,c("resec_fish_pval","resec_pval")]
sel_correlation$id=row.names(sel_correlation)
sel2=melt(sel_correlation, id.vars = "id")
my.pvalue.list<-list("Ileal Resection"=sel_correlation$resec_pval, "Ileal Resection FISH"=sel_correlation$resec_fish_pval)
```

```{r ,echo=FALSE, warning=FALSE, eval=FALSE }
ggplot(sel2, aes(value, fill=variable, color=variable)) + geom_histogram(bins=20,position="dodge", alpha =0.7) + theme_bw()
```

- QQplot function from: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R

```{r,eval = FALSE}
library(lattice)
qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
    
    
    #error checking
    if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
    if(!(class(pvalues)=="numeric" || 
         (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
        stop("pvalue vector is not numeric, can't draw plot")
    if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
    if (already.transformed==FALSE) {
        if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
    } else {
        if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
    }
    
    
    grp<-NULL
    n<-1
    exp.x<-c()
    if(is.list(pvalues)) {
        nn<-sapply(pvalues, length)
        rs<-cumsum(nn)
        re<-rs-nn+1
        n<-min(nn)
        if (!is.null(names(pvalues))) {
            grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
            names(pvalues)<-NULL
        } else {
            grp=factor(rep(1:length(pvalues), nn))
        }
        pvo<-pvalues
        pvalues<-numeric(sum(nn))
        exp.x<-numeric(sum(nn))
        for(i in 1:length(pvo)) {
            if (!already.transformed) {
                pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
            } else {
                pvalues[rs[i]:re[i]] <- pvo[[i]]
                exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
            }
        }
    } else {
        n <- length(pvalues)+1
        if (!already.transformed) {
            exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
            pvalues <- -log10(pvalues)
        } else {
            exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
        }
    }
    
    
    #this is a helper function to draw the confidence interval
    panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
        require(grid)
        conf.points = min(conf.points, n-1);
        mpts<-matrix(nrow=conf.points*2, ncol=2)
        for(i in seq(from=1, to=conf.points)) {
            mpts[i,1]<- -log10((i-.5)/n)
            mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
            mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
            mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
        }
        grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
    }
    
    #reduce number of points to plot
    if (should.thin==T) {
        if (!is.null(grp)) {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                      exp.x = round(exp.x, thin.exp.places),
                                      grp=grp))
            grp = thin$grp
        } else {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                      exp.x = round(exp.x, thin.exp.places)))
        }
        pvalues <- thin$pvalues
        exp.x <- thin$exp.x
    }
    gc()
    
    prepanel.qqunif= function(x,y,...) {
        A = list()
        A$xlim = range(x, y)*1.02
        A$xlim[1]=0
        A$ylim = A$xlim
        return(A)
    }
    
    #draw the plot
    xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
           prepanel=prepanel, scales=list(axs="i"), pch=pch,
           panel = function(x, y, ...) {
               if (draw.conf) {
                   panel.qqconf(n, conf.points=conf.points, 
                                conf.col=conf.col, conf.alpha=conf.alpha)
               };
               panel.xyplot(x,y, ...);
               panel.abline(0,1);
           }, par.settings=par.settings, ...
    )
}


```

```{r ,echo=FALSE, warning=FALSE, eval=FALSE }
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
```

