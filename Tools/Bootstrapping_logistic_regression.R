library (psych)
library(plyr)
setwd("~/Desktop/Rebuttal_Trishla/Female_health/")
phenotypes <- read.table("../pheno_2.txt", sep="\t", header=T, row.names = 1)
aro_index <- read.csv("./annotation_all_markers.txt", sep='\t', row.names = 1, header = T)
#remove antibiotic users
phenos=phenotypes[phenotypes$Antibiotics_merged<1,]


## Add second phenotypes with IBS and Menstruations

pheno_male=phenos[phenos$Sex=="Male",]
pheno_female=phenos[phenos$Sex=="Female",]
pheno_female$Sex=NULL
pheno_male$Sex=NULL
pheno_female$Antibiotics_merged=NULL
pheno_male$Antibiotics_merged=NULL
pheno_male$Menstruations_present=NULL
pheno_female2=pheno_female[,c(3,1,2,4)]
pheno_male2=pheno_male[,c(3,1,2)]
pheno_female3=pheno_female2
pheno_female3$Menstruations_present=NULL

markers_RAW <- read.csv("../Female_health/ARG_all_recode_LLD.tsv", sep='\t', header=T)

all_gf=ddply(markers_RAW,"ARO.Accession",numcolwise(sum)) 
row.names(all_gf) <- all_gf$ARO.Accession
all_gf$ARO.Accession <- NULL
all_gf_pre=all_gf
# replace !=0 -> 1
all_gf_pre[all_gf_pre>0]=1 
#############################################################################################################################################################################################
##################################################################################Select only males or females ##############################################################################
#############################################################################################################################################################################################
#results=matrix(nrow=1000, ncol=10)

########MODIFY (uncomment) HERE FOR FEMALE OR MALE ANALYSIS ##################################
#names(pheno_male2)=c("1_IBS", "2_RD", "3_Age")
names(pheno_female3)=c("1_IBS", "2_RD", "3_Age")


########MODIFY (uncomment) HERE FOR FEMALE OR MALE ANALYSIS ##################################
#gf_males=merge (pheno_male2,as.data.frame(t(all_gf_pre)), by="row.names")
gf_males=merge (pheno_female3,as.data.frame(t(all_gf_pre)), by="row.names")

# ncases=24 ncontrols=449
IBS_males=gf_males[gf_males$`1_IBS`==1,]
CON_males=gf_males[gf_males$`1_IBS`==0,]
results=matrix(nrow=1000, ncol=12)
for (t in 1:1000){
 # subset with replacement (bootstrap)
 # Subset ARO table 
 sub_IBS=IBS_males[sample(nrow(IBS_males),size=24,replace=TRUE),]
 sub_CON=CON_males[sample(nrow(CON_males),size=449,replace=TRUE),]
 data=as.data.frame(t(merge (as.data.frame(t(sub_IBS)), as.data.frame(t(sub_CON)), by="row.names")))
 #colnames(data)=data[1,]
 names(data) <- lapply(data[1, ], as.character)
 data=data[-1,]
 data$Row.names=NULL
 phe=data[,1:3]
 aros=data[,-c(1:3)]
 # Filter 5% 
 data_f=aros[,colSums(aros!=0) > round(nrow (aros) * 0.05) ]
 results[t,1]=ncol(data_f)
 results_2=matrix(nrow=ncol(data_f), ncol=11)
 results_2[,1]=colnames(data_f)
 data_f=cbind(phe,data_f)
 runs=1
 for (a in 4:ncol(data_f)){
   
   test_data=data_f[,c(1,2,3,a)]
   results_2[runs,2]=colnames(test_data)[4]
   test_data=as.matrix(test_data)
   mode(test_data) = "numeric"
   test_data=data.frame(test_data)
   test_data=test_data[,c(4,1,2,3)]
   names(test_data)=c("AB","IBS","RD","Age")
   my_out=glm(AB ~. , family = binomial(link="logit"), data=test_data)
   my_out2=summary(my_out)
   results_2[runs,3]=my_out2$coefficients[2,4]
   runs=runs+1
 # Number of AB present 
 # test lm  
 # Number of significant hits 
 # Count 8 AB significance 
  }
  results_2=as.data.frame(results_2)
  results_2$qval=p.adjust(as.numeric(as.character(results_2$V3)),method = "fdr")
  for (i in 1:nrow(results_2)){
    if (results_2[i,1]=="ARO:3002639" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,4]=1
    } else if (results_2[i,1]=="ARO:3002660" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,5]=1
    } else if (results_2[i,1]=="ARO:3000237" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,6]=1
    } else if(results_2[i,1]=="ARO:3002174" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,7]=1
    } else if (results_2[i,1]=="ARO:3002628" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,8]=1
    } else if (results_2[i,1]=="ARO:3002923" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,9]=1
    } else if (results_2[i,1]=="ARO:3002818" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,10]=1
    } else if (results_2[i,1]=="ARO:3000027" && as.numeric(as.character(results_2[i,12])) < 0.05){
        results[t,11]=1
    } 
  }
  results[t,2]=nrow(results_2[as.numeric(as.character(results_2$V3))<0.05,])
  results[t,3]=nrow(results_2[results_2$qval<0.05,])
  
}
