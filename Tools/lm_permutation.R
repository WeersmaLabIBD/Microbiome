## Permutations: split dataset and perform lm analyses ##

# Creator: Arnau Vich 
# Year: 2018
# Input taxonomical data with some metadata (taxa and metadata in columns x each sample a row)
# This script was made to test which is the probability to get a significant result when subseting a cohort
# Useful to detect any batch effect.
# Needs optimitzation and FDR-correction, so use it for exploratory analyses. 

## Taxonomy Analysis Function ##

# Creator: Arnau Vich | Paula Sureda
# Year: 2017

#Import data
my_data=read.table("taxa.txt", header = T, row.names = 1, sep="\t")

#Select the data that you need (phenos + taxa)
my_data_sel=my_data[c(1,4,5,6,54,56:534)]

#Subset two cohorts
LLD_controls=my_data_sel[my_data_sel$myc=="LLD",]
MTC_controls=my_data_sel[my_data_sel$myc=="MTC",]

#Create empty matrix for results
results_controls=matrix(, nrow = 479,ncol=1002)

names_taxa=colnames(LLD_controls)[6:484]
x=0

#Run 1000 permutations, first we will split one cohort into 132 samples and we will test them against the remaining samples.
# lm function is implemented for the test and Age, Sex, BMI, sequenced reads are added as a covariates

for (i in 1:1000){
  test_set=LLD_controls[sample(nrow(LLD_controls), 132), ]
  remove=row.names(test_set)
  control_set=LLD_controls[!row.names(LLD_controls) %in% remove,]
  test_set$myc="LLD_subset"
  new_test=rbind(control_set,test_set)
  x=0
  for (a in 6:484){
    x=x+1
    result=lm(new_test[,a]~myc+Sex+BMI+AgeAtFecalSampling+PFReads,data = new_test)
    re2=summary(result)
    results_controls[x,i]=re2$coefficients[2,4]
  }
}
row.names(results_controls)=names_taxa

#Repeat the same but now the permutation will be use for testing the same number of samples of two different cohorts. 
results_mibs=matrix(, nrow = 479,ncol=1002)
names_taxa=colnames(LLD_controls)[6:484]

x=0
for (i in 1:1000){
  test_set=LLD_controls[sample(nrow(LLD_controls), 132), ]
  new_test=rbind(test_set,MTC_controls)
  x=0
  for (a in 6:484){
    x=x+1
    result=lm(new_test[,a]~myc+Sex+BMI+AgeAtFecalSampling+PFReads,data = new_test)
    re2=summary(result)
    results_mibs[x,i]=re2$coefficients[2,4]
  }
}
row.names(results_mibs)=names_taxa

#Finally we can plot an histogram of the p-values distributions
xx=melt(results_controls)
xx2=xx[complete.cases(xx), ]
xx2$color="N"
xx2[xx2$value<0.06,]$color="S"
ggplot(xx2, aes(x=value, fill=color))+ geom_histogram(binwidth = 0.005) + theme_classic() + ylim(c(0,100000)) + scale_fill_manual(values = c("grey45", "#F8766D"))
