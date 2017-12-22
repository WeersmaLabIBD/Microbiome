#!bash/bin/Rscript

#Creator Shixian
#this script is to extract the species, calculate the mean and SD value, also count the number of samples for each species

#set two args
args=commandArgs(T)

if(length(args)!=2){
  cat("Usage:composition_analysis.r Dir+input_metaphlan_txt Dir+output \n
      -----But you need to check the input file form if it is match the script \n")
  quit()
}

library('stringr')

inputDir=args[1]
outputDir=args[2]



input=read.csv(file = inputDir, header=T,sep = '\t')

#this line depends
input=input[-1,]

#write the col.names to output file
ID=input[,1]
n=length(ID)
write.table(input[0,],file = outputDir, sep = ',')

#extract species according to how many '__' in the name
for(i in 2:n){
  check_row=str_count(input[i,1], pattern = '__')
  if(check_row==7){
    write.table(input[i,],file = outputDir, append = T, col.names = F, sep = ",", row.names = F)
  }
  
}

#open the output from above
input3=read.table(file = outputDir, header=T, sep = ',')
input2=input3[,-1]

#caculate mean and SD
row_mean=apply(input2,1,mean)
row_SD=apply(input2,1,var)

#add mean and SD value to new columns 
input3$mean=row_mean
input3$SD=row_SD

#caculate how many samples are present
n=length(input2[,1])
m=length(input2[1,])
presented_sample=c()

for(i in 1:n){
  non_zero=m-length(which(input2[i,]==0))
  presented_sample=append(presented_sample,non_zero)
}

input3$presented_sample=presented_sample

write.table(input3, file = 'last_count.csv',sep = ',',row.names = F)
