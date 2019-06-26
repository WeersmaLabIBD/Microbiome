options = commandArgs(trailingOnly = TRUE)
file = options[1]
out=paste(file,"microbiome.txt",sep=".")

get_taxonomy_table = function(otu_table, replace_string){
  otu_notax = as.matrix(otu_table[,-ncol(otu_table)])
  taxonomy = sub(replace_string,"",otu_table[,ncol(otu_table)])
  dnew = aggregate(otu_notax ~ as.factor(taxonomy),FUN = sum)
  rownames(dnew) = as.character(dnew[,1])
  dnew = dnew[,-1,drop=F]
  dnew
  return(dnew)
}
otus = read.table(file,header=T,row.names=1,sep="\t",as.is = T,check.names=F)
colnames(otus)[1]=paste("Sample",colnames(otus)[1],sep="_")
metadata = data.frame(tax = c("species","genus","family","order","class","phylum"),replace_string = c("; D_7.*","; D_6.*","; D_5.*","; D_4.*","; D_3.*","; D_2.*"))
result = list()
for (i in 1:nrow(metadata)){
  taxonomy_table = get_taxonomy_table(otus,replace_string = metadata[i,2])
  result[[i]] = taxonomy_table
}
final_table =rbind(result[[1]],result[[2]],result[[3]],result[[4]],result[[5]],result[[6]])
colnames(final_table)=colnames(otus)[1]

write.table(final_table,file = out,sep="\t")
