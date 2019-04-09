# Author: Alexander Kurilshikov

cutoff_presence = 0.1
get_taxonomy_table = function(otu_table, replace_string,cutoff = 0.1){
  otu_notax = as.matrix(otu_table[,-ncol(otu_table)])
  taxonomy = sub(replace_string,"",otu_table[,ncol(otu_table)])
  dnew = aggregate(otu_notax ~ as.factor(taxonomy),FUN = sum)
  rownames(dnew) = as.character(dnew[,1])
  dnew = dnew[,-1]
  dnew = t(dnew)
  dnew
  filter = dnew[,(colSums(dnew > 0) > 0.1*nrow(dnew))]
  return(dnew)
}
otus = read.table("otu_table.tsv",header=T,row.names=1,sep="\t",as.is = T)
metadata = data.frame(tax = c("genus","family","order","class"),replace_string = c("; D_6.*","; D_5.*","; D_4.*","; D_3.*"))
result = list()
for (i in 1:nrow(metadata)){
  taxonomy_table = get_taxonomy_table(otus,replace_string = metadata[i,2])
  result[[i]] = taxonomy_table
}
final_table = cbind(result[[1]],result[[2]],result[[3]],result[[4]])
final_table = final_table[,(colSums(final_table>0)>cutoff_presence * nrow(final_table))]
final_ln = log(final_table)
final_table[final_table == "-Inf"] = NA
write.table(final_table,file = "microbes.txt",sep="\t")
