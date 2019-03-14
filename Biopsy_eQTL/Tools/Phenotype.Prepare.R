#============================================================================================================================
# this script is to set the same order of samples of expression table as genotype fam
#============================================================================================================================

#Usage: Phenotype.Prepare.R Normalized_CD/UC.txt CD/UC.plink.fam

options = commandArgs(trailingOnly = TRUE)
table=options[1]
fam=options[2]

normalized=read.table(table,header = T,stringsAsFactors = F,sep = "\t",check.names = F)
normalized_table=as.data.frame(t(normalized[,-1]))
n=normalized$`-`
colnames(normalized_table)=n

coupling=read.table("coupling_file.txt",header = F,check.names = F,stringsAsFactors = F,sep = "\t")
colnames(coupling)=c("exome","biopsy")

normalized_table=normalized_table[rownames(normalized_table) %in% coupling$biopsy,]
coupling=coupling[coupling$biopsy %in% rownames(normalized_table),]

normalized_table=normalized_table[order(rownames(normalized_table)),]
coupling=coupling[order(coupling$biopsy),]
rownames(normalized_table)=coupling$exome

fam=read.table(fam,sep = " ",header = F,check.names = F,stringsAsFactors = F)
normalized_table=normalized_table[order(match(rownames(normalized_table), fam$V2)),]
normalized_table=normalized_table+100
pheno=normalized_table+100

write.table(normalized_table,file = "Reordered.phenotype.txt",sep = "\t",row.names = T,quote = F)
write.table(pheno,file = "Pheno.txt",sep = " ",row.names = F,col.names = F,quote = F)
