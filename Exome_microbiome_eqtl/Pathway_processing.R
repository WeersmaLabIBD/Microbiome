###################################################################################################################
############################################# pick cluster ########################################################
###################################################################################################################



# IBD
# import raw pathway data and select most prevelant ones with present rate more than 25%
# calculate relative abundance
# this is used for picking clusters, not used for mbQTL analyses

ibd_path_raw=read.table(file = "IBD_humann2_pathways_uniref90_082017.txt",sep = "\t",header = T,comment.char = "", check.names = FALSE,row.names = 1)
ibd_path=ibd_path_raw[grep("__",rownames(ibd_path_raw),invert=T),]

ibd_path_ID=as.data.frame(unlist(colnames(ibd_path)),stringsAsFactors = F)
colnames(ibd_path_ID)="path_ID"

ibd_id_change=read.table(file = "IBD_allID_change.txt",sep = "\t",header = T)
ibd_id_change$path_ID = sapply(strsplit(as.character(ibd_id_change$Original),'_'), "[", 1)
ibd_id_change$path_ID=paste(ibd_id_change$path_ID,"_kneaddata_merged_Abundance",sep = "")
pathways=ibd_path[,ibd_path_ID$path_ID %in% ibd_id_change$path_ID]

pathways = pathways[rowSums(pathways > 0) > ncol(pathways)*0.25,]
pathways=apply(pathways,2,function(x){
  
  x=x/rowSums(pathways)
  x
})

pathways=as.data.frame(pathways)
pathways[pathways == 0] = NA
pathways = log(pathways)

pathways=t(pathways)
pathways=as.data.frame(pathways,stringsAsFactors = F)
rownames(pathways)=ibd_id_change$Sample_ID

ibd_covariate=read.table("eqtl_IBD_covariate.txt",sep = "\t",header = T)
ibd_covariate=ibd_covariate[order(rownames(ibd_covariate)),]
pathways=pathways[order(rownames(pathways)),]

ibd_pathways_correct = apply(pathways,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = ibd_covariate[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  if(ncol(covariate.subset)==ncol(covariate.subset.matrix)){
    covariate.subset = covariate.subset[,apply(covariate.subset.matrix,2,sd) !=0,drop = FALSE]
  }
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

# in fact, as.dist did nothing here, the distance matrix is exactly based on 1-r2, as.dist is to keep half of the matirx
r_matrix=as.matrix(cor(ibd_pathways_correct,method = "spearman",use="pairwise.complete.obs"))
r2_matrix=as.matrix(1-abs(r_matrix)^2)
dist_matrix=as.dist(r2_matrix)
clusters=hclust(dist_matrix)
plot(clusters,cex=0.1)
clusters_path = data.frame(colnames(ibd_pathways_correct),cutree(clusters,h=0.7))
clusters_path[,1]=as.character(clusters_path[,1])
colnames(clusters_path)=c("path","cluster_index")

features=matrix(nrow = max(clusters_path[,2]),ncol = 2)
colnames(features)=c("feature","cluster_index")

for(i in 1:max(clusters_path[,2])){
  
  cluster.matrix = r2_matrix[,clusters_path[clusters_path[,2]==i,1],drop=F]
  cluster.matrix =as.data.frame(cluster.matrix )
  cluster.matrix = cluster.matrix[colnames(cluster.matrix),]
  cluster.matrix =as.data.frame(cluster.matrix )
  
  if(ncol(cluster.matrix)>2){
  centroid = rowMeans(cluster.matrix)                       
  
  out = data.frame(value = centroid)
  representative=rownames(out)[out$value==min(out$value)]
  
  if(length(representative)>1){
    
    features[i,1]=as.character(representative[1])
    features[i,2]=i
    
  }else{
    
    features[i,1]=representative
    features[i,2]=i
    
  }
  }else {
    
    features[i,1]=clusters_path[clusters_path[,2]==i,1][1]
    features[i,2]=i
  }
}
features=as.data.frame(features)

# pick out representative features
ibd_pathways_correct=ibd_pathways_correct[,colnames(ibd_pathways_correct) %in% features$feature]

corrected_data = as.data.frame(t(ibd_pathways_correct))
corrected_data2 = cbind(rownames(corrected_data),corrected_data)
colnames(corrected_data2)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_data2),
                   Gene = rownames(corrected_data2),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_data2, file = "IBD_path_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "IBD_path_numeric.txt.annot",sep="\t",row.names=F,quote = F)

#generate presence/absence tables
# treat pathway as binary traits

binary_data = apply(pathways,2,function(x){as.integer(!is.na(x))})
dimnames(binary_data) = dimnames(pathways)
binary_data = binary_data[,colSums(binary_data)>0.25*nrow(binary_data)&colSums(binary_data)<0.9*nrow(binary_data)]
binary_data = as.data.frame(binary_data)+100

r_matrix=as.matrix(cor(binary_data,method = "spearman",use="pairwise.complete.obs"))
r2_matrix=as.matrix(1-abs(r_matrix)^2)
dist_matrix=as.dist(r2_matrix)
clusters=hclust(dist_matrix)
plot(clusters,cex=0.1)
clusters_path = data.frame(colnames(binary_data),cutree(clusters,h=0.7))
clusters_path[,1]=as.character(clusters_path[,1])
colnames(clusters_path)=c("path","cluster_index")

features=matrix(nrow = max(clusters_path[,2]),ncol = 2)
colnames(features)=c("feature","cluster_index")
features=as.data.frame(features)
features$feature=as.character(features$feature)

for(i in 1:max(clusters_path[,2])){
  
  cluster.matrix = r2_matrix[,clusters_path[clusters_path[,2]==i,1],drop=F]
  cluster.matrix =as.data.frame(cluster.matrix )
  cluster.matrix = cluster.matrix[colnames(cluster.matrix),]
  cluster.matrix =as.data.frame(cluster.matrix )
  
  if(ncol(cluster.matrix)>2){
    
    centroid = rowMeans(cluster.matrix)                      
    out = data.frame(value = centroid)
    representative=rownames(out)[out$value==min(out$value)]
    
    if(length(representative)>1){
      
      features[i,1]=as.character(representative[1])
      features[i,2]=i
      
    }else{
    
      features[i,1]=representative
      features[i,2]=i
    
    }
  }else {
    
    features[i,1]=clusters_path[clusters_path[,2]==i,1][1]
    features[i,2]=i
  }
}
features=as.data.frame(features)

# pick out representative features
binary_data=binary_data[,colnames(binary_data) %in% features$feature]
binary_data=as.data.frame(t(binary_data),stringsAsFactors = F)

binary_data = cbind(paste0(rownames(binary_data),".binary"),binary_data)
colnames(binary_data)[1] = "-"

binary_annot = data.frame(platform = "RDP",
                          HT12v3.ArrayAddress = binary_data[,1],
                          Gene = binary_data[,1],
                          Chr = 4,
                          ChrStart = 1000,
                          ChrEnd = 1000)

write.table(binary_data, file = "IBD_path_binary.txt",sep="\t",row.names = F,quote = F)
write.table(binary_annot,file = "IBD_path_binary.txt.annot",sep="\t",row.names=F,quote = F)

# LLD
# import raw pathway data and select most prevelant ones with present rate more than 25%
# calculate relative abundance

lld_path_raw=read.table(file = "LLD_humann2_Uniref90_092017.txt",sep = "\t",header = T,comment.char = "", check.names = FALSE,row.names = 1)
lld_path=lld_path_raw[grep("__",rownames(lld_path_raw),invert=T),]

lld_path_ID=as.data.frame(unlist(colnames(lld_path)),stringsAsFactors = F)
colnames(lld_path_ID)="path_ID"

lld_id_change=read.table(file = "LLD_allID_change.txt",sep = "\t",header = T)
lld_id_change$path_ID = sapply(strsplit(as.character(lld_id_change$ID),'_'), "[", 1)
lld_id_change$path_ID=paste(lld_id_change$path_ID,"_kneaddata_merged_Abundance",sep = "")
pathways=lld_path[,lld_path_ID$path_ID %in% lld_id_change$path_ID]

pathways = pathways[rowSums(pathways > 0) > ncol(pathways)*0.25,]
pathways=apply(pathways,2,function(x){
  
  x=x/rowSums(pathways)
  x
})

pathways=as.data.frame(pathways)
pathways[pathways == 0] = NA
pathways = log(pathways)

pathways=t(pathways)
pathways=as.data.frame(pathways,stringsAsFactors = F)
rownames(pathways)=lld_id_change$exome

lld_covariate=read.table("eqtl_LLD_covariate.txt",sep = "\t",header = T)

pathways=pathways[order(rownames(pathways)),]
lld_covariate=lld_covariate[order(rownames(lld_covariate)),]

lld_pathways_correct = apply(pathways,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = lld_covariate[!is.na(x),,drop = FALSE]
  covariate.subset.matrix = data.matrix(covariate.subset)
  if(ncol(covariate.subset)==ncol(covariate.subset.matrix)){
    covariate.subset = covariate.subset[,apply(covariate.subset.matrix,2,sd) !=0,drop = FALSE]
  }
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

# in fact, as.dist did nothing here, the distance matrix is exactly based on 1-r2, as.dist is to keep half of the matirx
r_matrix=as.matrix(cor(lld_pathways_correct,method = "spearman",use="pairwise.complete.obs"))
r2_matrix=as.matrix(1-abs(r_matrix)^2)
dist_matrix=as.dist(r2_matrix)
clusters=hclust(dist_matrix)
plot(clusters,cex=0.1)
clusters_path = data.frame(colnames(lld_pathways_correct),cutree(clusters,h=0.7))
clusters_path[,1]=as.character(clusters_path[,1])
colnames(clusters_path)=c("path","cluster_index")

features=matrix(nrow = max(clusters_path[,2]),ncol = 2)
colnames(features)=c("feature","cluster_index")

for(i in 1:max(clusters_path[,2])){
  
  cluster.matrix = r2_matrix[,clusters_path[clusters_path[,2]==i,1],drop=F]
  cluster.matrix =as.data.frame(cluster.matrix )
  cluster.matrix = cluster.matrix[colnames(cluster.matrix),]
  cluster.matrix =as.data.frame(cluster.matrix )
  
  if(ncol(cluster.matrix)>2){
    centroid = rowMeans(cluster.matrix)                       
    
    out = data.frame(value = centroid)
    representative=rownames(out)[out$value==min(out$value)]
    
    if(length(representative)>1){
      
      features[i,1]=as.character(representative[1])
      features[i,2]=i
      
    }else{
      
      features[i,1]=representative
      features[i,2]=i
      
    }
  }else {
    
    features[i,1]=clusters_path[clusters_path[,2]==i,1][1]
    features[i,2]=i
  }
}
features=as.data.frame(features)

# pick out representative features
lld_pathways_correct=lld_pathways_correct[,colnames(lld_pathways_correct) %in% features$feature]

corrected_data = as.data.frame(t(lld_pathways_correct))
corrected_data2 = cbind(rownames(corrected_data),corrected_data)
colnames(corrected_data2)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_data2),
                   Gene = rownames(corrected_data2),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_data2, file = "LLD_path_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "LLD_path_numeric.txt.annot",sep="\t",row.names=F,quote = F)

#generate presence/absence tables
# treat pathway as binary traits

binary_data = apply(pathways,2,function(x){as.integer(!is.na(x))})
dimnames(binary_data) = dimnames(pathways)
binary_data = binary_data[,colSums(binary_data)>0.25*nrow(binary_data)&colSums(binary_data)<0.9*nrow(binary_data)]
binary_data = as.data.frame(binary_data)+100

r_matrix=as.matrix(cor(binary_data,method = "spearman",use="pairwise.complete.obs"))
r2_matrix=as.matrix(1-abs(r_matrix)^2)
dist_matrix=as.dist(r2_matrix)
clusters=hclust(dist_matrix)
plot(clusters,cex=0.1)
clusters_path = data.frame(colnames(binary_data),cutree(clusters,h=0.7))
clusters_path[,1]=as.character(clusters_path[,1])
colnames(clusters_path)=c("path","cluster_index")

features=matrix(nrow = max(clusters_path[,2]),ncol = 2)
colnames(features)=c("feature","cluster_index")
features=as.data.frame(features)
features$feature=as.character(features$feature)

for(i in 1:max(clusters_path[,2])){
  
  cluster.matrix = r2_matrix[,clusters_path[clusters_path[,2]==i,1],drop=F]
  cluster.matrix =as.data.frame(cluster.matrix )
  cluster.matrix = cluster.matrix[colnames(cluster.matrix),]
  cluster.matrix =as.data.frame(cluster.matrix )
  
  if(ncol(cluster.matrix)>2){
    
    centroid = rowMeans(cluster.matrix)                      
    out = data.frame(value = centroid)
    representative=rownames(out)[out$value==min(out$value)]
    
    if(length(representative)>1){
      
      features[i,1]=as.character(representative[1])
      features[i,2]=i
      
    }else{
      
      features[i,1]=representative
      features[i,2]=i
      
    }
  }else {
    
    features[i,1]=clusters_path[clusters_path[,2]==i,1][1]
    features[i,2]=i
  }
}
features=as.data.frame(features)

# pick out representative features
binary_data=binary_data[,colnames(binary_data) %in% features$feature]
binary_data=as.data.frame(t(binary_data),stringsAsFactors = F)

binary_data = cbind(paste0(rownames(binary_data),".binary"),binary_data)
colnames(binary_data)[1] = "-"

binary_annot = data.frame(platform = "RDP",
                          HT12v3.ArrayAddress = binary_data[,1],
                          Gene = binary_data[,1],
                          Chr = 4,
                          ChrStart = 1000,
                          ChrEnd = 1000)

write.table(binary_data, file = "LLD_path_binary.txt",sep="\t",row.names = F,quote = F)
write.table(binary_annot,file = "LLD_path_binary.txt.annot",sep="\t",row.names=F,quote = F)
