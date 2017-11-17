# Creator: Arnau Vich
# Date: 2016
# Usage: copy / paste in R console


library (ggplot2)
library (vegan)

#Importa data
#Taxa: rows -> samples, columns -> taxonomy 
taxa=read.table("taxa.txt", sep = "\t", header = T, row.names = 1)
#Phenoypes
pheno=read.table("pheno.txt", sep = "\t", header = T, row.names = 1)

#Function to filter non-rendundant taxa
filtering_taxonomy = function(x){
     x = rev(x)
     result = c()
     for(i in 1:length(x)){
         rmstr = grep("t__",x[i])
         if (length(rmstr) >0) next
         check_prev = grep(x[i], result[length(result)])
         if (length(check_prev) == 0) result = append(result,x[i])
     }
     
     rev(result)
}

colnames_taxa2_filtered = filtering_taxonomy(colnames(taxa))
taxa2= taxa[,colnames_taxa2_filtered]

beta<- vegdist(taxa2, method="bray")
mypcoa=cmdscale(beta, k = 5)
mypcoa2=as.data.frame(mypcoa)
df=merge(pheno,mypcoa2,by="row.names")

# Calculate centroids as the mean of the first principal components  
centroids <- aggregate(cbind(V1,V2)~myc,df,mean)

#Color according phenotypes
df$mc="black"
df[df$myc=="MTC" ,]$mc <- "red"
df[df$myc=="LLD" ,]$mc <- "black"

#Make plot 
ggplot(df,aes(V1,-V2, color=mc))  +  geom_point (alpha=0.55) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_point(data=centroids2,size=25, colour=c("black","red"), alpha=0.9) +  scale_color_identity ("Datasets", breaks=df$color, labels=df$Cohort, guide="legend")

#Add ellipse
ggplot(df,aes(V1,-V2, color=mc))  +  geom_point (alpha=0.55) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_point(data=centroids,size=10, colour=c("black","red"), alpha=0.9) + stat_ellipse() +  scale_color_identity ("Datasets", breaks=df$color, labels=df$Cohort, guide="legend")
