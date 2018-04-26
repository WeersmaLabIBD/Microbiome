### Clustering dendrogram PCoA ###

#Creator: Paula Sureda | Arnau Vich
#Year: 2017

#Input files
#
#1.Taxonomy table at specific level (species, order, family...)
#2.Category table: file contaning in the first column sample_ID and in the second column the category name

# Example tax_level_table
#
#SID    SP1   SP2    SP3
#Sample1  23.5    3.5    5.5
#Sample2    10.9   43.3   34.6
#Sample3    50     10     30

# Example category_table
#
#SID    category
#Sample1  cat1
#Sample2  cat1
#Sample3  cat2



clustering_dendrogram <- function(tax_level_table, category_table) {
  #Packages needed
  library(vegan)
  library(gclus)
  library(data.table)

  ## Calculate distance matrix - Bray method
  beta <- vegdist(tax_level_table, method="bray")
  # Calculate AVERAGE distance in beta
  caver <- hclust(beta, method="aver")
  # Reorder the matrix generated    
  caver1 <- reorder.hclust(caver, beta)
  # Get the labels(sample identification) in the same order as are in caver1    
  order_caver1 = as.data.frame(caver1$labels[c(caver1$order)])  
  setDT(order_caver1, keep.rownames = TRUE)[]

  # Give names to the columns and use ID column to put rownames
  colnames(order_caver1)[1] <- "num"
  colnames(order_caver1)[2] <- "ID"

  # Before merge: create a new colum in category_table with the ID
  category_table$ID <- rownames(category_table)    
  
  ## Merge order_caver1 with the category_table        
  order_cat <- merge(category_table, order_caver1, by="ID")    
  rownames(order_cat) <- order_cat[,1]
  order_cat <- order_cat[,-1]

  # Calculate the number of categories
  category_number <- nlevels(order_cat[,1]) #First column because in merge the category table is the first one
  # Create a new column to assign a number to each category
  order_cat$category <- as.integer(order_cat[,1])
  # Create a new column to colour the plot depending on the category (level)        
  order_cat$color = "none"     
 
  # Create a palette of colors depending on the number of categories
  my_palette <- matrix(brewer.pal(category_number,"Paired"))

  #Assign one different color to each different category
  for (i in 1:category_number){
    order_cat[order_cat$category == i,]$color = my_palette[i,1]
  }

  # Create a new variable to get the same value to plot as barplot        
  order_cat$values <- "5"
  order_cat$values <- as.numeric(as.character(order_cat$values))

  # Sort the table by "num" column (same order as caver1)    
  order_cat1 <- order_cat[order(as.numeric(order_cat$num)),] 

  ### Cluster dendrogram  
  ## Instruction to combine plots: get two plots in the same view    
  par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)

  # Create the dendrogram clustering plot (based on caver1 distances) 
  plot(caver1, hang=-1, labels = FALSE, axes = FALSE, ylab = "", xlab="", sub="") 

  # Create a barplot colored by category (the samples are in the same order as in caver1)   
  cluster_dendrogram <- barplot(order_cat1$values, col=order_cat1$color, border = NA, yaxt="n")

  # Save the plots as pdf.file
  pdf("cluster_dendrogram_plot.pdf")
  print(cluster_dendrogram)
  dev.off()
}
