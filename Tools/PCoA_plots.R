### Principal coordinate analysis (PCoA) ###

#Creator: Paula Sureda | Arnau Vich

#Year: 2017

#Input files
#
#1.Taxonomy table at specific level (species, order, family...)
#2.pcoa_elements: number of principal components to calculate in the pcoa
#3.variable_table: metadata for coloring pcoa plot
#4.If variable_table is numeric number of the most abundant variants to color the different plots (if not: top_value = 1)


# Example tax_level_table
#
#SID       SP1    SP2    SP3
#Sample1   23.5    3.5    5.5
#Sample2   10.9   43.3   34.6
#Sample3   50     10     30

#Example variable_table
#
#SID          var1   var2
#Sample1      11.3    3.1
#Sample2      22.9   10.4 
#Sample3       3.2    7.3 


plot_pcoa <- function(tax_level_table, pcoa_elements, variable_table, top_value) {
  ##Required packages
  library(vegan)
  library(ggplot2)
  library(psych)
  library(RColorBrewer)

  ## Generate a distance matrix - Bray method
  beta <- vegdist(tax_level_table, method="bray")
  ## cmdscale -> multidimensional scaling of a data matrix (distance matrix) 
  my_pcoa <- as.data.frame(cmdscale(beta, k = pcoa_elements))
  
  # Variable that contains the ".pdf" string to save plots in pdf files 
  a <- ".pdf"  
  
  ## If the variable_table is numeric:
  if (is.numeric(variable_table[,1])) {
    # Variable that contains colors codes to color the plot 
    my_col=c("#0000FF","#62A4D1","#5BE55B","#FFF000", "#FF0000")

    # Merge the variable_table with my_pcoa table --> the tax_level_table and the variable_table need to have the same number of rows/samples
    numeric_new_table <- merge(variable_table, my_pcoa, by="row.names") 
    # First column as Rownames
    rownames(numeric_new_table) <- numeric_new_table[,1]
    numeric_new_table <- numeric_new_table[,-1]
          
    # Create new tables with new row number
    new_variable_table <- subset(numeric_new_table, select = 1:ncol(variable_table))
    new_pcoa <- numeric_new_table
    new_pcoa[1:ncol(variable_table)] <- NULL

    # If the category table has more than one column:
    if (ncol(variable_table)>1) {
      # Calculate the mean for each column of the variable_table
      summary_table <- describe(new_variable_table)
      # Find the rows with the highest mean values (more abundant variables): top_value
      top_tax <- summary_table[order(summary_table$mean, decreasing = T)[1:top_value],]
      # Save the names of the highest mean values
      total_tax_names <- cbind(rownames(top_tax))
      # First column as Rownames
      rownames(total_tax_names) <- total_tax_names[,1]

      # Transpose the category table to merge  
      t_variable_table <- as.data.frame(t(new_variable_table))    
      # Merge the t_variable_table with the total_tax_names (table with the highest mean values) 
      numeric_plot_table <- merge(total_tax_names, t_variable_table, by = "row.names")    
      # First column as Rownames
      rownames(numeric_plot_table) <- numeric_plot_table[,1]
      numeric_plot_table <- numeric_plot_table[,-1]
      # Remove repeated column (rownames)
      numeric_plot_table <- numeric_plot_table[,-1]

      # Transpose the new table to have the variables as columns  
      t_plot_table <- as.data.frame(t(numeric_plot_table))

      write.table(t_plot_table, file = "./pcoa_table.txt", sep = "\t", quote = F)  

      # For loop to get plots as number of variables in t_plot_table 
      for (jj in 1:ncol(t_plot_table)) {
        # Save the name of the variable to title the plot
        name_category <- colnames(t_plot_table)[jj]
        # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file
        name_pdf <- paste(name_category,a, sep = "")
        # Create the plot
        bla = ggplot(new_pcoa, aes(x=V1, y=V2, geom = "blank", colour = t_plot_table[,jj])) + geom_point() + scale_color_gradientn(colours = my_col, colnames(t_plot_table[,jj])) + theme_classic() + labs(x = "PCoA1", y = "PCoA2") + ggtitle(name_category)
          #x_axis: V1 from my_pcoa table (first PCoA element)
          #y_axis: V2 from my_pcoa table (second PCoA element)
          #colour: a different plot is generated for each column in t_plot_table 

        # Create the pdf file
        pdf(name_pdf)
        # Print the plot
        print(bla)
        # Empty the current device to create the next plot in the next loop
        dev.off()
      }
    }

    # The variable_table has only one column:
    else { 
      # Save the name of the variable to title the plot
      name_category <- colnames(new_variable_table)[1]
      # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file   
      name_pdf <- paste(name_category, a, sep = "")
      # Create the plot
      bla = ggplot(new_pcoa, aes(x=V1, y= V2, geom = "blank", colour = new_variable_table[,1])) + geom_point() + scale_color_gradientn(colours = my_col, colnames(variable_table[,1])) + theme_classic() + labs(x="PCoA1", y = "PCoA2") + ggtitle(name_category)
        #x_axis: V1 from my_pcoa table (first PCoA element)
        #y_axis: V2 from my_pcoa table (second PCoA element)
        #colour: colored by the column in the variable_table

      # Create the pdf file
      pdf(name_pdf)
      # Print the plot   
      print(bla)
      # Empty the current device  
      dev.off()
    }
  }
    
  ## The variable_table is categoric:          
  else { 
    # Merge the variable_table with my_pcoa table --> needed the same num of rows in both tables
    categoric_plot_table <- merge(variable_table, my_pcoa, by="row.names") 
    # First column as Rownames
    rownames(categoric_plot_table) <- categoric_plot_table[,1]
    categoric_plot_table <- categoric_plot_table[,-1]

    write.table(categoric_plot_table, file = "./pcoa_table.txt", sep = "\t", quote = F)  

    # Calculate the number of categories
    category_number <- nlevels(categoric_plot_table[,1])
    # Create a new column to assign a number to each category
    categoric_plot_table$category <- as.integer(categoric_plot_table[,1])
    # Create a new column to colour the plot depending on the category (level)        
    categoric_plot_table$color = "none" 
    # Create a palette of colors depending on the number of categories
    my_palette <- matrix(brewer.pal(category_number,"Set1"))

    # For loop to assign one different color to each different category in a new column
    for (i in 1:category_number){
      categoric_plot_table[categoric_plot_table$category == i,]$color = my_palette[i,1]
    }
        
    # Save the name of the variable to title the plot
    name_category <- colnames(variable_table)[1]
    # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file  
    name_pdf <- paste(name_category, a, sep = "")
    
    # Create the plot
    bla = ggplot (categoric_plot_table, aes(x=V1, y=V2, geom="blank", colour=color)) + geom_point () + scale_color_identity("All_categories", breaks=categoric_plot_table$color, labels=categoric_plot_table$category, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
      #x_axis: V1 from my_pcoa table (first PCoA element)
      #y_axis: V2 from my_pcoa table (second PCoA element)
      #colour: colored by the new color column 
      
    # Create the pdf file
    pdf(name_pdf)
    # Print the plot
    print(bla)
    # Empty the current device
    dev.off()
  }
}
