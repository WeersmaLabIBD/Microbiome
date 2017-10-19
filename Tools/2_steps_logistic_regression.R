
### Logistic regresion - metadata ###

# Creators: Paula Sureda | Arnau Vich

# Year: 2017


#Script to perform logistic regression test
#Steps:
  #1.Metadata quality control: replace NA for median values and remove columns with only one level
  #2.Convert taxonomy table into presence/absence taxonomy table (0,1)
  #3.Model1: Logistic regression with all metadata
  #4.Model2: Select significant factors in step 3 and create a new model with the new factors
  #5.Test model1 and model2
  #6.Output
    #6.1: creates a new table for each combination and each factor
    #6.2: creates a new table for each combination and each factor with significant corrected p.value (<0.05)


#Input files
#
#1.Metadata table
#2.Taxonomy table
#3.column_number: number of the column in the metadata that contains the category factor - numeric value

#Output files
#
#A file for each factor that has an effect on taxonomy variance

#Example metadata_table
#
#SID        factor1     factor2   factor3
#Sample1      3.2         23        no
#Sample2      2.4          3        yes
#Sample3     10.3          5        yes

#Example taxonomy_table
#
#SID        tax1    tax2    tax3
#Sample1    0.01    1.34    10.2
#Sample2    5.6     0.56    50.2
#Sample3    3.2     6.2     2.34

#Example output
#
#Taxonomy      presence_cat1    absence_cat1      presence_cat2     absence_cat2      factor       effect        p-value      corrected p-value     model
#tax1               57                298               2                99             age         0.09          0.0001            0.002           model1
#tax2              125                230              10                91             age         1.3           0.003             0.01            model2
#tax3              335                 20              69                32             age        -0.79          0.00004           0.0009          model2


logistic_regression <- function(metadata_input, taxonomy_table, column_number) {

  #Package needed
  library (psych)

  ##Function to calculate nº of 0
  nzsum <- function(x){
    sum (x==0)
  }

  ##Function to calculate nº of non-0
  nsum <- function(x){
    sum (x!=0)
  }

  #Function to create a table for multiple combinations
  expand.grid.unique <- function(x, y, include.equals=FALSE){
    
    x <- unique(x)
    y <- unique(y)

    g <- function(i){
      z <- setdiff(y, x[seq_len(i-include.equals)])
      if(length(z)) cbind(x[i], z, deparse.level=0)
    }

    do.call(rbind, lapply(seq_along(x), g))
  }


  # Remove NA values
  # Convert categorical values to numeric
  for (i in 1:ncol(metadata_input)) {
    if (is.factor(metadata_input[,i]) & any(is.na(metadata_input[,i]))) {
      metadata_input[,i] <- as.integer(metadata_input[,i])
    }
  }

  # Replace NA values: median value
  for (ii in 1:ncol(metadata_input)){
    for (jj in 1:nrow(metadata_input)) {
      if (is.na(metadata_input[jj,ii])){
        x = describe(metadata_input[,ii])
        a = x$median
        metadata_input[jj,ii] = a
      }
    }
  }

  ##Remove columns with only one level
  metadata_input <- metadata_input[, sapply(metadata_input, function(col) length(unique(col))) > 1]

  #Create presence/absence taxonomy table
  p_a_table <- taxonomy_table

  for (i in 1:ncol(taxonomy_table)) {
    for (j in 1:nrow(taxonomy_table)) {
      if (taxonomy_table[j,i]>0) {
        p_a_table[j,i] = 1
      }
    }
  }

  
  #Multiple combinations
  llista = unique(as.vector(metadata_input[,column_number]))
  combination_list <- expand.grid.unique(llista,llista)

  matrix_list <- list()
  table_variables <- matrix(ncol = 2, nrow = ncol(p_a_table))

  for (aa in 1:nrow(combination_list)){

    new_metadata <- metadata_input[metadata_input[,1]==combination_list[aa,1] | metadata_input[,1] == combination_list[aa,2],]

    ##Remove columns with only one level
    new_metadata1 <- new_metadata[, sapply(new_metadata, function(col) length(unique(col))) > 1]


    #For each taxonomy
    for (x in 1:ncol(p_a_table)) {
      #Get column
      name_column <- colnames(p_a_table)[x]
      taxonomy <- subset(p_a_table, select = name_column)

      #Create a table for the model. Merge taxonomy column with metadata  
      model_table <- merge(taxonomy, new_metadata1, by = "row.names" ) 
      row.names(model_table) <- model_table[,1]
      model_table <- model_table[,-1]
      #Change taxonomy name
      colnames(model_table)[1] <- "Taxonomy" 

      #Sort
      #model_table <- model_table[ , order(names(model_table))]
  
      #Calculate model
      model <- glm(Taxonomy ~ . , family = binomial(link = "logit"), data = model_table)
  
      ##Calculate Anova
      anova_test <- anova(model, test = "Chisq")
  
      ##Keep significative variables for model2
      variables_model2 <- subset(anova_test, anova_test[,5]<0.05)
  
      list_variables_model1 <- as.vector(c(colnames(model_table)))
      list_variables_model1 <- paste(c(list_variables_model1), collapse=',' )


      #If there are significative variables
      if (nrow(variables_model2)>0) {
        # Save names of the significative variables: create a new metadata table with these variables (model_table2)
        matrix_variables_model2 <- as.data.frame(rownames(variables_model2))
        rownames(matrix_variables_model2) <- matrix_variables_model2[,1]
        t_model_table <- t(model_table)

        list_variables_model2 <- as.vector(matrix_variables_model2$`rownames(variables_model2)`)
        list_variables_model2 <- paste(c(list_variables_model2), collapse=', ' )

        t_model_table2 <- merge(matrix_variables_model2, t_model_table, by = "row.names")
        rownames(t_model_table2) <- t_model_table2[,1]
        t_model_table2[1:2] <- NULL

        model_table2 <- t(t_model_table2) 

        # Merge the new table with the taxonomy again (lost in the last merge)
        model_table2 <- merge(taxonomy, model_table2, by = "row.names" ) 
        row.names(model_table2) <- model_table2[,1]
        model_table2 <- model_table2[,-1]
    
        #Change taxonomy name
        colnames(model_table2)[1] <- "Taxonomy"
    
        # Change character to numeric
        for (ii in 1:ncol(model_table2))  {
          column_name2 <- colnames(model_table2)[ii]
      
          for (jjj in 1:ncol(model_table)) {
            column_name <- colnames(model_table)[jjj]
        
            if (column_name == column_name2 & is.numeric(model_table[,jjj])) {
              model_table2[,ii] <- as.numeric(as.character(model_table2[,ii]))
            }
          }
        }
    
        # Sort
        #model_table2 <- model_table2[ , order(names(model_table2))]
    
        ## Calculate model2 with the new table: contains only the significative variables in anova test
        model2 <- glm(Taxonomy ~ . , family = binomial(link="logit"), data = model_table2)
    
        # Test the two models
        anova_test <- anova(model, model2, test = "Chisq")
    
        # Model 2 and model 1 are equal: save model 1 results
        if (is.na(anova_test[2,5])) {
      
          summary_table <- summary(model)  
          coefficients_table <- as.data.frame(summary_table$coefficients)
      
          category1_samples <- subset(model_table, model_table[,2]==combination_list[aa,1])
          nsum1 = nsum(category1_samples$Taxonomy)
          nzsum1 = nzsum(category1_samples$Taxonomy)
      
          category2_samples <- subset(model_table, model_table[,2]==combination_list[aa,2])
          nsum2 = nsum(category2_samples$Taxonomy)
          nzsum2 = nzsum(category2_samples$Taxonomy)
      
          loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
          colnames(loop_matrix) <- c("Taxonomy","presence_cat1", "absence_cat1", "presence_cat2", "absence_cat2", "Variable","effect","p_value", "Model")
      
          a = 1
      
          for (jj in 1:nrow(coefficients_table)) {
            loop_matrix[a,1] = name_column
            loop_matrix[a,2] = nsum1
            loop_matrix[a,3] = nzsum1
            loop_matrix[a,4] = nsum2
            loop_matrix[a,5] = nzsum2
            loop_matrix[a,6] = rownames(coefficients_table)[jj]
            loop_matrix[a,7] = coefficients_table[jj,1]
            loop_matrix[a,8] = coefficients_table[jj,4]
            loop_matrix[a,9] = "model_1"
            a <- a + 1
          }
      
          loop_matrix <- na.omit(loop_matrix)
          loop_matrix1 <- loop_matrix
      
          for (kk in 1:nrow(loop_matrix)) {
            if (loop_matrix[kk,6]=="(Intercept)") {
              loop_matrix1 <- loop_matrix[-kk,]
            }
          }

          matrix_list[[x]] <- loop_matrix1
      
          table_variables[x,1] = name_column
          table_variables[x,2] = list_variables_model1
        }
    
        #Model 2 is not better than model 1: save model 1 results
        else if (anova_test[2,5]<0.05) {

          summary_table <- summary(model)  
          coefficients_table <- as.data.frame(summary_table$coefficients)
      
          category1_samples <- subset(model_table, model_table[,2]==combination_list[aa,1])
          nsum1 = nsum(category1_samples$Taxonomy)
          nzsum1 = nzsum(category1_samples$Taxonomy)
      
          category2_samples <- subset(model_table, model_table[,2]==combination_list[aa,2])
          nsum2 = nsum(category2_samples$Taxonomy)
          nzsum2 = nzsum(category2_samples$Taxonomy)
      
          loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
          colnames(loop_matrix) <- c("Taxonomy","presence_cat1", "absence_cat1", "presence_cat2", "absence_cat2", "Variable","effect","p_value", "Model")
      
          a = 1
      
          for (jj in 1:nrow(coefficients_table)) {
        
            loop_matrix[a,1] = name_column
            loop_matrix[a,2] = nsum1
            loop_matrix[a,3] = nzsum1
            loop_matrix[a,4] = nsum2
            loop_matrix[a,5] = nzsum2
            loop_matrix[a,6] = rownames(coefficients_table)[jj]
            loop_matrix[a,7] = coefficients_table[jj,1]
            loop_matrix[a,8] = coefficients_table[jj,4]
            loop_matrix[a,9] = "model_1"
            a <- a + 1  
          }
      
          loop_matrix <- na.omit(loop_matrix)
          loop_matrix1 <- loop_matrix
      
          for (kk in 1:nrow(loop_matrix)) {
            if (loop_matrix[kk,6]=="(Intercept)") {
              loop_matrix1 <- loop_matrix[-kk,]
            }
          }
      
          matrix_list[[x]] <- loop_matrix1
      
          table_variables[x,1] = name_column
          table_variables[x,2] = list_variables_model1
        }
    
        ##Model 2 is better than model 1: save model 2 results
        else {
          summary_table2 <- summary(model2)  

          #Save coefficients results
          coefficients_table <- as.data.frame(summary_table2$coefficients)

          #Get number presence/absence of the taxonomy for each category
          category1_samples <- subset(model_table, model_table[,2]==combination_list[aa,1])
          nsum1 = nsum(category1_samples$Taxonomy)
          nzsum1 = nzsum(category1_samples$Taxonomy)
      
          category2_samples <- subset(model_table, model_table[,2]==combination_list[aa,2])
          nsum2 = nsum(category2_samples$Taxonomy)
          nzsum2 = nzsum(category2_samples$Taxonomy)

          loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
          colnames(loop_matrix) <- c("Taxonomy","presence_cat1", "absence_cat1", "presence_cat2", "absence_cat2", "Variable","effect","p_value", "Model")

          a = 1
          #Save in a new matrix:
          for (jj in 1:nrow(coefficients_table)) {
            loop_matrix[a,1] = name_column    #name taxonomy
            loop_matrix[a,2] = nsum1          #presence taxonomy in category1
            loop_matrix[a,3] = nzsum1         #absence taxonomy in category1
            loop_matrix[a,4] = nsum2          #presence taxonomy in category2
            loop_matrix[a,5] = nzsum2         #absence taxonomy in category2
            loop_matrix[a,6] = rownames(coefficients_table)[jj]   #name of the variable
            loop_matrix[a,7] = coefficients_table[jj,1]    #effect value
            loop_matrix[a,8] = coefficients_table[jj,4]    #p_value
            loop_matrix[a,9] = "model_2"
            a <- a + 1  
          }

          loop_matrix <- na.omit(loop_matrix)  #remove empty rows (NA values)
          loop_matrix1 <- loop_matrix

          for (kk in 1:nrow(loop_matrix)) {   #remove (Intercept) results
            if (loop_matrix[kk,6]=="(Intercept)") {
              loop_matrix1 <- loop_matrix[-kk,]
            }
          }

          matrix_list[[x]] <- loop_matrix1   #Save the new matrix in a list of matrix

          table_variables[x,1] = name_column
          table_variables[x,2] = list_variables_model2
        }
      }
  
      ## If are not significative variables in anova test: keep model1 results
      else {

        summary_table <- summary(model)  
        coefficients_table <- as.data.frame(summary_table$coefficients)

        category1_samples <- subset(model_table, model_table[,2]==combination_list[aa,1])
        nsum1 = nsum(category1_samples$Taxonomy)
        nzsum1 = nzsum(category1_samples$Taxonomy)
      
        category2_samples <- subset(model_table, model_table[,2]==combination_list[aa,2])
        nsum2 = nsum(category2_samples$Taxonomy)
        nzsum2 = nzsum(category2_samples$Taxonomy)

        loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
        colnames(loop_matrix) <- c("Taxonomy","presence_cat1", "absence_cat1", "presence_cat2", "absence_cat2", "Variable","effect","p_value", "Model")

        a = 1

        for (jj in 1:nrow(coefficients_table)) {

          loop_matrix[a,1] = name_column
          loop_matrix[a,2] = nsum1
          loop_matrix[a,3] = nzsum1
          loop_matrix[a,4] = nsum2
          loop_matrix[a,5] = nzsum2
          loop_matrix[a,6] = rownames(coefficients_table)[jj]
          loop_matrix[a,7] = coefficients_table[jj,1]
          loop_matrix[a,8] = coefficients_table[jj,4]
          loop_matrix[a,9] = "model_1"
          a <- a + 1  
        }
    
        loop_matrix <- na.omit(loop_matrix)
        loop_matrix1 <- loop_matrix
    
        for (kk in 1:nrow(loop_matrix)) {
          if (loop_matrix[kk,6]=="(Intercept)") {
            loop_matrix1 <- loop_matrix[-kk,]
          }
        }
    
        matrix_list[[x]] <- loop_matrix1
    
        table_variables[x,1] = name_column
        table_variables[x,2] = list_variables_model1
      }


      #Save in the same matrix all the results
      all_matrix <- as.data.frame(do.call(rbind, matrix_list))

      colnames(table_variables) <- c("Taxonomy", "Variables")
      write.table(table_variables, file = "./table_variables_logistic_regression.txt", sep = "\t", quote = F)
    }
  
    #Split by Variable   
    split_matrix <- split(all_matrix, all_matrix$Variable)

    for (bb in 1:length(split_matrix)){
      #Correct by p_values
      matrix <- as.data.frame(split_matrix[bb])
      p_value <- as.vector(matrix[,8])
  
      corrected_pvalues <- p.adjust(p_value, method = "fdr")
  
      #Add a new column with the new p_values 
      matrix <- cbind(matrix, corrected_pvalues)
  
      name <- colnames(matrix)[6]
      name <- paste(name, combination_list[aa,1], "vs", combination_list[aa,2], sep = "")
      nc <- paste(name, ".txt", sep ="")
      assign(nc, matrix)
      final_name_matrix <- paste("./", nc, sep = "")
  
      write.table(matrix, file = final_name_matrix, quote = F, sep = "\t")
  
      ## Filtering significant results
      #Filter by the new p_values
      filtered_matrix <- subset(matrix, matrix[,10]<0.05)
  
      name2 <- colnames(filtered_matrix)[6]
      name2 <- paste(name2, combination_list[aa,1],"vs",combination_list[aa,2], sep ="")
      nc2 <- paste(name2, "_filtered.txt", sep ="")
  
      assign(nc2, filtered_matrix)
      final_name_matrix2 <- paste("./", nc2, sep = "")
  
      #Not print empty tables
      if (nrow(filtered_matrix)>0) {
        write.table(filtered_matrix, file = final_name_matrix2, quote = F, sep = "\t" ) 
      }
    }
  }  
}
