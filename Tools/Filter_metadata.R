
## Create filtered metadata file ##

# Creator: Paula Sureda | Arnau Vich

# Year: 2017

##USAGE##

# Copy or import function to R #

#Input files
#
#1.metadata_table
#2.num_column_reads: number of the column in the metadata that contains the reads count - numeric value
#3.reads_threshold: threshold (value) to filter metadata_table based on reads count - numeric value
#4.na_threshold: percentage of NA (num max of tolerated NA) - numeric value
#5.If filter_var is missing the script does not filter that columns with only one level. Write "yes" if you want to filter also columns with only one level

#Example metadata_table
#
#SID        factor1     factor2   factor3
#Sample1      3.2         23        no
#Sample2      2.4          3        yes
#Sample3     10.3          5        yes

#Example filtered_metadata
#
#SID      factor1     factor3
#Sample1     3.2        no
#Sample3    10.3        yes


filter_metadata <- function(metadata_input, num_column_reads, reads_threshold, na_threshold, filter_var) {
  
  #Calculate threshold based on number of samples

  na_percentage <- (na_threshold*nrow(metadata_input))/100

  if(missing(filter_var)) {
    
    #Create a table with reads_threshold
    filtered_metadata <- as.data.frame(metadata_input[with(metadata_input, metadata_input[,num_column_reads]>reads_threshold),])
  
    ##Remove factors/variables with > na_threshold (number max of tolerated NA)
    cond_na <- sapply(filtered_metadata, function(col) sum(is.na(col)) < na_percentage)
    filtered_metadata <- filtered_metadata[, cond_na, drop = FALSE]
      
    #Write filtered_metadata
    write.table(filtered_metadata, file = "./filtered_metadata.txt", quote = F, sep="\t")
  }

  
  else {

    #Create a table with reads_threshold
    filtered_metadata <- as.data.frame(metadata_input[with(metadata_input, metadata_input[,num_column_reads]>reads_threshold),])

    ##Remove factors/variables with > nq_threshold (number max of tolerated NA)
    cond_na <- sapply(filtered_metadata, function(col) sum(is.na(col)) < na_percentage)
    filtered_metadata <- filtered_metadata[, cond_na, drop = FALSE]
      
    ##Remove columns with only one level
    filtered_metadata <- filtered_metadata[, sapply(filtered_metadata, function(col) length(unique(col))) > 1]
      
    #Write filtered_metadata
    write.table(filtered_metadata, file = "./filtered_metadata.txt", quote = F, sep="\t")
  }    
}
