QMP Absolute quantification project using FISH data
=====================================

*Creator: Anne Boddeke* 
*Year: 2018* 

1. Select the samples 
-----------------------------------------------
```
# Read table
FISH_seq_samples_all <- read.table("metaphlanmerged.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)

#### Correct and simplify IDs
### Remove "_metaphlan"
colnames(FISH_seq_samples_all) <- gsub("_metaphlan", "", colnames(FISH_seq_samples_all))
### Remove data of M4 samples
FISH_seq_samples <- FISH_seq_samples_all[,!grepl("_M4$",names(FISH_seq_samples_all))]

### Rename p001_M1 IDs to P1 IDs
## Remove "_M1"
colnames(FISH_seq_samples) <- gsub("_M1", "", colnames(FISH_seq_samples))
## Substitute "p00" for "P0" and "p0" for "P" (as some have a name with "p0") 
colnames(FISH_seq_samples) <- gsub("p00", "P0", colnames(FISH_seq_samples))
colnames(FISH_seq_samples) <- gsub("p0", "P", colnames(FISH_seq_samples))

#### Names are correct now ####

# Remove sample 39 (no FISH data available)
FISH_seq_samples$P39=NULL

# So now there is a subset of only the data of the samples for this preject: FISH_seq_samples

# To easily use it next times, make a txt file 
write.table(FISH_seq_samples, "~/Desktop/Thesis/FISH_seq_samples.txt", sep="\t")


### Add FISH data to the table
# Import the tabel with the counts
FISH_counts <- read.table("RiseUp_Results_Total_Counts.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
# Transpose the sequencing table
t_FISH_seq_samples <- t(FISH_seq_samples)
# Merge tables 
FISH_seq_counts <- merge(t_FISH_seq_samples, FISH_counts, by="row.names")
rownames(FISH_seq_counts) = FISH_seq_counts$Row.names
FISH_seq_counts$Row.names=NULL

# To easily use it next times, make a txt file
write.table(FISH_seq_counts, "~/Desktop/Thesis/FISH_seq_and_counts.txt", sep="\t")
```

2. Make taxonomical tables (species and genus)
-----------------------------------------------
```
library (vegan)

# Read table
FISH_seq_samples <- read.table("FISH_seq_samples.txt", header = T, sep="\t", as.is=TRUE, row.names = 1, check.names = F)
# Transpose table
FISH_seq_samples_all <- t(FISH_seq_samples)
# Look at different taxonomies
taxas = colnames(FISH_seq_samples_all)
mini = head(taxas)
strsplit(mini, "\\|")


###### Splitting taxonomy function ######
### Input you have to give are your taxonomy table and the numbre corresponding to the taxonomical level you want to make a table of
## Taxonomy levels are numeric, in taxonomy name k_|p_|c_|o_|f_|g_|s_|t_, species is 7, order is 4, phylum is 2
# In the end, a .txt table is exported, make sure to give it a new, correct name (a.k.a. put the taxonomical level in the name)

split_taxonomy_level <- function(taxonomy_table,taxonomy_level) {
  
  #Transpose taxonomy table (taxonomies in rows)
  t_tax_table <- as.data.frame(t(taxonomy_table))
  #Create a table to save the taxonomies at specific level. Same size as the original one
  loop_table <- as.data.frame(matrix(nrow = nrow(t_tax_table) , ncol = ncol(t_tax_table)))
  
  ## For each row/taxonomy in the table:
  for (i in 1:nrow(t_tax_table)) {
    # If the taxonomy has the desired number of levels: 
    if (count.fields(textConnection(row.names(t_tax_table[i,])), sep="|") == taxonomy_level){
      #print (paste0("Species found: ", row.names(tax_table[i,]))) ##Loop check
      
      # Only the rows that meet the condition are filled, the rest get NA values
      loop_table[i,] = t_tax_table[i,]
      
    }
  }
  
  #Taxonomy names as rownames (same as in the original transposed table)
  row.names(loop_table) = row.names(t_tax_table)
  #Sample names as colnames (same as in the original table)
  colnames(loop_table) = colnames(t_tax_table)
  
  ##Remove all rows with NA values  
  level_table <- na.omit(loop_table)
  
  ##For each row/taxonomy in the new level table
  for (i in 1:nrow(level_table)){
    #Save only the last part of the taxonomy name = taxonomy level name  
    name_taxonomy <- rownames(level_table)[i]
    new_taxonomy_name <- unlist(strsplit(name_taxonomy, split = "|", fixed= TRUE))[taxonomy_level]
    #Save the new name as rowname
    rownames(level_table)[i] <- new_taxonomy_name
  }
  
  ##Transpose the new level table to get taxonomies in columns and samples in rows
  t_level_table <- as.data.frame(t(level_table))
  
  #Export the new table
  write.table(t_level_table, file = "./level_table.txt", quote = F, sep = "\t")
}

# Perform the function for the different levels
split_taxonomy_level(FISH_seq_samples_all, 2)
# Be sure to do it one level at a time and change the name of the generated table immediately
```

3. Calculate Diversity and perform correlation diversity - counts
-----------------------------------------------
```
# Install and open needed packages
library(vegan)
library(ggplot2)

# Read in level table, for now we'll use species or phylum
# Not yet coupled to counts
species_table <- read.table("FISH_s_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)
# Transpose the table for the diversity function
t_species_table <- as.data.frame(t(species_table))
################## Because outlyer, remove P51 ###################
t_species_table$P51=NULL
species_table1 <- as.data.frame(t(t_species_table))
# Perform diversity function, gives Shannon index values
alpha <- as.data.frame(diversity(species_table1,index="shannon"))
# Write a .txt file of this table
write.table(alpha, file = "~/Desktop/Thesis/Alpha_Shannon_Species_2.txt", quote = F, sep = "\t")

### Visualize the diversity
## Open needed packages
library(dplyr)
library(ggplot2)
## Open the needed table if it's not open yet
alpha_FISH <- read.table("Alpha_Shannon_Species_FISH.txt", header = T, sep="\t", as.is=TRUE, row.names = 1, check.names = F)
## Make the table ready for the plot function
colnames(alpha_FISH) <- "diversity"
bp <- ggplot(data = alpha_FISH, aes(x = "", y = diversity)) + geom_boxplot() + ggtitle("Species alpha diversity of FISH samples") + xlab("") + ylab("Diversity") + theme(plot.title = element_text(size = 20, face = "bold") + axis.title.y = element_text(size = 18)) 
# Save your plot
jpeg(filename = "Boxplot Alpha Diversity Species FISH.jpg")
print(bp)
dev.off()

################################### Correlation ###################################

# Read the table with counts
FISH_counts_1 <- read.table("RiseUp_Results_Total_Counts.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
FISH_counts_2 <- as.data.frame(t(FISH_counts_1))
FISH_counts_2$P51=NULL
FISH_counts <- as.matrix(FISH_counts_1)

# Merge the alpha and counts tables
counts_alpha1 <- merge(FISH_counts, alpha, by="row.names")
rownames(counts_alpha1) <- counts_alpha1[,1]
counts_alpha <- counts_alpha1[,-1]
colnames(counts_alpha)[2] <- "Diversity"

# Perform correlation analysis
library(psych)

for (i in 1:ncol(counts_alpha)){
  
  if (is.numeric(counts_alpha[,i])){
    counts_alpha[,i] <- counts_alpha[,i]
  }
  
  else {
    counts_alpha[,i] <- as.numeric(factor(counts_alpha[,i]))
  }
}

##Calculate the correlation SPEARMAN
cor_alpha_FISH <- corr.test(counts_alpha, method = "spearman", adjust = "fdr")

#Correlation table
corr_alpha_FISH <- as.matrix(cor_alpha_FISH$r)
corr_alpha_FISH[is.na(corr_alpha_FISH)] <- 0
#P_values table
fdr_corr_alpha_FISH <- as.matrix(cor_alpha_FISH$p)
fdr_corr_alpha_FISH[is.na(fdr_corr_alpha_FISH)] <- 1

##### Extra writing tables etc. #####
write.table(corr_alpha_FISH, file = "./Correlation Species Alpha Diversity and counts Spearman.txt", quote = F, sep = "\t")
write.table(fdr_corr_alpha_FISH, file = "./Significance Correlation Species Alpha Diversity and counts Spearman.txt", quote = F, sep = "\t")


################# Visualize your correlation #################
library(ggplot2)
cor_plot <- ggplot(counts_alpha, aes(x = Bacteria_per_gram, y = Diversity)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE, colour="#CC0000") + labs(title = "Correlation Species Alpha Diversity and Bacteria per gram", x = "Bacterial counts per gram feces", y = "Diversity")
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Correlation Species Diversity and Counts 2.jpg")
print(cor_plot)
dev.off()
```


4. Perform correlation number of species - counts
-----------------------------------------------
```
##### Count number of species per sample #####
# Open the FISH species table
species_table <- read.table("FISH_s_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)
################## Because outlyer, remove P51 ###################
t_species_table <- as.data.frame(t(species_table))
t_species_table$P51=NULL
species_table1 <- as.data.frame(t(t_species_table))
### Count the number of species per sample
## Make a presence - absence table (0 if not present, 1 if present)
p_a_table <- species_table1

for (i in 1:ncol(species_table1)) {
  for (j in 1:nrow(species_table1)) {
    if (species_table1[j,i]>0) {
      p_a_table[j,i] = 1
    }
  }
}
## Count the number of times 1 is present per row
# This number (a.k.a. the number of species per sample) is shown in a new, final column
p_a_table$no_calls <- rowSums(p_a_table == 1)
# Make a new table with just the number of species (set row- and column names)
num_species <- as.data.frame(p_a_table$no_calls)
row.names(num_species) <- rownames(p_a_table)
colnames(num_species) <- "no_species"
# Table with the number of species per sample is ready now, called num_species
# Write a .txt table for future use
write.table(num_species, file = "~/Desktop/Thesis/FISH seq data number of species.txt", quote = F, sep = "\t")


### Make a table with counts and number of species for correlation
# Open the counts table in a format suitable for correlation
FISH_counts_1 <- read.table("RiseUp_Results_Total_Counts.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
FISH_counts_2 <- as.data.frame(t(FISH_counts_1))
FISH_counts_2$P51=NULL
FISH_counts <- as.matrix(t(FISH_counts_2))
# Merge the num_species and counts tables
counts_species1 <- merge(FISH_counts, num_species, by="row.names")
rownames(counts_species1) <- counts_species1[,1]
counts_species <- counts_species1[,-1]

############# Do correlation analysis #############
# Perform correlation analysis
library(psych)

for (i in 1:ncol(counts_species)){
  
  if (is.numeric(counts_species[,i])){
    counts_species[,i] <- counts_species[,i]
  }
  
  else {
    counts_species[,i] <- as.numeric(factor(counts_species[,i]))
  }
}

##Calculate the correlation SPEARMAN
cor_species_FISH <- corr.test(counts_species, method = "spearman", adjust = "fdr")

#Correlation table
corr_species_FISH <- as.matrix(cor_species_FISH$r)
corr_species_FISH[is.na(corr_species_FISH)] <- 0
#P_values table
fdr_corr_species_FISH <- as.matrix(cor_species_FISH$p)
fdr_corr_species_FISH[is.na(fdr_corr_species_FISH)] <- 1

##### Extra writing tables etc. #####
write.table(corr_species_FISH, file = "./Correlation Number of Species and counts Spearman.txt", quote = F, sep = "\t")
write.table(fdr_corr_species_FISH, file = "./Significance Correlation Number of Species and counts Spearman.txt", quote = F, sep = "\t")


################# Visualize your correlation #################
library(ggplot2)
cor_plot <- ggplot(counts_species, aes(x = Bacteria_per_gram, y = no_species)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE, colour="#CC0000") + labs(title = "Correlation Number of Species and Bacteria per gram", x = "Bacterial counts per gram feces", y = "Number of Species")
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Correlation Number of Species and Counts 2.jpg")
print(cor_plot)
dev.off()
```

5. Perform correlation relative abundances F.prau and E.coli FISH - Sequencing
-----------------------------------------------
```
### The correlation between the relative abundance of F.prau and E.coli will be looked at
### Especially the correlation F.prau FISH - F.prau seq and E.coli FISH - E.coli seq


########################## F.prau ########################## 
### Get the relative abundances as estimated by FISH
# Open the table with relative abundances of F.prau counted by FISH
rel_FISH <- read.table("Rel_Abu_FISH.txt", sep = "\t", header = T, row.names = 1, check.names = F)
# Select only the F.prau data and change "," to "."
rel_fprau_FISH <- rel_FISH
rel_fprau_FISH$E.coli=NULL
rel_fprau_FISH$F.prau <- as.numeric(gsub(",", ".", rel_fprau_FISH$F.prau))
colnames(rel_fprau_FISH) <- "rel_fprau_FISH"
# Remove sample P51
t_rel_fprau_FISH <- as.data.frame(t(rel_fprau_FISH))
t_rel_fprau_FISH$P51
rel_fprau_FISH1 <- as.data.frame(t(t_rel_fprau_FISH))

### Get the relative abundances as estimated by sequencing
# Open the sequencing species table (called 'FISH' as FISH data is also available of these samples)
species_table <- read.table("FISH_s_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)
################## Because outlyer, remove P51 ###################
t_species_table <- as.data.frame(t(species_table))
t_species_table$P51=NULL
species_table1 <- as.data.frame(t(t_species_table))

# Select only the F.prau column (with correct row-/colnames)
fprau_seq <- as.data.frame(species_table1$s__Faecalibacterium_prausnitzii)
row.names(fprau_seq) <- rownames(species_table1)
colnames(fprau_seq) <- "rel_fprau_seq"

### Make merged table with all data for correlation analysis
fprau_FISH_seq1 <- merge(rel_fprau_FISH1, fprau_seq, by="row.names")
rownames(fprau_FISH_seq1) <- fprau_FISH_seq1[,1]
fprau_FISH_seq <- fprau_FISH_seq1[,-1]

# For more easy future use, write a .txt file of this table
write.table(fprau_FISH_seq, file = "~/Desktop/Thesis/Relative abundance F.prau FISH and Seq 2.txt", quote = F, sep = "\t")


########### Do correlation analysis ###########
# Perform correlation analysis
library(psych)

for (i in 1:ncol(fprau_FISH_seq)){
  
  if (is.numeric(fprau_FISH_seq[,i])){
    fprau_FISH_seq[,i] <- fprau_FISH_seq[,i]
  }
  
  else {
    fprau_FISH_seq[,i] <- as.numeric(factor(fprau_FISH_seq[,i]))
  }
}

##Calculate the correlation SPEARMAN
cor_fprau_FISH_seq <- corr.test(fprau_FISH_seq, method = "spearman", adjust = "fdr")

#Correlation table
corr_fprau_FISH_seq <- as.matrix(cor_fprau_FISH_seq$r)
corr_fprau_FISH_seq[is.na(corr_fprau_FISH_seq)] <- 0
#P_values table
fdr_corr_fprau_FISH_seq <- as.matrix(cor_fprau_FISH_seq$p)
fdr_corr_fprau_FISH_seq[is.na(fdr_corr_fprau_FISH_seq)] <- 1

##### Extra writing tables etc. #####
write.table(corr_fprau_FISH_seq, file = "./Correlation Relative abundance F.prau FISH and sequencing Spearman.txt", quote = F, sep = "\t")
write.table(fdr_corr_fprau_FISH_seq, file = "./Significance Correlation Relative abundance F.prau FISH and sequencing Spearman.txt", quote = F, sep = "\t")


########## Visualize your correlation ###########
library(ggplot2)
cor_plot <- ggplot(fprau_FISH_seq, aes(x = rel_fprau_FISH, y = rel_fprau_seq)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE, colour="#CC0000") + labs(title = "Correlation Relative abundance F.prau FISH and sequencing", x = "Relative abundance F.prau by FISH", y = "Relative abundance F.prau by sequencing")
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Correlation Relative Abundance F.prau by FISH and sequencing 2.jpg")
print(cor_plot)
dev.off()



########################## E.coli ########################## 
# Open the FISH relative abundances and select only the E.coli data
rel_ecoli_FISH <- rel_FISH
rel_ecoli_FISH$F.prau=NULL
rel_ecoli_FISH$E.coli <- as.numeric(gsub(",", ".", rel_ecoli_FISH$E.coli))
colnames(rel_ecoli_FISH) <- "rel_ecoli_FISH"
# Remove sample P51
t_rel_ecoli_FISH <- as.data.frame(t(rel_ecoli_FISH))
t_rel_ecoli_FISH$P51
rel_ecoli_FISH1 <- as.data.frame(t(t_rel_ecoli_FISH))

### Get the relative abundances as estimated by sequencing
# Open the sequencing species table (called 'FISH' as FISH data is also available of these samples)
species_table <- read.table("FISH_s_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)
################## Because outlyer, remove P51 ###################
t_species_table <- as.data.frame(t(species_table))
t_species_table$P51=NULL
species_table1 <- as.data.frame(t(t_species_table))

# Select only the E.coli column
ecoli_seq <- as.data.frame(species_table1$s__Escherichia_coli)
row.names(ecoli_seq) <- rownames(species_table1)
colnames(ecoli_seq) <- "rel_ecoli_seq"


### Make merged table with all data for correlation analysis
ecoli_FISH_seq1 <- merge(rel_ecoli_FISH1, ecoli_seq, by="row.names")
rownames(ecoli_FISH_seq1) <- ecoli_FISH_seq1[,1]
ecoli_FISH_seq <- ecoli_FISH_seq1[,-1]

# For more easy future use, write a .txt file of this table
write.table(ecoli_FISH_seq, file = "~/Desktop/Thesis/Relative abundance E.coli FISH and Seq.txt", quote = F, sep = "\t")


########### Do correlation analysis ###########
# Perform correlation analysis
library(psych)

for (i in 1:ncol(ecoli_FISH_seq)){
  
  if (is.numeric(ecoli_FISH_seq[,i])){
    ecoli_FISH_seq[,i] <- ecoli_FISH_seq[,i]
  }
  
  else {
    ecoli_FISH_seq[,i] <- as.numeric(factor(ecoli_FISH_seq[,i]))
  }
}

##Calculate the correlation SPEARMAN
cor_ecoli_FISH_seq <- corr.test(ecoli_FISH_seq, method = "spearman", adjust = "fdr")

#Correlation table
corr_ecoli_FISH_seq <- as.matrix(cor_ecoli_FISH_seq$r)
corr_ecoli_FISH_seq[is.na(corr_ecoli_FISH_seq)] <- 0
#P_values table
fdr_corr_ecoli_FISH_seq <- as.matrix(cor_ecoli_FISH_seq$p)
fdr_corr_ecoli_FISH_seq[is.na(fdr_corr_ecoli_FISH_seq)] <- 1

##### Extra writing tables etc. #####
write.table(corr_ecoli_FISH_seq, file = "./Correlation Relative abundance E.coli FISH and sequencing Spearman.txt", quote = F, sep = "\t")
write.table(fdr_corr_ecoli_FISH_seq, file = "./Significance Correlation Relative abundance E.coli FISH and sequencing Spearman.txt", quote = F, sep = "\t")


########## Visualize your correlation ###########
library(ggplot2)
cor_plot_2 <- ggplot(ecoli_FISH_seq, aes(x = rel_ecoli_FISH, y = rel_ecoli_seq)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE, colour="#CC0000") + labs(title = "Correlation Relative abundance E.coli FISH and sequencing", x = "Relative abundance E.coli by FISH", y = "Relative abundance E.coli by sequencing")
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Correlation Relative Abundance E.coli by FISH and sequencing.jpg")
print(cor_plot_2)
dev.off()
```

6. Perform correlation sequencing depth - counts
-----------------------------------------------
```
##### Read in sequencing depth data #####
seq_depth_1 <- read.table("QR_Riboflavin.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
### I cheated a bit, as I got a .csv file which I exported to Excel
  # Then I deleted some columns ("PDO Number", "Product", "# Lanes in Aggregation")
  # Also I changed column names so they did'n contain spaces anymore (replace " " for "_")
  # The current read in table contains many many useless columns, so:
# Select only the columns with data in them:
seq_depth <- seq_depth_1[,1:4]

# Make Samples the rownames
rownames(seq_depth) = seq_depth$Sample
# Remove the column where the Samples were and the individual ID column
seq_depth$Sample=NULL
seq_depth$Individual_ID=NULL
# Replace commas
seq_depth$Total_Reads = gsub(",", "", seq_depth$Total_Reads)
seq_depth$Total_Reads = as.numeric(seq_depth$Total_Reads)
seq_depth$PF_Reads = gsub(",", "", seq_depth$PF_Reads)
seq_depth$PF_Reads = as.numeric(seq_depth$PF_Reads)

# Temporarily transform the dataset for selection of samples
t_seq_depth <- t(seq_depth)
t_seq_depth_df <- as.data.frame(t_seq_depth)
### Remove data of M4 samples
library(dplyr)
t_seq_final <- select(t_seq_depth_df, ends_with("_M1"))
### Rename p001_M1 IDs to P1 IDs
## Remove "_M1"
colnames(t_seq_final) <- gsub("_M1", "", colnames(t_seq_final))
## Substitute "p00" for "P0" and "p0" for "P" (as some have a name with "p0") 
colnames(t_seq_final) <- gsub("p00", "P0", colnames(t_seq_final))
colnames(t_seq_final) <- gsub("p0", "P", colnames(t_seq_final))
# Remove column of "P37" (as we don't have counts of this sample)
t_seq_final$P37=NULL
# Transform back!!!
seq_depth_sub <- as.data.frame(t(t_seq_final))
# Select only the total reads
seq_depth_sub$PF_Reads=NULL

# To easily use it next times, make a txt file with just the total reads per individual
write.table(seq_depth_sub, "~/Desktop/Thesis/seq_depth.txt", sep="\t")


##### Read in sequencing depth data #####
seq_depth_1 <- read.table("seq_depth.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
t_seq_depth <- as.data.frame(t(seq_depth_1))
t_seq_depth$P51=NULL
seq_depth <- as.data.frame(t(t_seq_depth))

##### Combine with counts ##### 
### Make a table with counts and total reads (seq_depth) for correlation
# Open the counts table in a format suitable for correlation
FISH_counts_1 <- read.table("RiseUp_Results_Total_Counts.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
FISH_counts_2 <- as.data.frame(t(FISH_counts_1))
FISH_counts_2$P51=NULL
FISH_counts <- as.data.frame(t(FISH_counts_2))

# Merge the num_species and counts tables
counts_depth1 <- merge(FISH_counts, seq_depth, by="row.names")
rownames(counts_depth1) <- counts_depth1[,1]
counts_depth <- counts_depth1[,-1]

############# Do correlation analysis #############
# Perform correlation analysis
library(psych)

for (i in 1:ncol(counts_depth)){
  
  if (is.numeric(counts_depth[,i])){
    counts_depth[,i] <- counts_depth[,i]
  }
  
  else {
    counts_depth[,i] <- as.numeric(factor(counts_depth[,i]))
  }
}

##Calculate the correlation SPEARMAN
cor_counts_depth <- corr.test(counts_depth, method = "spearman", adjust = "fdr")

#Correlation table
corr_counts_depth <- as.matrix(cor_counts_depth$r)
corr_counts_depth[is.na(corr_counts_depth)] <- 0
#P_values table
fdr_corr_counts_depth <- as.matrix(cor_counts_depth$p)
fdr_corr_counts_depth[is.na(fdr_corr_counts_depth)] <- 1

##### Extra writing tables etc. #####
write.table(corr_counts_depth, file = "./Correlation Sequencing Depth and counts Spearman.txt", quote = F, sep = "\t")
write.table(fdr_corr_counts_depth, file = "./Significance Correlation Sequencing Depth and counts Spearman.txt", quote = F, sep = "\t")

################# Visualize your correlation #################
library(ggplot2)
cor_plot <- ggplot(counts_depth, aes(x = Bacteria_per_gram, y = seq_depth)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE, colour="#CC0000") + labs(title = "Correlation Sequencing Depth and Bacteria per gram", x = "Bacterial counts per gram feces", y = "Total reads")
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Correlation Sequencing Depth and Counts 2.jpg")
print(cor_plot)
dev.off()
```

7. Prepare estimated reads data for visualization
-----------------------------------------------
```
# Open table and clear up names and select needed samples
e_reads_table_all <- read.table("estimated_reads_metaphlan.txt", sep = "\t", header = T, row.names = 1, check.names = F)

#### Correct and simplify IDs
### Remove "_metaphlan"
colnames(e_reads_table_all) <- gsub("_metaphlan", "", colnames(e_reads_table_all))
### Remove data of M4 samples
e_reads_table <- e_reads_table_all[,!grepl("_M4$",names(e_reads_table_all))]

### Rename p001_M1 IDs to P1 IDs
## Remove "_M1"
colnames(e_reads_table) <- gsub("_M1", "", colnames(e_reads_table))
## Substitute "p00" for "P0" and "p0" for "P" (as some have a name with "p0") 
colnames(e_reads_table) <- gsub("p00", "P0", colnames(e_reads_table))
colnames(e_reads_table) <- gsub("p0", "P", colnames(e_reads_table))

#### Names are correct now ####

# Remove sample 39 (no FISH data available)
e_reads_table$P39=NULL

# So now there is a subset of only the data of the samples for this preject: e_reads_table
# IDs are columns, bacteria are rows

# To easily use it next times, make a txt file 
write.table(e_reads_table, "~/Desktop/Thesis/estimated_reads_table_corrected.txt", sep="\t")

#################################### Make genus table ####################################
library (vegan)

# Read table 
e_reads_table <- read.table("estimated_reads_table_corrected.txt", header = T, sep="\t", as.is=TRUE, row.names = 1, check.names = F)
# Transpose table
e_reads_table_all <- t(e_reads_table)
# Look at different taxonomies
taxas = colnames(e_reads_table_all)
mini = head(taxas)
strsplit(mini, "\\|")


###### Splitting taxonomy function ######
### Input you have to give are your taxonomy table and the numbre corresponding to the taxonomical level you want to make a table of
## Taxonomy levels are numeric, in taxonomy name k_|p_|c_|o_|f_|g_|s_|t_, species is 7, order is 4, phylum is 2
# In the end, a .txt table is exported, make sure to give it a new, correct name (a.k.a. put the taxonomical level in the name)

split_taxonomy_level <- function(taxonomy_table,taxonomy_level) {
  
  #Transpose taxonomy table (taxonomies in rows)
  t_tax_table <- as.data.frame(t(taxonomy_table))
  #Create a table to save the taxonomies at specific level. Same size as the original one
  loop_table <- as.data.frame(matrix(nrow = nrow(t_tax_table) , ncol = ncol(t_tax_table)))
  
  ## For each row/taxonomy in the table:
  for (i in 1:nrow(t_tax_table)) {
    # If the taxonomy has the desired number of levels: 
    if (count.fields(textConnection(row.names(t_tax_table[i,])), sep="|") == taxonomy_level){
      #print (paste0("Species found: ", row.names(tax_table[i,]))) ##Loop check
      
      # Only the rows that meet the condition are filled, the rest get NA values
      loop_table[i,] = t_tax_table[i,]
      
    }
  }
  
  #Taxonomy names as rownames (same as in the original transposed table)
  row.names(loop_table) = row.names(t_tax_table)
  #Sample names as colnames (same as in the original table)
  colnames(loop_table) = colnames(t_tax_table)
  
  ##Remove all rows with NA values  
  level_table <- na.omit(loop_table)
  
  ##For each row/taxonomy in the new level table
  for (i in 1:nrow(level_table)){
    #Save only the last part of the taxonomy name = taxonomy level name  
    name_taxonomy <- rownames(level_table)[i]
    new_taxonomy_name <- unlist(strsplit(name_taxonomy, split = "|", fixed= TRUE))[taxonomy_level]
    #Save the new name as rowname
    rownames(level_table)[i] <- new_taxonomy_name
  }
  
  ##Transpose the new level table to get taxonomies in columns and samples in rows
  t_level_table <- as.data.frame(t(level_table))
  
  #Export the new table
  write.table(t_level_table, file = "./level_table.txt", quote = F, sep = "\t")
}

# Perform the function for the genus level
split_taxonomy_level(e_reads_table_all, 6)

#### So now there is a table with estimated reads at genus level of the FISH data
#### This table is called 'FISH_g_e_reads_table'
```

8. Make one figure including composition, estimated reads, FISH counts and sequecing depth
-----------------------------------------------
```
# Open needed packages
library("psych")
library("reshape2") 
library("ggplot2")

### Choose five random 'normal' samples, plus outlier P51
sample(1:79, 5)
# [1] 74 61 29 28 66 Are the chosen numbers:
## So use P28, P29, P61, P66, P74 and outlier P51

###################### Relative abundance composition analysis ######################
# NOTE: BEWARE to execute this code sample per sample!!!

##### Perform composision analysis on these samples, on genus level
# Open genus table
genus_table <- read.table("FISH_g_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

################# P28 #################
# Select sample P28
genus_results <- as.data.frame(t(genus_table))
genus_P28 <- as.data.frame(genus_results$P28)
# Create an extra column for the bacterial names
genus_results$bacteria=row.names(genus_results)
row.names(genus_results)=NULL
rownames(genus_P28) = genus_results$bacteria
colnames(genus_P28) = "P28"

###############
# The code from here on works with "genus_table", so change the name of your sample table to "genus_table"
genus_table <- as.data.frame(genus_P28)
### Make a table with the four most common phyla and the other phyla as 'other'
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P28, decreasing = T)[1:4],]
## Sum the means of the top 4
sum_top4_genus <- sum(top4_genus)
## Create a category for the other (non-top-four) means
others_genus <- (100 - sum_top4_genus)

###############
### For subsetting the top 4 bacteria WITH bacteria names
# Create an extra column for the bacteria names
genus_table$bacteria=row.names(genus_table)
# Remove the row names (since bacteria names are in separate column now)
row.names(genus_table)=NULL
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P28, decreasing = T)[1:4],]
## Change rownames back to bacteria genus
row.names(top4_genus) = top4_genus$bacteria
## Remove extra column
top4_genus$bacteria=NULL

# Create one table with both the top four and the 'others' category 
genus_results_P28 <- rbind(top4_genus, others_genus)
rownames(genus_results_P28) [5] <- "Others"
# Final table is 'genus_results_P28' 

### Write a .txt table of the Final (1x) table WITH CORRECT NAME
write.table(genus_results_P28, file = "~/Desktop/Thesis/genus_results_P28.txt", sep = "\t", quote = F)


################# P29 #################
# Open genus table
genus_table <- read.table("FISH_g_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Select sample P29
genus_results <- as.data.frame(t(genus_table))
genus_P29 <- as.data.frame(genus_results$P29)
# Create an extra column for the bacterial names
genus_results$bacteria=row.names(genus_results)
row.names(genus_results)=NULL
rownames(genus_P29) = genus_results$bacteria
colnames(genus_P29) = "P29"

###############
# The code from here on works with "genus_table", so change the name of your sample table to "genus_table"
genus_table <- as.data.frame(genus_P29)
### Make a table with the four most common phyla and the other phyla as 'other'
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P29, decreasing = T)[1:4],]
## Sum the means of the top 4
sum_top4_genus <- sum(top4_genus)
## Create a category for the other (non-top-four) means
others_genus <- (100 - sum_top4_genus)

###############
### For subsetting the top 4 bacteria WITH bacteria names
# Create an extra column for the bacteria names
genus_table$bacteria=row.names(genus_table)
# Remove the row names (since bacteria names are in separate column now)
row.names(genus_table)=NULL
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P29, decreasing = T)[1:4],]
## Change rownames back to bacteria genus
row.names(top4_genus) = top4_genus$bacteria
## Remove extra column
top4_genus$bacteria=NULL

# Create one table with both the top four and the 'others' category 
genus_results_P29 <- rbind(top4_genus, others_genus)
rownames(genus_results_P29) [5] <- "Others"
# Final table is 'genus_results_P29' 

### Write a .txt table of the Final (1x) table WITH CORRECT NAME
write.table(genus_results_P29, file = "~/Desktop/Thesis/genus_results_P29.txt", sep = "\t", quote = F)

################# P61 #################
# Open genus table
genus_table <- read.table("FISH_g_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Select sample P61
genus_results <- as.data.frame(t(genus_table))
genus_P61 <- as.data.frame(genus_results$P61)
# Create an extra column for the bacterial names
genus_results$bacteria=row.names(genus_results)
row.names(genus_results)=NULL
rownames(genus_P61) = genus_results$bacteria
colnames(genus_P61) = "P61"

###############
# The code from here on works with "genus_table", so change the name of your sample table to "genus_table"
genus_table <- as.data.frame(genus_P61)
### Make a table with the four most common phyla and the other phyla as 'other'
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P61, decreasing = T)[1:4],]
## Sum the means of the top 4
sum_top4_genus <- sum(top4_genus)
## Create a category for the other (non-top-four) means
others_genus <- (100 - sum_top4_genus)

###############
### For subsetting the top 4 bacteria WITH bacteria names
# Create an extra column for the bacteria names
genus_table$bacteria=row.names(genus_table)
# Remove the row names (since bacteria names are in separate column now)
row.names(genus_table)=NULL
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P61, decreasing = T)[1:4],]
## Change rownames back to bacteria genus
row.names(top4_genus) = top4_genus$bacteria
## Remove extra column
top4_genus$bacteria=NULL

# Create one table with both the top four and the 'others' category 
genus_results_P61 <- rbind(top4_genus, others_genus)
rownames(genus_results_P61) [5] <- "Others"
# Final table is 'genus_results_P61' 

### Write a .txt table of the Final (1x) table WITH CORRECT NAME
write.table(genus_results_P61, file = "~/Desktop/Thesis/genus_results_P61.txt", sep = "\t", quote = F)

################# P66 #################
# Open genus table
genus_table <- read.table("FISH_g_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Select sample P66
genus_results <- as.data.frame(t(genus_table))
genus_P66 <- as.data.frame(genus_results$P66)
# Create an extra column for the bacterial names
genus_results$bacteria=row.names(genus_results)
row.names(genus_results)=NULL
rownames(genus_P66) = genus_results$bacteria
colnames(genus_P66) = "P66"

###############
# The code from here on works with "genus_table", so change the name of your sample table to "genus_table"
genus_table <- as.data.frame(genus_P66)
### Make a table with the four most common phyla and the other phyla as 'other'
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P66, decreasing = T)[1:4],]
## Sum the means of the top 4
sum_top4_genus <- sum(top4_genus)
## Create a category for the other (non-top-four) means
others_genus <- (100 - sum_top4_genus)

###############
### For subsetting the top 4 bacteria WITH bacteria names
# Create an extra column for the bacteria names
genus_table$bacteria=row.names(genus_table)
# Remove the row names (since bacteria names are in separate column now)
row.names(genus_table)=NULL
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P66, decreasing = T)[1:4],]
## Change rownames back to bacteria genus
row.names(top4_genus) = top4_genus$bacteria
## Remove extra column
top4_genus$bacteria=NULL

# Create one table with both the top four and the 'others' category 
genus_results_P66 <- rbind(top4_genus, others_genus)
rownames(genus_results_P66) [5] <- "Others"
# Final table is 'genus_results_P66' 

### Write a .txt table of the Final (1x) table WITH CORRECT NAME
write.table(genus_results_P66, file = "~/Desktop/Thesis/genus_results_P66.txt", sep = "\t", quote = F)

################# P74 #################
# Open genus table
genus_table <- read.table("FISH_g_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Select sample P74
genus_results <- as.data.frame(t(genus_table))
genus_P74 <- as.data.frame(genus_results$P74)
# Create an extra column for the bacterial names
genus_results$bacteria=row.names(genus_results)
row.names(genus_results)=NULL
rownames(genus_P74) = genus_results$bacteria
colnames(genus_P74) = "P74"

###############
# The code from here on works with "genus_table", so change the name of your sample table to "genus_table"
genus_table <- as.data.frame(genus_P74)
### Make a table with the four most common phyla and the other phyla as 'other'
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P74, decreasing = T)[1:4],]
## Sum the means of the top 4
sum_top4_genus <- sum(top4_genus)
## Create a category for the other (non-top-four) means
others_genus <- (100 - sum_top4_genus)

###############
### For subsetting the top 4 bacteria WITH bacteria names
# Create an extra column for the bacteria names
genus_table$bacteria=row.names(genus_table)
# Remove the row names (since bacteria names are in separate column now)
row.names(genus_table)=NULL
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P74, decreasing = T)[1:4],]
## Change rownames back to bacteria genus
row.names(top4_genus) = top4_genus$bacteria
## Remove extra column
top4_genus$bacteria=NULL

# Create one table with both the top four and the 'others' category 
genus_results_P74 <- rbind(top4_genus, others_genus)
rownames(genus_results_P74) [5] <- "Others"
# Final table is 'genus_results_P74' 

### Write a .txt table of the Final (1x) table WITH CORRECT NAME
write.table(genus_results_P74, file = "~/Desktop/Thesis/genus_results_P74.txt", sep = "\t", quote = F)

################# P51 #################
# Open genus table
genus_table <- read.table("FISH_g_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Select sample P51
genus_results <- as.data.frame(t(genus_table))
genus_P51 <- as.data.frame(genus_results$P51)
# Create an extra column for the bacterial names
genus_results$bacteria=row.names(genus_results)
row.names(genus_results)=NULL
rownames(genus_P51) = genus_results$bacteria
colnames(genus_P51) = "P51"

###############
# The code from here on works with "genus_table", so change the name of your sample table to "genus_table"
genus_table <- as.data.frame(genus_P51)
### Make a table with the four most common phyla and the other phyla as 'other'
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P51, decreasing = T)[1:4],]
## Sum the means of the top 4
sum_top4_genus <- sum(top4_genus)
## Create a category for the other (non-top-four) means
others_genus <- (100 - sum_top4_genus)

###############
### For subsetting the top 4 bacteria WITH bacteria names
# Create an extra column for the bacteria names
genus_table$bacteria=row.names(genus_table)
# Remove the row names (since bacteria names are in separate column now)
row.names(genus_table)=NULL
## Find the four most abundant phyla
top4_genus <- genus_table [order(genus_table$P51, decreasing = T)[1:4],]
## Change rownames back to bacteria genus
row.names(top4_genus) = top4_genus$bacteria
## Remove extra column
top4_genus$bacteria=NULL

# Create one table with both the top four and the 'others' category 
genus_results_P51 <- rbind(top4_genus, others_genus)
rownames(genus_results_P51) [5] <- "Others"
# Final table is 'genus_results_P51' 

### Write a .txt table of the Final (1x) table WITH CORRECT NAME
write.table(genus_results_P51, file = "~/Desktop/Thesis/genus_results_P51.txt", sep = "\t", quote = F)


##################################################################
##### Make one plot with 6 stacked bargraphs
### Read all the tables made by code above
P28 <- read.table("genus_results_P28.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
P29 <- read.table("genus_results_P29.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
P61 <- read.table("genus_results_P61.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
P66 <- read.table("genus_results_P66.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
P74 <- read.table("genus_results_P74.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
P51 <- read.table("genus_results_P51.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)

# Merge tables 
### Make one table that contains the composition of all above samples
## all = T makes sure that during merging no genus is lost
genus_composition_1 <- merge(P28, P29, by = "row.names", all = T)
row.names(genus_composition_1)=genus_composition_1$Row.names
genus_composition_1$Row.names=NULL

genus_composition_2 <- merge(genus_composition_1, P61, by = "row.names", all = T)
row.names(genus_composition_2)=genus_composition_2$Row.names
genus_composition_2$Row.names=NULL

genus_composition_3 <- merge(genus_composition_2, P66, by = "row.names", all = T)
row.names(genus_composition_3)=genus_composition_3$Row.names
genus_composition_3$Row.names=NULL

genus_composition_4 <- merge(genus_composition_3, P74, by = "row.names", all = T)
row.names(genus_composition_4)=genus_composition_4$Row.names
genus_composition_4$Row.names=NULL

genus_composition_5 <- merge(genus_composition_4, P51, by = "row.names", all = T)
row.names(genus_composition_5)=genus_composition_5$Row.names
genus_composition_5$Row.names=NULL

### Write a .txt table WITH CORRECT NAME
write.table(genus_composition_5, file = "~/Desktop/Thesis/FISH genus composition samples.txt", sep = "\t", quote = F)

############# Make the barplot #############
## Read in the table with the data
genus_composition <- read.table("FISH genus composition samples.txt", sep = "\t", header = T, row.names = 1, check.names = F)

## Adjust the table for making a stacked barplot
# Improve rownames
rownames(genus_composition) <- gsub("g__", "", rownames(genus_composition))
# Create an extra column for the bacteria names
genus_composition$bacteria=row.names(genus_composition)
# Remove the row names (since bacteria names are in separate column now)
row.names(genus_composition)=NULL
# Make a new table that puts the two groups in one column
my_table = melt(genus_composition)

### Make a stacked barplot ###
FISH_genus_barplot <- ggplot (my_table, aes(x=variable, y=value)) + geom_bar (aes(fill = bacteria), stat = "identity") + theme_classic() + xlab("Sample") + ylab("Relative abundance") + theme(legend.text=element_text(size=12)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold")) + ggtitle("Genus composition \n individual samples") + theme(plot.title = element_text(size=20, face="bold"))
# Have a look!
print(FISH_genus_barplot)
# Warning message that two rows containing missing values were removed: this is because the top four contained different phyla for the two groups. Ignore
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Genus composition individual samples FISH.jpg")
print(FISH_genus_barplot)
dev.off()


############################################ Estimated number of reads analysis ############################################
# Read the genus FISH table of estimated reads
genus_e_reads_table <- read.table("FISH_g_e_reads_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

### Select only the samples you want
### We want P28, P29, P61, P66, P74 and outlier P51
# Transform dataset
t_g_e_reads <- as.data.frame(t(genus_e_reads_table))
# Select columns
library(dplyr)
selection_tge_reads <- select(t_g_e_reads, P28, P29, P61, P66, P74, P51)
#Transform back
selection_g_e_reads <- as.data.frame(t(selection_tge_reads))
### Table 'selection_g_e_reads' now has the estimated number of reads per genus
### for only our samples, ID as rows, genus as columns

# Calculate Total number of estimated reads, a.k.a. the sum per row
tot_e_reads <- as.data.frame(rowSums(selection_g_e_reads))
colnames(tot_e_reads)[1] <- "total_e_reads"
tot_e_reads$ID=row.names(tot_e_reads)

##### Barplot of total estimated reads #####
# Make barplot using ggplot2
tot_e_reads_bp <- ggplot (tot_e_reads, aes(x = ID, y = total_e_reads)) + geom_bar(stat = "identity", fill = "purple") + theme_classic() + xlab("Sample") + ylab("Estimated total \n number of reads") + theme(legend.text=element_text(size=15)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold")) + ggtitle("Estimated total \n number of reads") + theme(plot.title = element_text(size=20, face="bold")) + scale_x_discrete(limits=c("P28", "P29", "P61", "P66", "P74", "P51"))
print(tot_e_reads_bp)


############################################ Combine the two graphs ############################################
# http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/

install.packages("ggpubr")
library(ggpubr)

# Figure with relative abundance genus level and estimated total reads
total_figure <- ggarrange(FISH_genus_barplot, tot_e_reads_bp, labels = c("A", "B"), ncol = 1, nrow = 2, align = "hv", legend = "right", common.legend = T)
print(total_figure)

# EXPORT THIS NOW
jpeg(filename = "Genus composition FISH samples + estimated reads.jpg")
print(total_figure)
dev.off()

############# Add FISH counts and sequencing depth #############

# Open needed packages
library("psych")
library("reshape2") 
library("ggplot2")
library("dplyr")

########################## Select samples FISH counts and sequencing depth ##############################
###################### FISH counts ######################
# Open table with FISH counts
FISH_counts <- read.table("RiseUp_Results_Total_Counts.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)

# Select samples P28, P29, P61, P66, P74 and outlier P51
t_FISH <- as.data.frame(t(FISH_counts))
t_selection_FISH <- select(t_FISH, P28, P29, P61, P66, P74, P51)
selection_FISH_counts <- as.data.frame(t(t_selection_FISH))
# Improve row- and colnames
selection_FISH_counts$ID=row.names(selection_FISH_counts)
### The table to work with from now on is 'selection_FISH_counts' ###
# Write a table for combining figures
write.table(selection_FISH_counts, file = "~/Desktop/Thesis/Fish counts random samples.txt", sep = "\t", quote = F)


###################### Sequencing depth ######################
# Open table with total reads
seq_depth <- read.table("seq_depth.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Select samples P28, P29, P61, P66, P74 and outlier P51
t_depth <- as.data.frame(t(seq_depth))
t_selection_depth <- select(t_depth, P28, P29, P61, P66, P74, P51)
selection_seq_depth <- as.data.frame(t(t_selection_depth))
# Improve row- and colnames
selection_seq_depth$ID=row.names(selection_seq_depth)
### The table to work with from now on is 'selection_seq_depth' ###
# Write a table for combining figures
write.table(selection_seq_depth, file = "~/Desktop/Thesis/Sequencing depth random samples.txt", sep = "\t", quote = F)


########################## Make the barplots ##########################
###################### FISH counts ######################
# Make barplot using ggplot2
FISH_counts_bp <- ggplot (selection_FISH_counts, aes(x = ID, y = Bacteria_per_gram)) + geom_bar(stat = "identity", fill = "deepskyblue2") + theme_classic() + xlab("Sample") + ylab("Total number of bacteria \n per gram feces") + theme(legend.text=element_text(size=15)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold")) + ggtitle("Total bacterial count per gram feces \n by FISH") + theme(plot.title = element_text(size=20, face="bold")) + scale_x_discrete(limits=c("P28", "P29", "P61", "P66", "P74", "P51"))
# Have a look!
print(FISH_counts_bp)
# Warning message that two rows containing missing values were removed: this is because the top four contained different phyla for the two groups. Ignore
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Total FISH counts random samples.jpg")
print(FISH_counts_bp)
dev.off()

###################### Sequencing depth ######################
# Make barplot using ggplot2
seq_depth_bp <- ggplot (selection_seq_depth, aes(x = ID, y = Total_Reads)) + geom_bar(stat = "identity", fill = "chartreuse2") + theme_classic() + xlab("Sample") + ylab("Total number of reads") + theme(legend.text=element_text(size=15)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold")) + ggtitle("Sequencing depth") + theme(plot.title = element_text(size=20, face="bold")) + scale_x_discrete(limits=c("P28", "P29", "P61", "P66", "P74", "P51"))
# Have a look!
print(seq_depth_bp)
# Warning message that two rows containing missing values were removed: this is because the top four contained different phyla for the two groups. Ignore
# Export as jpeg (or pdf, or whatever you like) MAKE SURE TO CHANGE THE FILE NAME
jpeg(filename = "Sequencing depth random samples.jpg")
print(seq_depth_bp)
dev.off()


###### Total figure as above PLUS FISH counts and sequencing depth ######
# Open two tables and make accompanying barplots
selection_FISH_counts <- read.table("Fish counts random samples.txt", sep = "\t", header = T, row.names = 1, check.names = F)
selection_seq_depth <- read.table("Sequencing depth random samples.txt", sep = "\t", header = T, row.names = 1, check.names = F)
FISH_counts_bp <- ggplot (selection_FISH_counts, aes(x = ID, y = Bacteria_per_gram)) + geom_bar(stat = "identity", fill = "deepskyblue2") + theme_classic() + xlab("Sample") + ylab("Total number of bacteria \n per gram feces") + theme(legend.text=element_text(size=15)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold")) + ggtitle("Total bacterial count per gram feces \n by FISH") + theme(plot.title = element_text(size=20, face="bold")) + scale_x_discrete(limits=c("P28", "P29", "P61", "P66", "P74", "P51"))
seq_depth_bp <- ggplot (selection_seq_depth, aes(x = ID, y = Total_Reads)) + geom_bar(stat = "identity", fill = "chartreuse2") + theme_classic() + xlab("Sample") + ylab("Total number of reads") + theme(legend.text=element_text(size=15)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold")) + ggtitle("Sequencing depth") + theme(plot.title = element_text(size=20, face="bold")) + scale_x_discrete(limits=c("P28", "P29", "P61", "P66", "P74", "P51"))

# Make total figure with relative abundance genus level, estimated total reads, FISH counts, and sequencing depth
total_figure <- ggarrange(FISH_genus_barplot, tot_e_reads_bp, FISH_counts_bp, seq_depth_bp, labels = c("A", "B", "C", "D"), hjust = 1, ncol = 1, nrow = 4, align = "hv", heights = 5, legend = "top", common.legend = T)
print(total_figure)

# EXPORT THIS 
pdf("total_figure.pdf", width = 8, height = 11)
ggarrange(FISH_genus_barplot, tot_e_reads_bp, FISH_counts_bp, seq_depth_bp, labels = c("A", "B", "C", "D"), hjust = 1, ncol = 1, nrow = 4, align = "hv", heights = 5, legend = "top", common.legend = T)
dev.off()
```


