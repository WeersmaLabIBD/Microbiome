QMP Absolute quantification project
=====================================

*Creator: Anne Boddeke* 
*Year: 2018* 

1. Read table 
-----------------------------------------------
```
# Read table
IBD_samples_all <- read.table("IBD_taxonomy_metaphlan2.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
```

2. Correct and simplify the IDs of the file
---------------------------------------------
```
#### Correct and simplify IDs
### Remove "_metaphlan"
colnames(IBD_samples_all) <- gsub("_metaphlan", "", colnames(IBD_samples_all))
### Rename G IDs to IBDFEC IDs
## Read table with ID conversions
my_id=read.table("./ids.txt", header = T)
## Correct the ID table for missing sample G105251
my_id_correct <- my_id[-c(62),]
# Select data with old IDs
IBD_samples_subset_old_id <- subset(IBD_samples_all, select = c(1:105)) 
# In subset, correct IDs
colnames(IBD_samples_subset_old_id) <- my_id_correct$Classic
# Merge with remaining data
IBD_samples_subset_remain <- subset(IBD_samples_all, select = c(106:544))
IBD_samples_correct_id <- cbind(IBD_samples_subset_old_id, IBD_samples_subset_remain)
# Just because, change "UMCGFECDNA" to "IBDFEC" -> numbers remain the same
colnames(IBD_samples_correct_id) <- gsub("UMCGFECDNA", "IBDFEC", colnames(IBD_samples_correct_id))
# IBD_samples_correct_id is now the table with all samples with correct ids ALL starting with "IBDFEC" 
```

3. Select a subset of patients 
--------------------------------
```
### Select subset for project
# Import table with IDs of the subset
my_subset_ids <- read.table("./ID subset project.txt", header = T)
# Transpose the IBD_samples_correct_id table so it's in the same configuration as the ID table
IBD_samples_correct_id_trans <- t(IBD_samples_correct_id)
IBD_samples_correct_id_trans <- as.data.frame(t(IBD_samples_correct_id))
### Of the subset IDs table, make the IDs the rownames
# Make an extra column
my_subset_ids$keep="KEEP"
# Make IDs the rownames
rownames(my_subset_ids) = my_subset_ids$New_ID
# Remove the column where the IDs were
my_subset_ids$New_ID=NULL
# Merge the tables 
IBD_sample_subset_final = merge(my_subset_ids, IBD_samples_correct_id_trans, by="row.names")
# Remove the "KEEP" column
IBD_sample_subset_final$keep=NULL
# So now there is a subset of only the data of the samples for this preject: IBD_samples_subset
```


