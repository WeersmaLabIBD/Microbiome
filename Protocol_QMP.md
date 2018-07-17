QMP Absolute quantification project: Quantitative Microbiome Profile (QMP)
=====================================

*Creator: Anne Boddeke* 
*Year: 2018* 

For steps on sample selection and data preparation etc., see: https://github.com/BoddekeAM/hello-world

NOTES:
- This is the QMP script of Vandeputte, modified for out data
- We don't do copy number correction
- Input are estimated no. of reads and FISH counts

1. Read table 
-----------------------------------------------
```
# Open input tables
estimated_reads_table <- read.table("FISH_s_e_reads.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
e_reads_1 <- as.data.frame(estimated_reads_table)
e_reads <- e_reads_1[-c(62, 63),] # Remove P37_M1 and P37_M4
rownames(e_reads)[64] <- "P39_M1"

FISH_counts <- read.table("RiseUp_Tot_Counts_M1M4.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
FISH_counts <- as.data.frame(FISH_counts)
FISH_counts$Tot_count <- as.numeric(gsub(",", ".", FISH_counts$Tot_count))
colnames(FISH_counts) <- "Bacteria_per_gram"
FISH_counts$ID <- rownames(FISH_counts)
FISH_counts <- as.data.frame(FISH_counts[-c(22, 34, 44, 124),]) # Remove P11_M4, P20_M4, P26_M4 and P71_M4
FISH_counts$ID=NULL

# Input files to use are 'e_reads' and 'FISH_counts'
```

2. The function
-----------------------------------------------
```
### TRY OUT: See below the function for notes on this script and what happened while executing it. ###

####################################################################################################################

# rarefaction to even sampling depth #
# author: Doris Vandeputte           #
######################################
# this script doesn't include copy number correction, a function for copy number correction is included in RDP classifier 2.12 
# this script uses function rarefy_even_depth from phyloseq 1.20.0, it needs package phyloseq to be installed and loaded in order to work.
# with cnv_corrected_abundance_table: a copy number variation corrected abundance table with sample-identifiers as rows, copy number corrected taxa-abundances as columns
# with cell_counts_table: a table with sample-identifiers as rows, cell counts as columns 

rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) 
{
  try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
  cell_counts_table = t(cell_counts_table[order(row.names(cnv_corrected_abundance_table)),]) # make sure the order of the samples is the same in both files  
  sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
  rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
  {
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(write.table(QMP, file = "./QMP_table.txt", quote = F, sep = "\t"))
  
}

rarefy_even_sampling_depth(e_reads, FISH_counts)

# Open the generated table
QMP_s_table <- read.table("QMP_table.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
```

# NOTES 
- 1.  First run: rownames didn't match. Turns out some data is missing and mislabelled. As follows:
         - M1 sample of P39 was mislabelled: in FISH I called it P39_M1; in e_reads it's called 39_M2
         - In the estimated reads data, missing are P26_M4 and P71_M4
         - In the FISH data, missing are P37_M1 and P37_M4
     Solution: rename P39_M2 to P39_M1 in e_reads data; remove missing samples in both datasets

- 2.  The function returns the QMP to you, but my input tables are too big so it can't print the whole thing
     Also I just want it to write a table for me.
     So I try to change the end of the function so that it will not return QMP but make and write a table

- 3.  Now the function ends with writing a table which it does not automatically open
     So after executing the function, read in the table.

- 4.  Do QMP again for the estimated reads tabel of all taxonomical levels
```
       #Open estimated reads table with all taxonomical levels
        # Read table with estimated reads at M1 and M4
        e_reads_all <- read.table("estimated_reads_metaphlan.txt", header = T, sep="\t", as.is=TRUE, row.names = 1)
        #### Correct and simplify IDs
        ### Remove "_metaphlan"
        colnames(e_reads_all) <- gsub("_metaphlan", "", colnames(e_reads_all))
        ### Rename p001_M1 IDs to P1_M1 IDs
        ## Substitute "p00" for "P0" and "p0" for "P" (as some have a name with "p0") 
        colnames(e_reads_all) <- gsub("p00", "P0", colnames(e_reads_all))
        colnames(e_reads_all) <- gsub("p0", "P", colnames(e_reads_all))
        #### Names are correct now ####
        # Improve format and remove samples with missing data
        e_reads_1 <- as.data.frame(t(e_reads_all))
        e_reads <- e_reads_1[-c(62, 63),] # Remove P37_M1 and P37_M4
        rownames(e_reads)[64] <- "P39_M1"
        ##### Final table is 'e_reads'; This table has all taxonomical levels, correct row/colnames and no missing data
        write.table(e_reads, file = "./e_reads_all_levels_corrected.txt", quote = F, sep = "\t")
      # Perform function
      rarefy_even_sampling_depth(e_reads, FISH_counts)
```



