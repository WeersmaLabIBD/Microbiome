Description 
====================================

*2_steps_logistic_regression.R*  : Logistic regression on presence / absence using feature selection

*DUDes_2_mpa.py* : convert DUDes (metagenomic tool) output to mpa (Metaphlan like) output

*Filtering_metadata.R* : filter a phenotypic data by min sequencing depth, number of NAs allowed, columns with one unique value

*Filter_taxonomy.R* : filter taxonomy table by prevalence and minimum relative abundance and/or a list of id's 

*PCoA_plots.R* : script to create principal coordinate analysis plots

*PCoA_with_centroids.R* : add centroids to the PCoA plot

*Shannon_index_calculation.R* : calculates richness ( or other eveness/richness measurement) and creates violin plots

*add_missing_taxa_levels.py* : based on mpa-like taxonomies, adds the missing taxonomical levels. 

*annotate_metacyc.py* : connects to MetaCYC database and download relevant information of microbial pathways

*bamTofastq.sh* : file formating conversion bam -> fastq (SLURM cluster)

*cluster_dendogram.R* : creates a dendogram and colors base on user pre-defined categories

*composition_barplots.R* : barplots of microbial composition

*convert_Braken_to_relative_abundances.py* : from Braken output to relative abundances

*create_samplesheet.sh* : samplesheet for qtl pipeline

*create_scripts_metagenomics.sh* : create SLURM scripts for our metagenomic pipeline (2017)

*extract_info_logs_Maaslin.sh* : extracts standard errors from Maaslin log files. 

*humann2_config_local.sh* : HumaNn2 local configuration

*lm_permutation.R*: randomly subset a table and perform linear models on it (several permutations)

*normalize_data.R* : normalization steps for microbiome data

*parse_fastQC.R* : extract info from FASTQC output

*parse_log_files_metagenomic_pipe.sh* : extract summary of QC steps in our metagenomic pipeline (2017)

*pricess_reads_count.sh* : same as parse_log_files_metagenomic_pipe.sh

*select_last_taxa_level.py*: extract the lowest taxonomical level. 

*transform_taxaid_to_mpa_taxonomy.py* : transform taxaID from NCBI and transform it to mpa-like taxonomy names

