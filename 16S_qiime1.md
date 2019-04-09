---
title: "MiBioGen 16s pipeline"
author: "Alexander Kurilshikov"
update: "Shixian"
date: "April 09, 2019"
output: html_document
---

# Overall description

This is the pipeline for analysis of 16s based on qiime and silva database.

## 1. software install
*1.1 qiime installation*
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod a+x Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
source activate qiime1
print_qiime_config.py
```
*1.2 SILVA installation*
```
#download SILVA v.119 database
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_119_release.zip
unzip Silva_119_release.zip
#create custom config file
echo "\
pick_otus_reference_seqs_fp $PWD/Silva119_release/rep_set/97/Silva_119_rep_set97.fna
pynast_template_alignment_fp $PWD/Silva119_release/core_alignment/core_Silva119_alignment.fna
assign_taxonomy_reference_seqs_fp $PWD/Silva119_release/rep_set/97/Silva_119_rep_set97.fna
assign_taxonomy_id_to_taxonomy_fp $PWD/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt" > qiime_config
```
## 2. join fastq files and remove barcodes

```
ml numpy

sed 's/2:N:0/1:N:0/g' Barcode.fastq > Barcode.fixed.fq

join_paired_ends.py -f $PWD/forward_reads.fastq -r $PWD/reverse_reads.fastq -b $PWD/barcodes.fastq -o $PWD/fastq-join_joined
(output files will be used in split_library.py, make sure barcode.fq and reads.fq have the same order)

split_libraries_fastq.py -i fastqjoin.join.fastq -b fastqjoin.join_barcodes.fastq -o slout_q20/ -m map.txt --store_qual_scores -q 19 
(if error appears, add one parameter mentioned in err file; output files will be used in SILVA mapping)

```

## 3. SILVA mapping

```
PROCNUM=1
source activate qiime1
export QIIME_CONFIG_FP=$PWD/qiime_config
pick_closed_reference_otus.py -i FASTA_FILE -o RESULT_FOLDER -a -O $PROCNUM
cd RESULT_FOLDER
biom convert -i otu_table.biom -o otu_table.tsv --to-tsv --header-key taxonomy
cat otu_table.tsv|tail -n+2 |perl -pe "s/#OTU ID/OTU_ID/" > temp.tsv
mv temp.tsv otu_table.tsv
```

## 3. Getting taxonomies from OTU table

```
ml R

Rscript Taxa.R
```
