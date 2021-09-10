---
title: "16s pipeline"
author: "Alexander Kurilshikov"
update: "Shixian"
date: "April 09, 2019"
---

# Overall description

This is the pipeline for analysis of 16s based on qiime and silva database OTU-picking.

## 1. software install
*1.1 qiime installation*
```
# wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# chmod a+x Miniconda3-latest-Linux-x86_64.sh
# bash ./Miniconda3-latest-Linux-x86_64.sh

Note: Miniconda3 might have conflicts with qiime1, so try miniconda2:
https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

install https://pypi.org/project/matplotlib/1.4.3/; ml Python/2.7xxx; ml numpy; python setup.py build; python setup.py install
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

source activate qiime1
Or source activate /home/umcg-hushixian/miniconda2/envs/qiime1
ml numpy
print_qiime_config.py
```
*1.2 SILVA installation*
```
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_119_release.zip
unzip Silva_119_release.zip

echo "\
pick_otus_reference_seqs_fp $PWD/Silva119_release/rep_set/97/Silva_119_rep_set97.fna
pynast_template_alignment_fp $PWD/Silva119_release/core_alignment/core_Silva119_alignment.fna
assign_taxonomy_reference_seqs_fp $PWD/Silva119_release/rep_set/97/Silva_119_rep_set97.fna
assign_taxonomy_id_to_taxonomy_fp $PWD/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt" > qiime_config
```
*1.3 ea-utils installation*
```
download from https://expressionanalysis.github.io/ea-utils/
make && make install
export PATH=$PATH:ABSOLUTE PATH TO fastq-join
```
## 2. join fastq files and remove barcodes

```
ml numpy

sed 's/2:N:0/1:N:0/g' Barcode.fastq > Barcode.fixed.fq

join_paired_ends.py -f forward_reads.fastq -r reverse_reads.fastq -b Barcode.fixed.fq -o fastq-join_joined
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
## 4. Rarefaction
```
Rscript step0.2_run_rarefaction.R SEQUENCES.FASTA 40000
source activate qiime1
filter_fasta.py -f SEQUENCES.FASTA -o SEQ_RARIFIED.FASTA -s 2filter.ids
rm 2filter.ids
```

## 5. taxonomies from OTU table

```
ml R

Rscript Taxa.R
```

## 6. Assign to KO and metaboic pathways 
##    (use Tax4Fun for SILVA based)

*6.1 software installation*

```
# download Tax4FUn and SILVA annotation database. NOTE, here use SILVA 123 as an example
wget http://tax4fun.gobics.de/Tax4Fun/ReferenceData/SILVA123.zip
wget http://tax4fun.gobics.de/SilvaSSURef_123_NR.zip
wget http://tax4fun.gobics.de/Tax4Fun/Tax4Fun_0.3.1.tar.gz

unzip SilvaSSURef_123_NR.zip
unzip SILVA123.zip

ml R
R

install.packages("qiimer")
install.packages("Matrix")
install.packages("biom")
install.packages("Tax4Fun_0.3.1.tar.gz", repos = NULL, type = "source")

```

*6.2 check OTU table*

Tax4Fun has a really strict requirement of the input format, so please double-check.

Otherwise thousands of errors will happen. 

For example: "length of 'dimnames' [1] not equal to array extent" 

```
# Standard format of otu_table.txt, take care with taxnomy column, NO "k__, p__ ", etc !!!

# Constructed from biom file                            
#OTU ID r3.7    r7.67   r10.21  taxonomy
denovo0 1.0     0.0     0.0     Bacteria; Firmicutes; Clostridia; Clostridiales; Lachnospiraceae; ; 
denovo1 0.0     1.0     0.0     Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; ; 
denovo2 0.0     0.0     1.0     Bacteria; Firmicutes; Bacilli; Lactobacillales; Streptococcaceae; Streptococcus
denovo3 0.0     0.0     0.0     Bacteria; Firmicutes; Clostridia; Clostridiales; Lachnospiraceae; Lachnospira;
```

*6.3 get KO terms (relative abundance)*

```
R

library(Tax4Fun)

QIIMESingleData <- importQIIMEData("otu_table.txt")
Tax4FunOutput <- Tax4Fun(QIIMESingleData, "SILVA123", fctProfiling = TRUE, refProfile = "UProC", shortReadMode = TRUE, normCopyNo = TRUE)
KO_table = t(Tax4FunOutput$Tax4FunProfile)

write.table("ID\t", file="KO_table.txt",append = FALSE, quote = FALSE, sep="\t",eol = "", na = "NA", dec = ".", row.names = F,col.names = F)
write.table(KO_table, file="KO_table.txt",append = T, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
```

*6.4 get only metabolic terms (relative abundance)*

```
Tax4FunOutput <- Tax4Fun(QIIMESingleData, "SILVA123", fctProfiling = FALSE, refProfile = "UProC", shortReadMode = TRUE, normCopyNo = TRUE)
KO_table = t(Tax4FunOutput$Tax4FunProfile)

write.table("ID\t", file="KO_table_fct.txt",append = FALSE, quote = FALSE, sep="\t",eol = "", na = "NA", dec = ".", row.names = F,col.names = F)
write.table(KO_table, file="KO_table_fct.txt",append = T, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
```

## 7. Assign to metaboic pathways 
##    (use PICRUST for GreenGene based)

*7.1 software installation*

```
conda create -n picrust1 -c bioconda -c conda-forge picrust
conda activate picrust1
ml numpy
```

*7.2 normalize data*

```
normalize_by_copy_number.py --gg_version 13_5 -i $table/$line.taxonomy/otu_table.biom -o $out/$line.normalized.biom
```

*7.3 metabolic predictions*

```
predict_metagenomes.py -i $out/$line.normalized.biom -o $out/$line.predictions.biom -a $out/$line.nsti.tab
```

*7.4 get KEGG terms*

```
categorize_by_function.py -f -i $out/$line.predictions.biom -c KEGG_Pathways -l 3 -o $out/$line.KEGG.txt
```
