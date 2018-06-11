Collaboration Prof.Hector Rodriguez: Phenolic metabolites and Crohns
=====================================================================

Creator: Arnau Vich 
Year:2018
Description: Blastx aligment to set of AA sequences

```Step 1: Create database
---------------------------


export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/ncbi-blast-2.6.0+/bin/

./ncbi-blast-2.6.0+/bin/makeblastdb -in sequences.txt -parse_seqids -dbtype prot
```

```

Step 2: Create aligment script
-------------------------------


#!/bin/bash

for i in *bam; do 
        sid=${i%.bam}

        echo "#!/bin/bash" > "$sid".sh
        echo "#SBATCH --job-name="$sid"_metagenomes" >> "$sid".sh
        echo "#SBATCH --error=$sid.err" >> "$sid".sh
        echo "#SBATCH --output=$sid.out" >> "$sid".sh
        echo "#SBATCH --mem=35gb" >> "$sid".sh
        echo "#SBATCH --time=05:50:00" >> "$sid".sh
        echo "#SBATCH --cpus-per-task=5" >> "$sid".sh

        echo "ml picard" >> "$sid".sh
        echo "ml kneaddata" >> "$sid".sh

        echo "mkdir ./$sid/"  >> "$sid".sh
        echo "mkdir ./$sid/clean_reads/"  >> "$sid".sh
        echo "mkdir ./$sid/alignment/"  >> "$sid".sh
        echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/ncbi-blast-2.6.0+/bin/" >> "$sid".sh

        echo "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$sid.bam F=./$sid/clean_reads/"$sid".fastq1 F2=./$sid/clean_reads/"$sid".fastq2" >> "$sid".sh
        echo "kneaddata --input ./$sid/clean_reads/"$sid".fastq1 -t 6 -p 7 --input ./$sid/clean_reads/"$sid".fastq2 -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$sid/clean_reads/ --log ./$sid/clean_reads/$sid.log" >> "$sid".sh
        echo "cat ./"$sid"/clean_reads/*kneaddata_paired_1.fastq >> ./"$sid"/"$sid"_merged.fq " >> "$sid".sh
        echo "cat ./"$sid"/clean_reads/*kneaddata_paired_2.fastq >> ./"$sid"/"$sid"_merged.fq " >> "$sid".sh
        echo "less ./"$sid"/"$sid"_merged.fq | awk '{if(NR%4==1) {printf(\">%s\n\",substr(\$0,2));} else if(NR%4==2) print;}' > ./"$sid"/"$sid"_merged.fa">> "$sid".sh
        #echo "rm *.fq">> "$sid".sh
        echo "blastx -query ./"$sid"/"$sid"_merged.fa -db sequences.txt -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp\" -num_threads 5  -evalue 1e-5 -qcov_hsp_perc 50 > ./$sid/alignment/"$sid".txt" >> "$sid".sh
        echo "rm -r ./$sid/clean_reads/">> "$sid".sh
        echo "rm ./$sid/*.fa">> "$sid".sh
        echo "rm ./$sid/*.fq">> "$sid".sh
done
```

```

3.Example script
-------------------


#!/bin/bash
#SBATCH --job-name=IBDFEC0729_metagenomes
#SBATCH --error=IBDFEC0729.err
#SBATCH --output=IBDFEC0729.out
#SBATCH --mem=35gb
#SBATCH --time=05:50:00
#SBATCH --cpus-per-task=5
ml picard
ml kneaddata
mkdir ./IBDFEC0729/
mkdir ./IBDFEC0729/clean_reads/
mkdir ./IBDFEC0729/alignment/
export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/ncbi-blast-2.6.0+/bin/
java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=IBDFEC0729.bam F=./IBDFEC0729/clean_reads/IBDFEC0729.fastq1 F2=./IBDFEC0729/clean_reads/IBDFEC0729.fastq2
kneaddata --input ./IBDFEC0729/clean_reads/IBDFEC0729.fastq1 -t 6 -p 7 --input ./IBDFEC0729/clean_reads/IBDFEC0729.fastq2 -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./IBDFEC0729/clean_reads/ --log ./IBDFEC0729/clean_reads/IBDFEC0729.log
cat ./IBDFEC0729/clean_reads/*kneaddata_paired_1.fastq >> ./IBDFEC0729/IBDFEC0729_merged.fq 
cat ./IBDFEC0729/clean_reads/*kneaddata_paired_2.fastq >> ./IBDFEC0729/IBDFEC0729_merged.fq 
less ./IBDFEC0729/IBDFEC0729_merged.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ./IBDFEC0729/IBDFEC0729_merged.fa
blastx -query ./IBDFEC0729/IBDFEC0729_merged.fa -db sequences.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp" -num_threads 5  -evalue 1e-5 -qcov_hsp_perc 50 > ./IBDFEC0729/alignment/IBDFEC0729.txt
rm -r ./IBDFEC0729/clean_reads/
rm ./IBDFEC0729/*.fa
rm ./IBDFEC0729/*.fq

```
