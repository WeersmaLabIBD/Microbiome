#!/bin/bash

#Creator: Arnau Vich
#Year: 2017
#USAGE: bash create_scripts_metagenomics.sh
#Description: For each bam file in a directory creates a two steps script for metagenomic analyses in SLURM enviorament. 
#WARNING: Adapt parameters / variables before running the script!! 

for i in *bam; do 
sid=${i%.bam}

## Parameter for cluster job, modify if is needed. 
        echo "#!/bin/bash"  > "$sid".sh
        echo "#SBATCH --job-name="$sid"_metagenomes" >> "$sid".sh
        echo "#SBATCH --error="$sid".err" >> "$sid".sh
        echo "#SBATCH --output="$sid".out" >> "$sid".sh
        echo "#SBATCH --mem=70gb" >> "$sid".sh
        echo "#SBATCH --time=5:00:00" >> "$sid".sh
        echo "#SBATCH --cpus-per-task=6" >> "$sid".sh

## Add variables to the path and load modules
        echo "module load picard " >> "$sid".sh
        echo "module load Python" >> "$sid".sh
        echo "module load Bowtie2" >> "$sid".sh
        echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/" >> "$sid".sh
        ##IMPORTANT, set your own home directory in the cluster!!
	echo "export PATH=\$PATH:/home/umcg-avich/.local/bin"
        #echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/kneaddata/" >> "$sid".sh
        #echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/bin/">> "$sid".sh

## Create directories 

        echo "mkdir ./"$sid"" >> "$sid".sh
        echo "mkdir ./"$sid"/filtering_data/" >> "$sid".sh
        echo "mkdir ./"$sid"/clean_reads/" >> "$sid".sh
        echo "mkdir ./"$sid"/metaphlan/" >> "$sid".sh
        echo "mkdir ./"$sid"/humann2/" >> "$sid".sh
        echo "mkdir ./"$sid"/DUDes/" >> "$sid".sh

## BAM to FASTQ

        echo "java -jar "\${EBROOTPICARD}"/picard.jar SamToFastq I="$i" F=./"$sid"/filtering_data/"$sid".fastq1 F2=./"$sid"/filtering_data/"$sid".fastq2" >> "$sid".sh
        echo "echo "Starting reads cleaning"" >> "$sid".sh
## QC in old way now integrated in BioBakery (a.k.a kneaddata)

        echo "kneaddata --input ./"$sid"/filtering_data/"$sid".fastq1 -t 6 -p 7 --input ./"$sid"/filtering_data/"$sid".fastq2 -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./"$sid"/filtering_data/ --log ./"$sid"/clean_reads/"$sid".log">> "$sid".sh

# Join pair-end reads file -> Metaphlan and Humann2 don't use merged pair reads 

        echo "cat ./"$sid"/filtering_data/"$sid"_kneaddata_paired_1.fastq > ./"$sid"/filtering_data/"$sid"_kneaddata_merged.fastq " >> "$sid".sh
        echo "cat ./"$sid"/filtering_data/"$sid"_kneaddata_paired_2.fastq >> ./"$sid"/filtering_data/"$sid"_kneaddata_merged.fastq " >> "$sid".sh
        echo "mv ./"$sid"/filtering_data/*kneaddata_paired_1.fastq ./"$sid"/clean_reads/" >> "$sid".sh
        echo "mv ./"$sid"/filtering_data/*kneaddata_paired_2.fastq ./"$sid"/clean_reads/ " >> "$sid".sh
        echo "mv ./"$sid"/filtering_data/*kneaddata_merged.fastq ./"$sid"/clean_reads/">> "$sid".sh
        echo "rm -r ./"$sid"/filtering_data/" >> "$sid".sh
## Run metaphlan (check defaults for extra options). Align info in G96213_kneaddata_merge.fastq.bowtie2out.txt 
       echo "echo "Starting taxonomy classification using Metaphlan"" >> "$sid".sh
       echo "metaphlan2.py ./"$sid"/clean_reads/"$sid"_kneaddata_merged.fastq  --input_type multifastq --mpa_pkl /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200.pkl --nproc 6 -o ./"$sid"/metaphlan/"$sid"_metaphlan.txt --tmp_dir ./"$sid"/clean_reads/" >> "$sid".sh

## Run DUDes

        echo "echo "Starting taxonomy profiling using DUDes"" >> "$sid".sh
        echo "bowtie2 -x /groups/umcg-gastrocol/tmp03/metagenomic_tools/dudes_v0_07/custom_db/db_refseq_20052017 --no-unal --fast -p 6 -k 50 -1 ./"$sid"/clean_reads/"$sid"_kneaddata_paired_1.fastq -2 ./"$sid"/clean_reads/"$sid"_kneaddata_paired_2.fastq -S ./"$sid"/DUDes/"$sid"_output.sam " >> "$sid".sh
        echo "python3 /groups/umcg-gastrocol/tmp03/metagenomic_tools/dudes_v0_07/DUDes.py -s ./"$sid"/DUDes/"$sid"_output.sam -d /groups/umcg-gastrocol/tmp03/metagenomic_tools/dudes_v0_07/custom_db/DUDES_refseq_db.npz -t 6 -m 50 -a 0.0005 -l strain -o ./"$sid"/DUDes/"$sid" " >> "$sid".sh
        echo "sbatch "$sid"_step2.sh" >> "$sid".sh

## Parameter for cluster job, modify if is needed. 
        echo "#!/bin/bash"  > "$sid"_step2.sh
        echo "#SBATCH --job-name="$sid"_metagenomes" >> "$sid"_step2.sh
        echo "#SBATCH --error="$sid".err" >> "$sid"_step2.sh
        echo "#SBATCH --output="$sid".out" >> "$sid"_step2.sh
        echo "#SBATCH --mem=70gb" >> "$sid"_step2.sh
        echo "#SBATCH --time=5:00:00" >> "$sid"_step2.sh
        echo "#SBATCH --cpus-per-task=6" >> "$sid"_step2.sh

## Add variables to the path and load modules
        echo "module load picard " >> "$sid"_step2.sh
        echo "module load Python" >> "$sid"_step2.sh
        echo "module load Bowtie2" >> "$sid"_step2.sh
        echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/" >> "$sid"_step2.sh
        echo "export PATH=\$PATH:/home/umcg-avich/.local/bin" >> "$sid"_step2.sh



## Run humann2
        echo "echo "Starting pathways prediction using Humann2"" >> "$sid"_step2.sh
        echo "humann2 --input ./"$sid"/clean_reads/"$sid"_kneaddata_merged.fastq --output ./"$sid"/humann2/ --taxonomic-profile ./"$sid"/metaphlan/"$sid"_metaphlan.txt --threads 6 --o-log ./"$sid"/clean_reads/"$sid".full.humann2.log --remove-temp-output --protein-database /groups/umcg-gastrocol/tmp03/metagenomic_tools/humann2-0.10.0/db/uniref/ --remove-temp-output" >> "$sid"_step2.sh
        # Use Uniref_50
        #echo "humann2 --input ./"$sid"/clean_reads/"$sid"_kneaddata_merged.fastq --output ./"$sid"/humann2/uniref_50 --taxonomic-profile ./"$sid"/metaphlan/"$sid"_metaphlan.txt --threads 6 --o-log ./"$sid"/clean_reads/"$sid".full.humann2.log --remove-temp-output --protein-database /groups/umcg-gastrocol/tmp03/metagenomic_tools/humann2-0.10.0/db/uniref_50/" >> "$sid".sh

## Remove cleaned reads. Commented to minimize the space/memory requeriments. Uncomment to keep cleaned reads
	    echo "mv ./"$sid"/clean_reads/*.log ./"$sid"/" >> "$sid"_step2.sh
        echo "rm ./"$sid"/DUDes/"$sid"_output.sam" >> "$sid"_step2.sh
        echo "rm -r ./"$sid"/clean_reads/" >> "$sid"_step2.sh

done
