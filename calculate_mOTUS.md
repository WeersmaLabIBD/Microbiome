#!/bin/bash

for i in p*; do 
sid=$i

## Parameter for cluster job, modify if is needed. 
        echo "#!/bin/bash"  > "$sid".sh
        echo "#SBATCH --job-name="$sid"_metagenomes" >> "$sid".sh
        echo "#SBATCH --error="$sid".err" >> "$sid".sh
        echo "#SBATCH --output="$sid".out" >> "$sid".sh
        echo "#SBATCH --mem=20gb" >> "$sid".sh
        echo "#SBATCH --time=5:59:00" >> "$sid".sh
        echo "#SBATCH --cpus-per-task=6" >> "$sid".sh

## Add variables to the path and load modules
        echo "module load picard " >> "$sid".sh
        echo "module load Python" >> "$sid".sh
        echo "module load Bowtie2" >> "$sid".sh
        #echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/" >> "$sid".sh
        #echo "export PATH=\$PATH:/home/umcg-avich/.local/bin" >> "$sid".sh
        echo "module load kneaddata" >> "$sid".sh
        #echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/bin/">> "$sid".sh


## QC in old way now integrated in BioBakery (a.k.a kneaddata)

        echo "kneaddata --input ./"$sid"/"$sid"_1.fastq -t 6 -p 7 --input ./"$sid"/"$sid"_2.fastq -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./"$sid"/ --log ./"$sid"/out.log">> "$sid".sh

# Join pair-end reads file -> Metaphlan and Humann2 don't use merged pair reads 

        echo "cat ./"$sid"/"$sid"/"$sid"_1_kneaddata_paired_1.fastq > ./"$sid"/"$sid"_kneaddata_merged.fastq " >> "$sid".sh
        echo "cat ./"$sid"/"$sid"/"$sid"_1_kneaddata_paired_2.fastq

# Use mOTUS_v2 for taxonomical classification

        echo "ml SAMtools/1.5-foss-2015b" >> "$sid".sh
        echo "ml BWA/0.7.15-foss-2015b" >> "$sid".sh
        echo "ml Python" >> "$sid".sh
        echo "mkdir ./"$sid"/mOTUs"
        echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/mOTUS/mOTUs_v2" >> "$sid".sh
#Counts
        echo "motus profile -f ./"$sid"/"$sid"_1_kneaddata_paired_1.fastq -r ./"$sid"/"$sid"_1_kneaddata_paired_2.fastq -o ./"$sid"/"$sid"_counts_profile.txt -q -t 6 -c -g 3 -y insert.scaled_counts" >> "$sid".sh
#Abundances
        echo "motus profile -f ./"$sid"/"$sid"_1_kneaddata_paired_1.fastq -r ./"$sid"/"$sid"_1_kneaddata_paired_2.fastq -o ./"$sid"/"$sid"_abundance_profile.txt -q -t 6 -g 3 -y insert.scaled_counts" >> "$sid".sh

done
