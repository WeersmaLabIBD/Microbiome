Growth rate calculation (Korem et al.)
=========================================

Creator: Arnau Vich

Year: 2015

Old cluster configuration
-----------------------------

```
#!/bin/bash

#load packages and Python
module load Python/2.7.9-foss-2015b
module load numpy/1.9.2-foss-2015b-Python-2.7.9 
module load pandas/0.16.2-foss-2015b-Python-2.7.9 
module load lmfit/0.8.3-foss-2015b-Python-2.7.9
module load dill/0.2.4-foss-2015b-Python-2.7.9

#Export paths 
export PATH=$PATH:/groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/GEM/
export PATH=$PATH:/groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/

#Array containing the name of folders to avoid 
declare -a avoid=(G82573 G82658 G82675 G82704 G82808 G83261 G83429 G83901 G84092)
```

Run commands
-----------

````
#Initiate names for some arrays that we will use later 
all_folders=()
keep=()
delete=()

#Path to input files and the files already processed 
dir=/groups/umcg-tifn/prm02/projects/LLD_MGS_results/cleaned_reads/
dir_done=/groups/umcg-weersma/prm02/raw_data/


# Check all folders in a directory with the raw data (input)
for d in "$dir"*;do
  	#folders=$(basename "$d")
    all_folders+=("$d")
done

# Add the path of the folders already analysed 
for del in ${avoid[@]}
do
	delete+=("$dir""$del")
	#keep=(${all_folders[@]/$delete})
	#tapes=(${tapes[@]//*$i*})
done

# Set the input files in array 
keep=( $(printf "%s\n" "${all_folders[@]}" "${delete[@]}" | sort | uniq -u) )


# Loop for each of the input files 
for d in "${keep[@]}";do
  	s="$d"/*.gz
    name=$(basename $s .tar.gz)
    check="$name".ptr
    		#echo "$check"
    # If we already process the file skip the analysis (it's present in the output folder )
    if [ -e /groups/umcg-weersma/prm02/rawdata/coverage_files_lld/"$check" ] || [-e /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/coverage_files_lld_2/"$check" ]; then
        echo "The file "$check" was already done" >> finish_test.txt
    else
    	my_path="$d"/"$check"
        #echo "I'm going to copy... "$my_path" " >> happy.txt 
        #Copy from the permanent storage to our working directory 
        cp "$d"/*.gz /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/
        #Uncompress and avoid multiple folders when uncompressing 
		tar -zxvf /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.tar.gz --strip-components 2 -C /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/
		#Check if uncompression finish correctly (samples that failed the file is corrupted/incomplete)
		if [ -e /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.um.1.fastq ] && [ -e /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.um.2.fastq ]; then
			s=groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*trim.60.single.um.fastq
		#name=$(basename $s .trim.60.single.um.fastq)
		# Map reads using GEM mapper 
			gem-mapper -I /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/index.gem -1 /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.um.1.fastq -2 /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.um.2.fastq -o /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/gem_mapping_files/"$name" -q ignore -d all -p --max-extendable-matches all --max-extensions-per-match 5 -v --threads 10
		# Remove temporal files
			rm /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.fastq
			rm /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.gz
			input_ptrc=/groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/gem_mapping_files/*.map
		# Run the coverage anaylisis and normalization 
			python PTRC.py -m $input_ptrc -outfol /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/coverage_files_lld_2/ -pe CA 
		# Remove temporal files 
			rm /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/gem_mapping_files/*.map
			rm /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/*.kill
			echo "File "$check" finish!!! " >> samples_complete_test.txt
		else
			echo "File "$check"  failed !!! " >> samples_failed.txt
			rm /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/*.kill
			rm /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.fastq
			rm /groups/umcg-weersma/scr01/Growth_rate/PTRC1.1/temporal/*.gz
		fi
	fi
done

## Ask data manager to move to prm02 storage
## Once is finish calculate the peak-to-trough ratio.
python PTRC.py -infol /groups/umcg-weersma/prm02/rawdata/coverage_files/coverage_files/ -preds generate -opreds /groups/umcg-gastrocol/tmp04/Growth_rate_IBD_LLD/IBD_LLD.preds -o /groups/umcg-gastrocol/tmp04/Growth_rate_IBD_LLD/IBD_LLD.csv -csv_output PTR
```
