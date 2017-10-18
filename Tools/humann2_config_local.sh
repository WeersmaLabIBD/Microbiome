 tar -xzvf humann2-0.10.0.tar.gz
 module load Python
cd humann2-0.10.0
python setup.py install --user
export PATH=$PATH:/home/umcg-$user/.local/bin 
humann2_config --update database_folders nucleotide /groups/umcg-gastrocol/tmp04/metagenomic_tools/chocophlan
humann2_config --update database_folders protein /groups/umcg-gastrocol/tmp04/metagenomic_tools/uniref
