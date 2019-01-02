Identify antibiotic resistance genes using ShortBred
========================================================

Creator: Arnau Vich 
Year: 2018 
Dependencies: ShortBred => https://bitbucket.org/biobakery/shortbred/wiki/Home

Create script to run in SLURM enviorament
-------------------------------------------

```{bash}
#!/bin/bash

for i in G*; do 
        sid=$i

        echo "#!/bin/bash" > "$sid".sh
        echo "#SBATCH --job-name="$sid"_metagenomes" >> "$sid".sh
        echo "#SBATCH --error=$sid.err" >> "$sid".sh
        echo "#SBATCH --output=$sid.out" >> "$sid".sh
        echo "#SBATCH --mem=25gb" >> "$sid".sh
        echo "#SBATCH --time=05:50:00" >> "$sid".sh
        echo "#SBATCH --cpus-per-task=6" >> "$sid".sh

        echo "module load Biopython/1.65-foss-2015b-Python-2.7.11" >> "$sid".sh

        echo "python ./ShortBred/shortbred_quantify.py --markers ./ShortBred/ShortBRED_CARD_2017_markers.faa --wgs ./data/"$i"/*.tar.gz --results ./data/"$i"/"$i"_out --tmp ./data/"$i"/"$i"_tmp --usearch ./usearch6.1.544 --threads 6" >> "$sid".sh
        echo "rm -r ./data/"$i"/*.tar.gz">> "$sid".sh
done

```
