out="/scratch/p282673/Inverton/PhaseFinder/IBD/"
path_fa="/scratch/p282673/Inverton/assembly_IBD_done/"
path_fq="/scratch/umcg-rgacesa/binning_LLD_IBS_IBD/fastqs_all"

for sample in $path_fa/c_*/*.fasta

do

name=$(echo ${sample%.fasta})
i=$(basename $name)
n=$(echo ${i#*_})


echo -e "#!/bin/bash" > $i.job.sh
echo -e "#SBATCH --job-name=$i" >> $i.job.sh
echo -e "#SBATCH --error=$i.err" >> $i.job.sh
echo -e "#SBATCH --output=$i.out" >> $i.job.sh
echo -e "#SBATCH --mem=64gb" >> $i.job.sh
echo -e "#SBATCH --time=5:00:00" >> $i.job.sh
echo -e "#SBATCH --cpus-per-task=6" >> $i.job.sh


echo -e "ml Python/2.7.11-foss-2016a" >> $i.job.sh
echo -e "ml Biopython/1.65-foss-2016a-Python-2.7.11" >> $i.job.sh
echo -e "ml BEDTools" >> $i.job.sh
echo -e "ml SAMtools" >> $i.job.sh
echo -e "ml Bowtie" >> $i.job.sh

echo -e "PATH=\$PATH:/home/p282673/Inverton/software/bin/" >> $i.job.sh
echo -e "PATH=\$PATH:/home/p282673/Inverton/software/EMBOSS/EMBOSS-6.6.0/emboss/" >> $i.job.sh

echo -e "python ~/Inverton/software/PhaseFinder-master/PhaseFinder.py locate -f $sample -t $i.einverted.tab -g 15 85 -p" >> $i.job.sh
echo -e "echo -e "=====test.einverted.tab done====="" >> $i.job.sh

echo -e "python ~/Inverton/software/PhaseFinder-master/PhaseFinder.py create -f $sample -t $i.einverted.tab -s 1000 -i $i.ID.fasta" >> $i.job.sh

echo -e "echo -e "=====mimic inverted repeats done====="" >> $i.job.sh

echo -e "ml Python/2.7.14-fosscuda-2018a" >> $i.job.sh

echo "python ~/Inverton/software/PhaseFinder-master/PhaseFinder.py ratio -i $i.ID.fasta -1 $path_fq/$n\_1.fastq -2 $path_fq/$n\_2.fastq  -p 16 -o $i.out" >> $i.job.sh

echo -e "echo -e "=====calculate span and pe reads done====="" >> $i.job.sh

echo -e "Create $i is finifshed"

done
                     
