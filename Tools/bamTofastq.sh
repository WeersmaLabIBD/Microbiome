#!/bin/bash

for sample in *.bam

do 
  i=$(echo ${sample%.bam})

  echo -e "#!/bin/bash" > $i.bamTofq.sh
  echo -e "#SBATCH --job-name=$i.bamtofastq" >> $i.bamTofq.sh
  echo -e "#SBATCH --error=$i.bamtofastq.err" >> $i.bamTofq.sh
  echo -e "#SBATCH --output=$i.bamtofastq.out" >> $i.bamTofq.sh
  echo -e "#SBATCH --mem=30gb" >> $i.bamTofq.sh
  echo -e "#SBATCH --time=9:59:00" >> $i.bamTofq.sh
  echo -e "#SBATCH --cpus-per-task=6" >> $i.bamTofq.sh

  echo -e "module load picard"  >> $i.bamTofq.sh
  echo -e "module load Python/3.4.1-foss-2015b" >> $i.bamTofq.sh

  echo -e "export PATH=\$PATH:~/.local/bin" >> $i.bamTofq.sh
  echo -e "outdir=\"/groups/umcg-weersma/tmp03/WES-project/rawdata/ngs\"" >> $i.bamTofq.sh

  echo -e "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$i.bam F=\${outdir}/${i}_1.fastq.gz F2=\${outdir}/${i}_2.fastq.gz" >> $i.bamTofq.sh
done
