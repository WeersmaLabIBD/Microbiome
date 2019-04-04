# Inverton identification from contigs 

```
Installing EMBOSS cost me 6 hours!!! Note, when try to configure software, make sure there is no space in the name of software

mkdir EMBOSS
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
tar zxfv EMBOSS-6.6.0.tar.gz
rm EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0
./configure --prefix=/home/hou/Software/bin
make
make install
```

```
ml Biopython
ml numpy/1.11.0-foss-2015b-Python-2.7.11
PATH=$PATH:~/bowtie-1.2.2-linux-x86_64/
ml SAMtools
ml BEDTools
PATH=$PATH:~/bin/
PATH=$PATH:/home/umcg-hushixian/EMBOSS-6.6.0/emboss/ 
ml pandas/0.16.2-foss-2015b-Python-2.7.9

out="/groups/umcg-weersma/tmp03/Inverton/test/"

for sample in /groups/umcg-gastrocol/prm03/rawdata/IBD_MGS/*.bam

name=$(echo ${sample%.bam})
i=$(basename $name)
mkdir $out/output/$i

java -jar /apps/software/picard/2.18.26-Java-1.8.0_74/picard.jar SamToFastq I=$sample F=$out/output/$i/$i.1.fastq.gz F2=$out/output/$i/$i.2.fastq.gz

python /groups/umcg-weersma/tmp03/Inverton/software/PhaseFinder.py locate -f $out/test.fa -t $out/test.einverted.tab -g 15 85 -p 

python /groups/umcg-weersma/tmp03/Inverton/software/PhaseFinder.py create -f $out/test.fa -t $out/test.einverted.tab -s 1000 -i $out/test.ID.fasta

python /groups/umcg-weersma/tmp03/Inverton/software/PhaseFinder.py ratio -i $out/test.ID.fa -1 $out/$i/$i.1.fastq.gz -2 $out/$i/$i.2.fastq.gz -p 16 -o $out/output/$i/

```

```

python /groups/umcg-weersma/tmp03/Inverton/software/PhaseFinder.py 

```
