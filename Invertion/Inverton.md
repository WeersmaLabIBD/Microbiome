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
ml BioPython
ml numpy/1.11.0-foss-2015b-Python-2.7.11
PATH=$PATH:~/bowtie-1.2.2-linux-x86_64/
ml SAMtools
PATH=$PATH:~/bin/
```

