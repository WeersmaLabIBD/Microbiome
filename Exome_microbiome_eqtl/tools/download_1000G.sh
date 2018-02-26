#!/bin/bash
#SBATCH --job-name=download_1000G
#SBATCH --error=download_1000G.err
#SBATCH --output=download_1000G.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load BCFtools

wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/chr{1-4,5-9,10-15,16-22X}.1kg.phase3.v5a.vcf.zip
for chrs in 1-4 5-9 10-15 16-22X; do unzip chr$chrs.1kg.phase3.v5a.vcf.zip -d chrs/; done
