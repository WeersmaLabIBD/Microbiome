#!/bin/bash
#SBATCH --job-name=LLD_split
#SBATCH --error=LLD_split.err
#SBATCH --output=LLD_split.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load BCFtools 

bcftools view -S ./LLD_ID_list.txt -o LLD_split.vcf ../qcpass.recode.vcf