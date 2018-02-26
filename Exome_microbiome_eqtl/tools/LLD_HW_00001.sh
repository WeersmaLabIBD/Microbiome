#!/bin/bash
#SBATCH --job-name=LLD_HW_00001
#SBATCH --error=LLD_HW_00001.err
#SBATCH --output=LLD_HW_00001.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load VCFtools 

vcftools --vcf ../LLD_split.vcf --recode --hwe 0.00001 --out LLD_HW_00001
