#!/bin/bash
#SBATCH --job-name=LLD_missing_genotype
#SBATCH --error=LLD_missing_genotype.err
#SBATCH --output=LLD_missing_genotype.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6


module load VCFtools

vcftools --vcf ../split_HWE/LLD_HW_00001.recode.vcf --recode --max-missing 0.99 --out LLD_split_HWE_missing