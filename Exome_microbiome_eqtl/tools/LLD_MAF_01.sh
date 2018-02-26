#!/bin/bash
#SBATCH --job-name=LLD_split_HWE_missing_MAF_01
#SBATCH --error=LLD_split_HWE_missing_MAF_01.err
#SBATCH --output=LLD_split_HWE_missing_MAF_01.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6


module load VCFtools

vcftools --vcf ../HWE_missing/LLD_split_HWE_missing.recode.vcf --recode --maf 0.01 --out LLD_split_HWE_missing_MAF_01;