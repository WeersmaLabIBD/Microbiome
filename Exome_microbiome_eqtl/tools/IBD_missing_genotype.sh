#!/bin/bash
#SBATCH --job-name=IBD_missing_genotype
#SBATCH --error=IBD_missing_genotype.err
#SBATCH --output=IBD_missing_genotype.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6


module load VCFtools

vcftools --vcf ../split_HWE/IBD_HWE_filtered.vcf.recode.vcf --recode --max-missing 0.99 --out IBD_split_HWE_missing

vcftools --vcf IBD_split_HWE_missing.recode.vcf --recode --maf 0.05 --out ../missing_MAF/IBD_split_HWE_missing_MAF_05
vcftools --vcf IBD_split_HWE_missing.recode.vcf --recode --maf 0.01 --out ../missing_MAF/IBD_split_HWE_missing_MAF_01