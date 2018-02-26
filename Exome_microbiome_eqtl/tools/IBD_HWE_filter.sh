#!/bin/bash
#SBATCH --job-name=HWE_filter
#SBATCH --error=HWE_filter.err
#SBATCH --output=HWE_filter.out
#SBATCH --mem=80gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load BEDTools/2.25.0-foss-2015b
module load BCFtools
module load VCFtools

# bedtools intersect -a ../IBD_split.vcf -b ../../LLD/split_HWE/LLD_HW_00001.recode.vcf -wa  -sorted> IBD_HWE_filtered.vcf

# bcftools isec -p isec_output ../IBD_split.vcf ../../LLD/split_HWE/LLD_HW_00001.recode.vcf

vcftools --vcf ../IBD_split.vcf --positions ../../LLD/split_HWE/LLD_HW_00001.recode.vcf --recode --out IBD_HWE_filtered.vcf