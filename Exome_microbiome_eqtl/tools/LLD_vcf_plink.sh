#!/bin/bash
#SBATCH --job-name=LLD_vcf_plink
#SBATCH --error=LLD_vcf_plink.err
#SBATCH --output=LLD_vcf_plink.out
#SBATCH --mem=60gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load plink/1.9-foss-2015b


 plink --vcf  /groups/umcg-weersma/tmp04/Shixian/exome/LLD/missing_MAF/LLD_split_HWE_missing_MAF_01.recode.vcf \
  --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out plink_input/LD_split_HWE_missing_MAF_01
