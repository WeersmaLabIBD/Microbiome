#!/bin/bash
#SBATCH --job-name=vcf_plink
#SBATCH --error=vcf_plink.err
#SBATCH --output=vcf_plink.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load plink/1.9-foss-2015b

for chr in {1..22} X

 do
 plink --vcf  chrs/chr$chr.1kg.phase3.v5a.vcf.gz \
  --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out chrs/kgp.chr$chr 
 done
