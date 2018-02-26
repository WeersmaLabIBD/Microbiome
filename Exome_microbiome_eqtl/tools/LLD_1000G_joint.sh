#!/bin/bash
#SBATCH --job-name=LLD_1000G_joint
#SBATCH --error=LLD_1000G_joint.err
#SBATCH --output=LLD_1000G_joint.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load BEDTools/2.25.0-foss-2015b
module load BCFtools
module load VCFtools
module load plink/1.9-foss-2015b

bgzip -c /groups/umcg-weersma/tmp04/Shixian/exome/LLD/MAF_LD/LLD_split_HWE_missing_MAF_01_indep.vcf > /groups/umcg-weersma/tmp04/Shixian/exome/LLD/MAF_LD/LLD_split_HWE_missing_MAF_01_indep.vcf.gz

tabix -p vcf /groups/umcg-weersma/tmp04/Shixian/exome/LLD/MAF_LD/LLD_split_HWE_missing_MAF_01_indep.vcf.gz

for chr in {{1..22},X}
do

tabix -h /groups/umcg-weersma/tmp04/Shixian/exome/LLD/MAF_LD/LLD_split_HWE_missing_MAF_01_indep.vcf.gz $chr > LLD_split_HWE_missing_MAF_01_indep.chr$chr.vcf
bgzip -c LLD_split_HWE_missing_MAF_01_indep.chr$chr.vcf > LLD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz
tabix -p vcf LLD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz

done

for chr in {{1..22},X}
do

bedtools intersect -a LLD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz \
                   -b /groups/umcg-weersma/tmp04/Shixian/exome/pca/1000G/chrs/chr$chr.1kg.phase3.v5a.vcf.gz -sorted -header > LLD_intersect_chr$chr.vcf
bedtools intersect -a /groups/umcg-weersma/tmp04/Shixian/exome/pca/1000G/chrs/chr$chr.1kg.phase3.v5a.vcf.gz \
                   -b LLD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz -sorted -header > 1000G_intersect_chr$chr.vcf
bgzip -c LLD_intersect_chr$chr.vcf >LLD_intersect_chr$chr.vcf.gz
tabix -p vcf LLD_intersect_chr$chr.vcf.gz

bgzip -c 1000G_intersect_chr$chr.vcf > 1000G_intersect_chr$chr.vcf.gz
tabix -p vcf 1000G_intersect_chr$chr.vcf.gz

vcf-merge 1000G_intersect_chr$chr.vcf.gz LLD_intersect_chr$chr.vcf.gz > intersect_chr$chr.vcf

done

for chr in {1..22} X 
do
plink --vcf  intersect_chr$chr.vcf \
--keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out intersect_chr$chr
done

cat intersect_chrX.fam > LLD_1000G.fam
cat intersect_chr{{1..22},X}.bim > LLD_1000G.bim
(echo -en "\x6C\x1B\x01"; tail -qc +4 intersect_chr{{1..22},X}.bed) > LLD_1000G.bed

/home/umcg-hushixian/bin/gcta64 --grm-bin LLD_1000G --pca 20 --out LLD_1000G --thread-num 10
(echo FID IID PC{1..20}; cat LLD_1000G.eigenvec) > LLD_1000G.pca
