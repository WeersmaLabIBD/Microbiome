# eQTL analysis between microbiome and exome data

creator: Shixian

year: 2018

----------------------------------------------------------------------

# 1. Get vcf file 

check the exome data on cluster: /groups/umcg-weersma/tmp04/Shixian/exome/ibd_weersma_wes.vcf.20180117.vcf.gz (with tabix file)

# 2. QC

In the vcf.gz file you will see the header "FILTER", this is about VQSR(variant quality score recalibration) score, which is based the comparsion between your data and know variants 
database on coverage(DP), QualByDepth and etc. It is something like machine learning process to get a balance between sensitivity and specifity using real data and known data. So we 
only keep those variants with "PASS".

sbatch qc_filter.sh

```
#!/bin/bash
#SBATCH --job-name=qc_filter
#SBATCH --error=qc_filter.err
#SBATCH --output=qc_filter.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load VCFtools


vcftools --gzvcf ibd_weersma_wes.vcf.20180117.vcf.gz --recode --remove-filtered-all --out qcpass
```

# 3. split QC pass vcf file into IBD and LLD 

Prepare IBD and LLD individules list: IBD_ID_list.txt and LLD_ID_list.txt
The format of list files only contains one column, listing all samples names. Vcftools, bcftools and plink, all of the software can do this stuff, but only bcftools can recognize
sample names with special symbols like space. That is why we choose bcftools rather than vcftools due to the existing space symbol in IBD sample names. After the splitting, we do
all the folowing analysis in both 2 corhorts separately.

```
208-2491
208-2512
208-2514
208-2515
208-2516
208-4892
208-4947
208-4950
CR 0165,1
CR 0171,1
CR 0172,1
CR 0173,1
CR 0178,1
CR 0220,1
CR 0275.1
CR 0345,1
CR 0362,1
```

mkdir IBD

mkdir LLD

cd IBD

sbatch IBD_split.sh

```
#!/bin/bash
#SBATCH --job-name=IBD_split
#SBATCH --error=IBD_split.err
#SBATCH --output=IBD_split.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load BCFtools

bcftools view -S ./IBD_ID_list.txt -o IBD_split.vcf ../qcpass.recode.vcf
```

cd LLD

sbatch LLD_split.sh

```
#!/bin/bash
#SBATCH --job-name=LLD_split
#SBATCH --error=LLD_split.err
#SBATCH --output=LLD_split.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load BCFtools

bcftools view -S ./LLD_ID_list.txt -o LLD_split.vcf ../qcpass.recode.vcf
```

# 4. Hardy-weinberg disequilibrium test

In natural condition, the allele frequency has reached the Hardy-weinberg banlance. So we do this step to test all the variants to see if they can satisfy this law, to remove 
potential sequencing or calling errors. We only do this in LLD cohort. And remove these variants appear in IBD cohort, too. I tried lots of tools including bedtools, vcftools and 
bcftools but only vcftools works out the stuff perfectly. Because other tools can not remove the exactly same variants between LLD and IBD_HWE_filtered files.

mkdir LLD/split_HWE

cd LLD/split_HWE

sbatch  LLD_HW_00001.sh

```
#!/bin/bash
#SBATCH --job-name=LLD_HW_00001
#SBATCH --error=LLD_HW_00001.err
#SBATCH --output=LLD_HW_00001.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load VCFtools

vcftools --vcf ../LLD_split.vcf --recode --hwe 0.00001 --out LLD_HW_00001

```

mkdir IBD/split_HWE

cd IBD/split_HWE

sbatch IBD_HWE_filter.sh

```
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

vcftools --vcf ../IBD_split.vcf --positions ../../LLD/split_HWE/LLD_HW_00001.recode.vcf --recode --out IBD_HWE_filtered.vcf

```

# 5. Variants calling rate

In this step, we remove those variats can not reach 99% calling rate among the whole corhort.

mkdir LLD/HWE_missing

cd LLD/HWE_missing

sbatch LLD_missing_genotype.sh

```
#!/bin/bash
#SBATCH --job-name=LLD_missing_genotype
#SBATCH --error=LLD_missing_genotype.err
#SBATCH --output=LLD_missing_genotype.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6


module load VCFtools

vcftools --vcf ../split_HWE/LLD_HW_00001.recode.vcf --recode --max-missing 0.99 --out LLD_split_HWE_missing
```

mkdir IBD/HWE_missing

cd IBD/HWE_missing

sbatch IBD_missing_genotype.sh

```
#!/bin/bash
#SBATCH --job-name=IBD_missing_genotype
#SBATCH --error=IBD_missing_genotype.err
#SBATCH --output=IBD_missing_genotype.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6


module load VCFtools

vcftools --vcf ../split_HWE/IBD_HWE_filtered.vcf.recode.vcf --recode --max-missing 0.99 --out IBD_split_HWE_missing
```

# 6. Minor allele frequency

In this step, we filter the variants based on thier MAF. By this way, we remove those low frequency or rare mutation in the whole cohort.

sbatch LLD_MAF_01.sh

```
#!/bin/bash
#SBATCH --job-name=LLD_split_HWE_missing_MAF_01
#SBATCH --error=LLD_split_HWE_missing_MAF_01.err
#SBATCH --output=LLD_split_HWE_missing_MAF_01.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6


module load VCFtools

vcftools --vcf ../HWE_missing/LLD_split_HWE_missing.recode.vcf --recode --maf 0.01 --out LLD_split_HWE_missing_MAF_01
```

the same to the IBD

```
vcftools --vcf IBD_split_HWE_missing.recode.vcf --recode --maf 0.01 --out ../missing_MAF/IBD_split_HWE_missing_MAF_01
```

# 7. LD analysis

From this step, we use PLINK. First of all, we shit the vcf file to PLINK format. Note the parameters used in this step.

mkdir -p IBD/MAF_LD/plink_input

cd IBD/MAF_LD

sbatch IBD_vcf_plink.sh

```
#!/bin/bash
#SBATCH --job-name=IBD_vcf_plink
#SBATCH --error=IBD_vcf_plink.err
#SBATCH --output=IBD_vcf_plink.out
#SBATCH --mem=60gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load plink/1.9-foss-2015b


 plink --vcf /groups/umcg-weersma/tmp04/Shixian/exome/IBD/missing_MAF/IBD_split_HWE_missing_MAF_01.recode.vcf \
  --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out plink_input/IBD_split_HWE_missing_MAF_01
```

module load plink/1.9-foss-2015b

plink --bfile plink_input/IBD_split_HWE_missing_MAF_01 --indep-pairwise 50 5 0.2 --out IBD_split_HWE_missing_MAF_01_indep

plink --bfile plink_input/IBD_split_HWE_missing_MAF_01 --extract IBD_split_HWE_missing_MAF_01_indep.prune.in --genome --make-bed --out IBD_split_HWE_missing_MAF_01_indep

mkdir -p LLD/MAF_LD/plink_input

cd LLD/MAF_LD

sbatch LLD_vcf_plink.sh

```
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
```

module load plink/1.9-foss-2015b

plink --bfile plink_input/LD_split_HWE_missing_MAF_01 --indep-pairwise 50 5 0.2 --out LLD_split_HWE_missing_MAF_01_indep

plink --bfile plink_input/LD_split_HWE_missing_MAF_01 --extract LLD_split_HWE_missing_MAF_01_indep.prune.in --genome --make-bed --out LLD_split_HWE_missing_MAF_01_indep

# 8. pca analysis

In this step, we need to convert the PLINK files (bim, fam and bed) to vcf file again...
We also need to mix our data with 1000 Genome Project.

cd IBD/MAF_LD

plink --bfile IBD_split_HWE_missing_MAF_01_indep --recode vcf --out IBD_split_HWE_missing_MAF_01_indep
 
cd LLD/MAF_LD

plink --bfile LLD_split_HWE_missing_MAF_01_indep --recode vcf --out LLD_split_HWE_missing_MAF_01_indep

mkdir -p exome/pca/1000G

mkdir -p exome/pca/IBD

mkdir -p exome/pca/LLD

sbatch download_1000G.sh

```
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
```

sbatch vcf_plink.sh

```

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

```

mkdir -p pca/IBD/joint_plik

cd pca/IBD/joint_plik

sbatch IBD_1000G_joint.sh

```
#!/bin/bash
#SBATCH --job-name=IBD_1000G_joint
#SBATCH --error=IBD_1000G_joint.err
#SBATCH --output=IBD_1000G_joint.out
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

module load BEDTools/2.25.0-foss-2015b
module load BCFtools
module load VCFtools

bgzip -c /groups/umcg-weersma/tmp04/Shixian/exome/IBD/MAF_LD/IBD_split_HWE_missing_MAF_01_indep.vcf > /groups/umcg-weersma/tmp04/Shixian/exome/IBD/MAF_LD/IBD_split_HWE_missing_MAF_01_indep.vcf.gz

tabix -p vcf /groups/umcg-weersma/tmp04/Shixian/exome/IBD/MAF_LD/IBD_split_HWE_missing_MAF_01_indep.vcf.gz

for chr in {{1..22},X}
do

tabix -h /groups/umcg-weersma/tmp04/Shixian/exome/IBD/MAF_LD/IBD_split_HWE_missing_MAF_01_indep.vcf.gz $chr > IBD_split_HWE_missing_MAF_01_indep.chr$chr.vcf
bgzip -c IBD_split_HWE_missing_MAF_01_indep.chr$chr.vcf > IBD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz
tabix -p vcf IBD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz

done

for chr in {{1..22},X}
do

bedtools intersect -a IBD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz \
                   -b /groups/umcg-weersma/tmp04/Shixian/exome/pca/1000G/chrs/chr$chr.1kg.phase3.v5a.vcf.gz -sorted -header > IBD_intersect_chr$chr.vcf
bedtools intersect -a /groups/umcg-weersma/tmp04/Shixian/exome/pca/1000G/chrs/chr$chr.1kg.phase3.v5a.vcf.gz \
                   -b IBD_split_HWE_missing_MAF_01_indep.chr$chr.vcf.gz -sorted -header > 1000G_intersect_chr$chr.vcf
bgzip -c IBD_intersect_chr$chr.vcf >IBD_intersect_chr$chr.vcf.gz
tabix -p vcf IBD_intersect_chr$chr.vcf.gz

bgzip -c 1000G_intersect_chr$chr.vcf > 1000G_intersect_chr$chr.vcf.gz
tabix -p vcf 1000G_intersect_chr$chr.vcf.gz

vcf-merge 1000G_intersect_chr$chr.vcf.gz IBD_intersect_chr$chr.vcf.gz > intersect_chr$chr.vcf

done

for chr in {1..22} X 
do
plink --vcf  intersect_chr$chr.vcf \
--keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out intersect_chr$chr
done

cat intersect_chrX.fam > IBD_1000G.fam
cat intersect_chr{{1..22},X}.bim > IBD_1000G.bim
(echo -en "\x6C\x1B\x01"; tail -qc +4 intersect_chr{{1..22},X}.bed) > IBD_1000G.bed

/home/umcg-hushixian/bin/gcta64 --grm-bin IBD_1000G --pca 20 --out IBD_1000G --thread-num 10
(echo FID IID PC{1..20}; cat IBD_1000G.eigenvec) > IBD_1000G.pca

```

mkdir -p pca/LLD/joint_plik

cd pca/LLD/joint_plik

sbatch LLD_1000G_joint.sh

```
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

```

put the pca_plot.R, LLD_1000G.pca and IBD_1000G.pca files in one folder, then run the R script

```

pdf('LLD_1000G.pdf')

pca <- read.table('LLD_1000G.pca', header = TRUE)
pop <- read.table('kgp.pop', header = TRUE)
lld=read.table("LLD_ID_list.txt",header = F)
lld$FID=0
lld$pop="LLD"
colnames(lld)=c("IID","FID","POP")
id=rbind(pop,lld)
pca=merge(id,pca)

p <- list()
p[[1]] <- ggplot(pca, aes(x=PC1, y=PC2))
p[[2]] <- ggplot(pca, aes(x=PC1, y=PC3))
p[[3]] <- ggplot(pca, aes(x=PC2, y=PC3))
p[[4]] <- ggplot(pca, aes(x=PC1, y=PC4))
p[[5]] <- ggplot(pca, aes(x=PC2, y=PC4))
p[[6]] <- ggplot(pca, aes(x=PC3, y=PC4))
for (i in 1:6) {
  print(p[[i]] + geom_point( aes(color=POP),size=1) +
          theme_bw())
}

dev.off()

pdf('IBD_1000G.pdf')

pca <- read.table('IBD_1000G.pca', header = TRUE)
pop <- read.table('kgp.pop', header = TRUE)
ibd=read.table("IBD_ID_list.txt",header = F)
ibd$FID=0
ibd$pop="IBD"
colnames(ibd)=c("IID","FID","POP")
id=rbind(pop,ibd)
pca=merge(id,pca)

p <- list()
p[[1]] <- ggplot(pca, aes(x=PC1, y=PC2))
p[[2]] <- ggplot(pca, aes(x=PC1, y=PC3))
p[[3]] <- ggplot(pca, aes(x=PC2, y=PC3))
p[[4]] <- ggplot(pca, aes(x=PC1, y=PC4))
p[[5]] <- ggplot(pca, aes(x=PC2, y=PC4))
p[[6]] <- ggplot(pca, aes(x=PC3, y=PC4))
for (i in 1:6) {
  print(p[[i]] + geom_point( aes(color=POP),size=1) +
          theme_bw())
}

dev.off()
``` 

remove those ancestry outliers

Rcript pca_outliers.Rcript

Then modify the remove_smaple_list mannully because of the special format

```
cd /groups/umcg-weersma/tmp04/Shixian/eqtl/IBD/autosomal_genotypes/

vcftools --remove IBD_pca_remove.txt --vcf IBD_split_HWE_missing_MAF_01_indep.vcf --recode --out IBD_final.vcf

bgzip -c IBD_final.vcf.recode.vcf > IBD_final.vcf.gz
tabix -p vcf IBD_final.vcf.gz

java -Xmx40G -jar /home/umcg-hushixian/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i ./autosomal_genotypes -I VCFFOLDER -o genotypes_trityper -O TRITYPER

cd /groups/umcg-weersma/tmp04/Shixian/eqtl/LLD/autosomal_genotypes/

vcftools --remove LLD_pca_remove.txt --vcf LLD_split_HWE_missing_MAF_01_indep.vcf --recode --out LLD_final.vcf

bgzip -c LLD_final.vcf.recode.vcf > LLD_final.vcf.gz
tabix -p vcf LLD_final.vcf.gz

java -Xmx40G -jar /home/umcg-hushixian/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i ./autosomal_genotypes -I VCFFOLDER -o genotypes_trityper -O TRITYPER
```

prepare eqtl_LLD_coupling.txt eqtl_LLD_taxa eqtl_LLD_taxa.annot

prepare eqtl_IBD_coupling.txt eqtl_IBD_taxa eqtl_IBD_taxa.annot

# 8. eQTL analysis

In both IBD and LLD folder, do the following stuff, NOTE: miQTL_cookbook-master scripts are different between IBD and LLD, check them before running

sh generate_xlm.sh

```
#!/bin/bash


j=0

awk '{print $1}' eqtl_LLD_taxa | while read line

do
  ((j=j+1)) 
  cat miQTL_cookbook-master/software/benchmark_templates/template_benchmark0.xml |perl -pe "s/COHORTNAME/LLD/" > $j\_1_benchmark.xml
  cat $j\_1_benchmark.xml | perl -pe "s/id968/$line/" > $j\_2_benchmark.xml
  cat $j\_2_benchmark.xml | perl -pe "s/genus.Alistipes.id.968.txt/$line.txt/" > $j\_3_benchmark.xml
  cat $j\_3_benchmark.xml | perl -pe "s/genus.Alistipes.id.968.txt.annot/$line.txt.annot/" > $j\_benchmark.xml
  rm $j\_1_benchmark.xml
  rm $j\_2_benchmark.xml
  rm $j\_3_benchmark.xml
done

rm 1_benchmark.xml
```

module load R

Rscript miQTL_cookbook-master/software/step3.1B_prepare_benchmark_data.R eqtl_LLD(IBD)_taxa

```
options = commandArgs(trailingOnly = TRUE)
input_taxonomy = options[1]
input_annotation = paste0(options[1],".annot")
output_folder = "taxa_benchmark_selection"
taxonomy_table = read.table(input_taxonomy,header=T,as.is = T,sep="\t",check.names = F)
colnames(taxonomy_table)[1] = "-"
annot_table = read.table(input_annotation,header=T,as.is = T,check.names = F)

list=read.table(file="eqtl_LLD_taxa",header=T,sep="\t",check.names = F)

taxa=as.vector(list[,1])

dir.create(output_folder)

for (i in taxa){
	tax = taxonomy_table[taxonomy_table[,1]==i,]
	annot = annot_table[annot_table[,2] == i,]
	write.table(tax,file = paste0(output_folder,"/",i,".txt"),sep="\t",quote=F,row.names=F)
	write.table(annot,file = paste0(output_folder,"/",i,".txt.annot"),sep="\t",quote=F,row.names=F)
}

```

check coupling file again to make sure it is correct

sbatch IBD(LLD)_eqtl.sh

```
#!/bin/bash
#SBATCH --job-name=LLD_eqtl
#SBATCH --error=LLD_eqtl.err
#SBATCH --output=LLD_eqtl.out
#SBATCH --mem=5gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

for i in *.xml

do

java -XX:ParallelGCThreads=5 -Xmx30G -jar eqtl-mapping-pipeline-1.4nZ/eqtl-mapping-pipeline.jar --mode metaqtl --settings $i

done
```

plot manhattan

mkdir IBD/manhattan

cd IBD/manhattan

sh LLD_manhattan.sh (also need Manhattan.R)

```
#!/bin/bash
#SBATCH --job-name=LLD_manhattan_input
#SBATCH --error=LLD_manhattan_input.err
#SBATCH --output=LLD_manhattan_input.out
#SBATCH --mem=15gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6

n=1
for i in $(ls /groups/umcg-weersma/tmp04/Shixian/eqtl/re-LLD/LLD_output/)
do
((n=n+1))
zcat /groups/umcg-weersma/tmp04/Shixian/eqtl/re-LLD/LLD_output/$i/eQTLs.txt.gz | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $1,$5}' > $n.LLD_manhattan_input.txt
done

cat *.LLD_manhattan_input.txt > LLD_manhattan_input.txt

module load R
Rscript Manhattan.R
```

```
options(scipen=200)
input=read.table(file="LLD_manhattan_input.txt",header=T, sep="\t", check.names=F)

colnames(input)=c("SNP","CHR","BP","P","taxa")

input$CHR=as.numeric(input$CHR)
input$BP=as.numeric(input$BP)
input$P=as.numeric(as.character(input$P))
#p=input$P
#input$P=p.adjust(p,method="fdr",length(p))

#significant=input[input$P<0.05,]

#significant=significant[order(significant$P,decreasing=F),]

#write.table(significant,file="IBD_significant_result.txt",row.names = F,sep = "\t")

library("qqman")

png("IBD_manhattan.png",width = 1000,height = 800)

manhattan(input,col=c("#A4D3EE","#DDA0DD"))

dev.off()
```

