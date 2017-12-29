#Pipeline for exome data analysis
--------------------------------

Creator: Shixian
Year: 2017

1.Installing the necessary resources and creating workdir
```
sh NGS_DNA-resources.sh
sh makestructure.sh
```
2.creating samplesheet from bam files, the outfile is Exomeproject.csv
```
sh creat_samplesheet.sh -{PATH of samples}
```
3.BamToFastq
```
# using bamTofastq.sh to create script for each bam file

sbatch bamTofastq.sh
for i in $(ll | grep bamTofq | awk '{print $9}')
do
  sbatch $i
done
```
4.Preparing and running NGS_DNA pipeline
```
# copy fastq data to destination
group="umcg-weersma"
tmpName="tmp03"
cp –rf samples.fastq.gz /groups/${group}/${tmpName}/WES-project/rawdata/ngs/

module load NGS_DNA

folder="Exomeproject"  
mkdir /groups/${group}/${tmpName}/WES-project/generatedscripts/$folder/
cd /groups/${group}/${tmpName}/WES-project/generatedscripts/$folder/
cp -rf samplesheet.csv ./
cp -rf ~/github/NGS_DNA/templates/generate_template.sh ./

# NOTE: the latest version of NGS_NDA is 3.4.3, but the Boxy cluster can not load it. We should install it in our {HOME} directory. umcg-weersma.csv must be added

# do some change in generate_template.sh:

1) module load NGS_DNA/3.4.3-beta
2) module load Python/2.7.11-foss-2015b

# errors with module load Python/2.7.11-foss-2015b ：
# Exception in thread "main" java.lang.Exception: Parameter 'batchIDList' used in steps [CopyPrmTmpData, CreateInhouseProjects] does not have a value in any of the parameter files
#        at org.molgenis.compute.model.impl.WorkflowImpl.findMatchingOutput(WorkflowImpl.java:146)
#        at org.molgenis.compute.generators.impl.EnvironmentGenerator.writeEnvironmentToFile(EnvironmentGenerator.java:51)
#        at org.molgenis.compute.generators.impl.EnvironmentGenerator.generate(EnvironmentGenerator.java:76)
#        at org.molgenis.compute.ComputeCommandLine.generate(ComputeCommandLine.java:318)
#        at org.molgenis.compute.ComputeCommandLine.execute(ComputeCommandLine.java:143)
#        at org.molgenis.compute.ComputeCommandLine.main(ComputeCommandLine.java:74)

# errors without module load Python/2.7.11-foss-2015b :
# dos2unix: converting file /groups/umcg-weersma/tmp03/WES-project/generatedscripts/Exomeproject//Exomeproject.csv to UNIX format ...
#  File "/home/umcg-hushixian/github/NGS_DNA//scripts/sampleSheetChecker.py", line 36
#    print "no hpo id's"
#                      ^
#SyntaxError: Missing parentheses in call to 'print'. Did you mean print(nt "no hpo id's")?
#  File "/home/umcg-hushixian/github/NGS_DNA//scripts/gender.py", line 10
#    if genderBool == 'true':
#                           ^
#TabError: inconsistent use of tabs and spaces in indentation
#ls: cannot access *.txt.tmp: No such file or directory
#Samplesize is 2

# 13982 [main] INFO org.molgenis.compute.ComputeCommandLine  - All scripts have been generated
# 13982 [main] INFO org.molgenis.compute.ComputeCommandLine  - You can find them in: /groups/umcg-weersma/tmp03/WES-project/generatedscripts/Exomeproject/scripts

# conclusion: run generate_template.sh with python 2.7 at first and run it again with python 3.5

3) RE-SET {tmpDirectory}, {workDir} and {group}

sh generate_template.sh
cd scripts

# you will find CopyPrmTmpData_0.sh and CreateInhouseProjects_0.sh 2 scripts.
# But CopyPrmTmpData_0.sh can not be run successfully because wrong location of rawdata in prm folder. Therefore, we need to do this job by several commands manually.

sh submit.sh
mv WES-project/rawdata/ngs/* ./projects/Exomeproject/run01/rawdata/ngs/
cd ./projects/Exomeproject/run01/rawdata/ngs/

# delete 0 kb fq files, leave md5 files 

rm -rf *.fq.gz

# remame our rawdata

for sample in *.md5
do
i=$(echo ${sample%.md5})
q=$(echo $i | awk -F "L_" '{print $2}')
  for n in *.fastq.gz
  do 
  m="$(echo ${n%.fastq.gz}).fq.gz"
  if [[  "$q"  ==  "$m"  ]]   
  then
    mv $n ./$i
  fi
  done
done

cd ../../../../projects/Exomeproject/run01/jobs/
```

#debug 

```
# before run all the scripts, we have to debug at first
# s09a_Manta_0.sh: bug 1. unrecognised parameter -a/groups/......  shift to -a /groups/...
#                  bug 2. add command "module load BEDTools/2.25.0-foss-2015b"

for i in s09a_Manta_*.sh
do
   sed -i 's;-a/groups/;-a /groups/;g' $i
   sed -i '11c module load BEDTools/2.25.0-foss-2015b' $i
done

# R package installation: 
 module load R
 install.packages("knitr")
 install.packages("stringr")
 install.packages("stringi")
 install.packages("magrittr")
 install.packages("reticulate")
 install.packages("Rcpp")
 install.packages("jsonlite")
 install.packages("highr")
 
# s13_CoverageCalculations_0.sh: bug . time is too short

for i in s13_CoverageCalculations_*.sh
do
  sed -i 's;--time=05:59:00;--time=96:59:00;g' $i
done
 
# s09e_MantaAnnotation_0.sh: bug . script is killed

for i in s09e_MantaAnnotation_*.sh
do
  sed -i 's;--mem 6gb;--mem 20gb;g' $i
done

```
```
# submit jobs

sh submit.sh

# keep track of the logs below
# molgenis.bookkeeping.log : recording if the scripts run successfully
# molgenis.bookkeeping.walltime : recoding the time cost by each scripts
# molgenis.submitted.log : recording sbatch job numbers
```

#data qulity control
Spike in the PhiX reads and Check the Illumina encoding

```
# The spike-in reads will be inserted in each sample; data will be converted to Illumina 1.9

sh s01_PrepareFastQ_0.sh

# input: 2015-08-13_Illumina_BAM_none_Lnone_214-1206_1.fq.gz 
#output: 2015-08-13_Illumina_BAM_none_Lnone_214-1206_1.phiX.recoded.fq.gz
```
  
Calculate QC metrics on raw data

```
sh s03_FastQC_0.sh

# input: 2015-08-13_Illumina_BAM_none_Lnone_214-1206_1.phiX.recoded.fq.gz
#output: 2015-08-13_Illumina_BAM_none_Lnone_214-1206_1.phiX.recoded_fastqc.html 
#        2015-08-13_Illumina_BAM_none_Lnone_214-1206_1.phiX.recoded_fastqc.zip

```

#data alignment 
Alignment + SortSam
```
# Burrows-Wheeler Aligner (BWA) is used to align the sequencing data to the reference genome. The method used is BWA mem.

sh s04_BwaAlignAndSortSam_0.sh

# input: 2015-08-13_Illumina_BAM_none_Lnone_214-1206_1.phiX.recoded.fq.gz 
#output: 2015-08-13_Illumina_BAM_none_Lnone_214-1206.sorted.bam
```
```
# merge those samples sequenced in multiple lanes

sh s05_MergeBam_0.sh

# input: 2015-08-13_Illumina_BAM_none_Lnone_214-1206.sorted.bam 
#output: 2015-08-13_Illumina_BAM_none_Lnone_214-1206.merged.bam
```
```
# Calculate more accurate base quality scores, the output of this step can be used as an argument in HaplotypeCaller(variant calling)

sh s06_BaseRecalibrator_0.sh

# input: 214-1206.merged.bam 
#output: 214-1206.merged.bam.recalibrated.table

```

#data deduplicates
Deleting duplicates(dedup)
```
# the BAM file is examined to locate duplicate reads

sh s07_MarkDuplicates_0.sh

# input: 214-1206.merged.bam 
#output: 214-1206.merged.dedup.bam
```

# Calculating dedup metrics

sh s08_Flagstat_0.sh

# input: 214-1206.merged.dedup.bam 
#output: 214-1206.merged.dedup.bam.flagstat

```

#SV calling 

```
# In this step, the progam Manta calls all types (DEL,DUP,INV,TRA,INS) from the merged BAM file. Output files are candidateSmallIndels, candidateSV and diploidSV along with information such as difference in length between REF and ALT alleles, type of structural variant end information about allele depth.

sh s09a_Manta_0.sh
```
# CoNVaDING (Copy Number Variation Detection In Next-generation sequencing Gene panels) was designed for small (single-exon) copy number variation (CNV) detection in high coverage next-generation sequencing (NGS) data
# But there is no control group, suitable for high coverage NGS, so it is not necessary
```
```
sh s09b_Convading_0.sh (skipped)
```
```
# XHMM to call copy number variation (CNV) from exome capture, this step is needed, but no control. it seems that it can not find "UMCG/All_Exon_v1" in /apps//data//Controls_Convading_XHMM//Controls_v1.txt

sh s09c_XHMM_0.sh (skipped)
sh s09d_DecisionTree_0.sh (skipped)
```
```
# Annotating manta output
sh s09e_MantaAnnotation_0.sh
# input: 

```

#gender check 
Determine gender
```
# Calculating the coverage on the non pseudo autosomal region and compare this to the average coverage on the complete genome predicts male or female well.

sh s10_GenderCalculate_0.sh

# input: 214-1206.merged.dedup.bam 
#output: 214-1206.merged.dedup.bam.nonAutosomalRegionChrX_hs_metrics
```
```
# gender check

sh s15_GenderCheck_0.sh

# input: 214-1206.merged.dedup.bam.nonAutosomalRegionChrX_hs_metrics
#output: 214-1206.chosenSex.txt

```

#alignment QC 
Calculate coverage per base and per target
```
# Calculates coverage per base and per target, the output will contain chromosomal position, coverage per base and gene annotation

sh s13_CoverageCalculations_0.sh
```
Calculate alignment QC metrics
```
# QC metrics are calculated for the alignment created in the previous steps

sh s14a_CollectMultipleMetrics_0.sh 
sh s14b_CollectHSMetrics_0.sh
sh s14c_CollectBamIndexMetrics_0.sh
sh s14d_CollectGCBiasMetrics_0.sh

```

#variant calling 
Variant discovery

```
# totally 26 batch scripts for each sample: 1-22, Xp, Xnp, Y,and MT. All the exome beds file are from exome_capture_kit: All_Exon_v1. That is why it is important to choose the right kit.
# reference genome and SNP are human_g1k_v37.fa and dbsnp_137.b37.vcf, 
# The b37 conventions were used by the 1KG Project. GATK and IGV use this name while the hg19 conventions were used by the UCSC genome browser.

sh s16a_VariantCalling_*.sh 
java -Xmx12g -jar /path/to/GATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \   ## analysis_type, name of tool to run
    -R human_g1k_v37.fa \   ##  Reference sequence file
    -I merged.dedup.bam \  ## Input file containing sequence data (BAM or CRAM)
    --BQSR merged.bam.calibrated.table \   ## Input covariates table file, for on-the-fly base quality score recalibration
    --dbsnp dbsnp_137.b37.vcf \   ## dbSNP file ##
    --newQuals \
    -o output.g.vcf.gz \
    -L captured.bed \
    --emitRefConfidence GVCF \   ## Mode for emitting reference confidence scores GVCF
                                 ## Reference model emitted with condensed non-variant blocks, i.e. the GVCF format.
    -ploidy 2  ##ploidy 1 in non autosomal chr X region in male
	
# input: 214-1206.merged.dedup.bam
#output: 214-1206.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.g.vcf 
#        It should be 26 vcf files for each sample
```
```
# next there will be a joint analysis over all the samples, combine all samples into 1 according to different batch (1-22, Xp, Xnp, Y, MT)

sh s16d_GenotypeVariants_0.sh

# input: 214-1206.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.g.vcf
#output: Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.genotyped.vcf
#        It should be 26 vcf files for whole project

```

#variant annotation
Annotation

```
# Annotating with SnpEff
# It annotates and predicts the effects of variants on genes (such as amino acid changes)

sh s17_SnpEff_0.sh

# input: Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.genotyped.vcf
#output: Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.snpeff.vcf
#        It should be 26 vcf files for whole project
```
```
# Annotating with VEP
# VEP determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions
sh ????????
```
```
# Annotating with CADD, GoNL, ExAC (CmdLineAnnotator)
# Data will be annotated with CADD, GoNL and ExAC

sh s18_CmdLineAnnotator_0.sh

# input: Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.snpeff.vcf
#output: Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.snpeff.exac.vcf
#        Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.snpeff.exac.gonl.vcf
#        Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.snpeff.exac.gonl.cadd.vcf
```
Merge batches
```
# Running GATK CatVariants to merge all the batch files created in the previous into 1

sh s19_MergeBatches_0.sh

# input: Exomeproject.${batchBed 1-22, Xp, Xnp, Y, MT}.variant.calls.snpeff.exac.gonl.cadd.vcf
#output: Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf
#        all batched and samples are in 1 file

```

#variant filteration 
Split indels and SNPs

```
# This step is necessary because the filtering of the vcf needs to be done seperately.

sh s20_SplitIndelsAndSNPs_0.sh

# input: Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf
#output: Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.indels.vcf
#        Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.snps.vcf
#        214-1206.annotated.indels.vcf
#        214-1206.annotated.snps.vcf
```
SNP and (b) Indel filtration
```
# Based on certain quality thresholds (based on GATK best practices) the SNPs and indels are filtered or marked as Pass.

sh s21a_SnpFiltration_0.sh
sh s21b_IndelFiltration_0.sh

# input: Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.indels.vcf
#        Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.snps.vcf
#        214-1206.annotated.indels.vcf
#        214-1206.annotated.snps.vcf
#output: Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.filtered.snps.vcf
#        Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.filtered.indels.vcf
#        214-1206.annotated.filtered.snps.vcf
#        214-1206.annotated.filtered.indels.vcf
```
Merge indels and SNPs
```
```
# Merge all the SNPs and indels into one file (per project) and merge SNPs and indels per sample.
sh s22_MergeIndelsAndSnps_0.sh
# input: Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.filtered.snps.vcf
#        Exomeproject.variant.calls.GATK.sorted.vcf.gz.vcf.annotated.filtered.indels.vcf
#        214-1206.annotated.filtered.snps.vcf
#        214-1206.annotated.filtered.indels.vcf
#output: 214-1206.final.vcf
#        Exomeproject.final.vcf

```

#Gavin analysis
Gavin split samples
```
# Tool that predict the impact of the SNP with the help of different databases (CADD etc)

sh s20_SplitIndelsAndSNPs_0.sh

# input: 214-1206.final.vcf
# processing files: 214-1206.GAVIN.RVCF.firstpass.vcf
#                   214-1206.GAVIN.toCadd.tsv
#                   214-1206.GAVIN.fromCadd.tsv
#                   214-1206.GAVIN.RVCF.final.vcf
#                   214-1206.GAVIN.RVCF.final.mergedWithOriginal.vcf
# final output: 214-1206.GAVIN.rlv.vcf
```


#gene network analysis
GeneNetwork
···
# Tool that ranks genes based on HPO ID's (focus on disease phenotype), there will be 2 extra INFO fields in the output vcf with a score and position.
# The Human Phenotype Ontology (HPO) aims to provide a standardized vocabulary of phenotypic abnormalities encountered in human disease. 

sh s23b_GeneNetwork_0.sh (skipped)

```
#pipeline QC 
QC for data analysis
```
# In silico concordance check

sh s26_InSilicoConcordance_0.sh

# tail -4 "//apps//data//inSilico/humanPhiX/InSilicoData.chrNC_001422.1.variant.calls.vcf" > "/groups/umcg-weersma//tmp03//tmp//Exomeproject/run01//InSilico.txt"
# output: inSilicoConcordance.txt
```
```
# collecting metrics

sh s27a_QCStats_0.sh

#output: 214-1206.total.qc.metrics.table
#        Exomeproject_qc.metricsList
```
```
# Generate quality control report

sh s27c_MultiQC_0.sh

# input: Exomeproject_qc.metricsList
#output: Exomeproject_multiqc_report.html
```