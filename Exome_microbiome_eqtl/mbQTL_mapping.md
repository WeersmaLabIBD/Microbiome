#  mbQTL analysis

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
