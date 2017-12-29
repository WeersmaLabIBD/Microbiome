#!/bin/bash

# $1=samples path

if [[ $# != 1 ]]
then
  break
  echo "Please input the path for samples"
fi

module load SAMtools

sam_path=$1

echo -e "externalSampleID,project,sequencer,sequencingStartDate,run,flowcell,lane,seqType,capturingKit,arrayFile,arrayID,barcode,externalFastQ_1,externalFastQ_2,Gender" > Exomeproject.csv

dir=$(ls -l $sam_path | awk '{print $9}' | grep bam)
cd $sam_path
for i in $dir
  do
  externalSampleID=$(echo $i | cut -d "." -f1)
  project="Exomeproject"
  sequencer="Illumina"
  m=$(samtools view -H $i | grep "DT:" |sed -n '1,1p'| awk -F '[\t:]' '{print $15}')
  sequencingStartDate=$(echo $m | cut -d "T" -f1)
  run="BAM"
  flowcell=""
  lane=""
  seqType="PE"
  capturingKit="UMCG\/All_Exon_v1"
  arrayFile=""
  arrayID=""
  barcode=$externalSampleID
  externalFastQ_1=$externalSampleID\_1.fastq.gz
  externalFastQ_2=$externalSampleID\_2.fastq.gz
  Gender="Unknown"
  
  echo -e "$externalSampleID,$project,$sequencer,$sequencingStartDate,$run,$flowcell,$lane,$seqType,$capturingKit,$arrayFile,$arrayID,$barcode,$externalFastQ_1,$externalFastQ_2,$Gender" >> Exomeproject.csv
done

  

