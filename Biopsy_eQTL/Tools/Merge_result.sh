echo -e "ExpressionGene\tChr\trsID\tPos\tMissingSample\tAllele1\tAllele0\tAllelFre\tBeta\tSE\tlogl_H1\tl_remle\tl_mle\tp_wald\tp_lrt\tp_score" > Merge.assoc.txt

for i in /groups/umcg-gastrocol/tmp04/Inflamed_eQTL/Previous_process/GEMMA_mixed_model/eQTL_CD/output/*.outcome.assoc.txt

do

name=$(basename $i)
expression=${name%.outcome.assoc.txt}
awk -v var="$expression" '{OFS="\t"}{if (NR!=1) print var,$0}' $i >> Merge.assoc.txt

done
