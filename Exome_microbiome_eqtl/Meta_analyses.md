# META ANALYSES

-----------------------------------------------------------------------------------------------------------------------------------------

To perform meta analyses, you need the following files, including:

1) Intersect.txt: overlapped taxonomies and pathways between IBD and LLD cohorts

2) IBD_output and LLD_output file: mbQTL mapping results from each cohort separately

3) meta_analysis.sh: to prepare files for "mental" software

4) change_allele.py: to check the consistency of ref and alt alleles between IBD and LLD (for example: to ensure allele A is ref in both IBD and LLD)

5) metal software installed in your home directory

6) meta_manhattan.sh, Manhattan_1.R and Manhattan_2.R: for plotting manhattan figure, not necessary

```
sbatch meta_analysis.sh

sbatch meta_manhattan.sh

```
