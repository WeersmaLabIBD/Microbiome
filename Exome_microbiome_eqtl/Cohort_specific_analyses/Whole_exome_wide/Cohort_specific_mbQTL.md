# Cohort-specific mbQTLs analyses
--------------------------

This analyses is aiming at identify IBD-influenced mbQTL.

Interaction linear model: bacterial data ~ disease status (0/1) * genotype dosage (0/1/2)

**To perform this analyses, you need the following files, including:**

1) IBD/LLD dosage file: genotype dosage

2) IBD/LLD numeric file: corrected bacterial data

3) IBD/LLD specific file: cohort specific mbQTLs (The criteria follows: discovery cohort p values < 6.83e-7, replication cohort p values > 0.05, significant IBD×genotype p values < 1.06e-3 (Bonferroni method, n=47 tests))

```
Rscript Cohort_specific.R
```
