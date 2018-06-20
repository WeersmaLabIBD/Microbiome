#!/usr/bin/python

from __future__ import division
import re
import os
import scipy.stats
import sys

input_name=sys.argv[1]
output_name=sys.argv[1]+".txt"

meta=open(output_name,'a+')
print("SNP", "RefAllele", "NonRefAllele", "P-value", "CorrelationEfficient", "STDERR", "N", "Probe", sep=" ", file=meta)

for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith('SNP'):
        pass
    else:
        ww=re.split(r'[ ]\s*', line)
        snp=ww[0]
        p=ww[4]
        coefficient=ww[5]
        se=ww[6]
        n=ww[7]
        probe=ww[8]
        allele_1=ww[1]
        allele_2=ww[2]
        alt=ww[3]
        if allele_1==alt:
            ref=allele_2
        elif allele_2==alt:
            ref=allele_1
        print(snp, ref, alt, p, coefficient, se, n, probe, sep=" ",file=meta)
meta.close()
