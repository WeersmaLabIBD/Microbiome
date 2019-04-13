# this script is to calculate FDR
# two approaches here: 1) Bonferroni 2) BH

import numpy as np

def p_adjust(pvalue, method):
    p = pvalue
    n = len(p)
    p0 = np.copy(p, order='K')
    nna = np.isnan(p)
    p = p[~nna]
    lp = len(p)
    if method == "bonferroni":
        p0[~nna] = np.fmin(1, lp * p)
    elif method == "fdr":
        i = np.arange(lp, 0, -1)
        o = (np.argsort(p))[::-1]
        ro = np.argsort(o)
        p0[~nna] = np.fmin(1, np.minimum.accumulate((p[o]/i*lp)))[ro]
    else:
        print "Method is not implemented"
        p0 = None
    return p0

import pandas as pd

data=pd.read_table("Merge.assoc.clean.txt", sep="\t")
pvalue=data['p_wald']
pvalue_array = np.array(pvalue)

correction=p_adjust(pvalue_array,"fdr")
data['FDR']=correction

data.to_csv('Merged.FDR.txt', sep='\t', index=False)
