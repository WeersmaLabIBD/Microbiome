#!/bin/python

#Creator: Arnau Vich

#Year: 2017

## Usage python DUDes_2_mpa.py file.out > file.mpa
## Converts DUDes output to mpa format (2 columns: Bacteria tab rel.abund).

import sys
import re
from collections import defaultdict

## Taxonomical levels in mpa format (used p.e in metaphlan2)
sk="k__"
ph="p__"
cl="c__"
od="o__"
fa="f__"
ge="g__"
sp="s__"
st="t__"
#file_name=sys.argv[1]
## Get file name and remove .out
#file_name = file_name.strip('.out')
file_name = "#SampleID"
input_file= open (sys.argv[1],'r')
count_line= 1
rel_abu="Relative_abundance"
## Print header
print "%s\t%s"  % (file_name, rel_abu)
with open(sys.argv[1]) as input_file:
	## Skip first 6 lines 
	for _ in xrange(6):
		next(input_file)
	for i, line in enumerate(input_file):
		line = line.strip('\n')
		each_line=line.split('\t')
		rank=each_line[1]
		my_tax=each_line[3]
		my_abu=each_line[4]
		my_tax=re.sub(' ', '_',my_tax)
		list_tax=my_tax.split('|')
		if rank == "superkingdom":
			new_my_tax=sk+list_tax[0]
		elif rank == "phylum":
			new_my_tax=sk+list_tax[0] + "|" + ph +	list_tax[1]
		elif rank == "class":		
			new_my_tax=sk+list_tax[0] + "|" + ph +	list_tax[1] + "|" + cl + list_tax[2]
		elif rank == "order":
			new_my_tax=sk+list_tax[0] + "|" + ph +	list_tax[1] + "|" + cl + list_tax[2] + "|" + od + list_tax[3]
		elif rank == "family": 
			new_my_tax=sk+list_tax[0] + "|" + ph +	list_tax[1] + "|" + cl + list_tax[2] + "|" + od + list_tax[3] + "|" + fa + list_tax[4]
		elif rank == "genus":
			new_my_tax=sk+list_tax[0] + "|" + ph +	list_tax[1] + "|" + cl + list_tax[2] + "|" + od + list_tax[3] + "|" + fa + list_tax[4] + "|" + ge + list_tax[5]
		elif rank == "species":
			new_my_tax=sk+list_tax[0] + "|" + ph +	list_tax[1] + "|" + cl + list_tax[2] + "|" + od + list_tax[3] + "|" + fa + list_tax[4] + "|" + ge + list_tax[5] + "|" + sp + list_tax[6]
		elif rank == "strain":
			new_my_tax=sk+list_tax[0] + "|" + ph +	list_tax[1] + "|" + cl + list_tax[2] + "|" + od + list_tax[3] + "|" + fa + list_tax[4] + "|" + ge + list_tax[5] + "|" + sp + list_tax[6] + "|" + st + list_tax[7]
		print "%s\t%s"  % (new_my_tax, my_abu)	
