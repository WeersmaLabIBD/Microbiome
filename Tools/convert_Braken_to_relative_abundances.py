#!/bin/python

#Creator: Arnau Vich

#Year: 2016

## USAGE: python step_2_add_missing_levels.py result_corr.txt > result_step_2.txt



import sys
import re
from collections import defaultdict

dom_kraken=0
dom_braken=0
dom_braken_norm=0
king_kraken=0
king_braken=0
king_braken_norm=0
phy_kraken=0
phy_braken=0
phy_braken_norm=0
cla_kraken=0
cla_braken=0
cla_braken_norm=0
ord_kraken=0
ord_braken=0
ord_braken_norm=0
fam_kraken=0
fam_braken=0
fam_braken_norm=0
gen_kraken=0
gen_braken=0
gen_braken_norm=0
sp_kraken=0
sp_braken=0
sp_braken_norm=0

header= "Name" "\t" "taxaid""\t""Taxa_mpa""\t""Kraken""\t""Braken""\t""Braken_norm""\t""Ref_length""\t""Rel_Abund_Kraken""\t""Rel_Abund_Braken""\t""Rel_Abund_Braken_norm"
print header

input_file= open (sys.argv[1],'r')
with open(sys.argv[1]) as input_file:
	next(input_file)
	for i, line in enumerate(input_file):
		each_line=line.split('\t')
		taxa_mpa=each_line[2]
		level=taxa_mpa.split("|")
		if len(level)==1:
			dom_kraken += float(each_line[3])
			dom_braken += float(each_line[4])
			dom_braken_norm += float(each_line[5])
		elif len(level)==2:
			king_kraken += float(each_line[3])
			king_braken += float(each_line[4])
			king_braken_norm += float(each_line[5])
		elif len(level)==3:
			phy_kraken += float(each_line[3])
			phy_braken += float(each_line[4])
			phy_braken_norm += float(each_line[5])
		elif len(level)==4:
			cla_kraken += float(each_line[3])
			cla_braken += float(each_line[4])
			cla_braken_norm += float(each_line[5])
		elif len(level)==5:
			ord_kraken += float(each_line[3])
			ord_braken += float(each_line[4])
			ord_braken_norm += float(each_line[5])
		elif len(level)==6:
			fam_kraken += float(each_line[3])
			fam_braken += float(each_line[4])
			fam_braken_norm += float(each_line[5])
		elif len(level)==7:
			gen_kraken += float(each_line[3])
			gen_braken += float(each_line[4])
			gen_braken_norm += float(each_line[5])
		elif len(level)==8:
			sp_kraken += float(each_line[3])
			sp_braken += float(each_line[4])
			sp_braken_norm += float(each_line[5])
with open(sys.argv[1]) as input_file:
	next(input_file)
	for i, line in enumerate(input_file):
		each_line=line.split('\t')
		taxa_mpa=each_line[2]
		level=taxa_mpa.split("|")
		if len(level)==1:
			rel_kraken= float(float(each_line[3])/dom_kraken) * 100
			rel_braken= float(float(each_line[4])/dom_braken) * 100
			rel_braken_norm= float(float(each_line[5])/dom_braken_norm) * 100
		elif len(level)==2:
			rel_kraken= float(float(each_line[3])/king_kraken) * 100
			rel_braken= float(float(each_line[4])/king_braken) * 100
			rel_braken_norm= float(float(each_line[5])/king_braken_norm) * 100
		elif len(level)==3:
			rel_kraken= float(float(each_line[3])/phy_kraken) * 100
			rel_braken= float(float(each_line[4])/phy_braken) * 100
			rel_braken_norm= float(float(each_line[5])/phy_braken_norm) * 100
		elif len(level)==4:
			rel_kraken= float(float(each_line[3])/cla_kraken) * 100
			rel_braken= float(float(each_line[4])/cla_braken) * 100
			rel_braken_norm= float(float(each_line[5])/cla_braken_norm) * 100
		elif len(level)==5:
			rel_kraken= float(float(each_line[3])/ord_kraken) * 100
			rel_braken= float(float(each_line[4])/ord_braken) * 100
			rel_braken_norm= float(float(each_line[5])/ord_braken_norm) * 100
		elif len(level)==6:
			rel_kraken= float(float(each_line[3])/fam_kraken) * 100
			rel_braken= float(float(each_line[4])/fam_braken) * 100
			rel_braken_norm= float(float(each_line[5])/fam_braken_norm) * 100
		elif len(level)==7:
			rel_kraken= float(float(each_line[3])/gen_kraken) * 100
			rel_braken= float(float(each_line[4])/gen_braken) * 100
			rel_braken_norm= float(float(each_line[5])/gen_braken_norm) * 100
		elif len(level)==8:
			rel_kraken= float(float(each_line[3])/sp_kraken) * 100
			rel_braken= float(float(each_line[4])/sp_braken) * 100
			rel_braken_norm= float(float(each_line[5])/sp_braken_norm) * 100
		name=str(each_line[0])
		taxaid=str(each_line[1])
		name_mpa=str(each_line[2])
		kraken_c=str(each_line [3])
		braken=str(each_line[4])
		norm_braken=str(each_line[5])
		seq_len=re.sub('\n', '', each_line[6])
		print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  % (name, taxaid, name_mpa, kraken_c, braken, norm_braken, seq_len,rel_kraken, rel_braken, rel_braken_norm)
