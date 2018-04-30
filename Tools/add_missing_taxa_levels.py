#!/bin/python

# Creator: Arnau Vich
# Year: 2016

##Create index file. From taxaid to mpa taxonomy like (metaphlan taxonomy)

## USAGE: python step_2_add_missing_levels.py result_corr.txt > result_step_2.txt


import sys
import re
from collections import defaultdict

header= "Name" "\t" "taxaid""\t""Taxa_mpa""\t""Kraken""\t""Braken""\t""Braken_norm""\t""Ref_length"
#Start dictionaries 
kraken={}
braken={}
braken_norm={}
new_kraken={}
new_braken={}
new_braken_norm={}
pattern_unknown = re.compile(".*__unknown")
#load input file in a reference dictionary
def tax_dict():
	input_file= open (sys.argv[1],'r')
	with open(sys.argv[1]) as input_file:
		#Skip header
		next(input_file)
		for i, line in enumerate(input_file):
			#Split columns tab separated
			each_line=line.split('\t')
			#assign taxonomy to 3 column (remember: python numeration is 0,1,2...)
			tax_mpa=each_line[2]
			#assign counts to variables
			kraken_counts=int(each_line[3])
			braken_counts=int(each_line[4])
			braken_norm_counts=float(each_line[5])
			##Create dictionaries: taxonomies as indexes and counts as values. 
			kraken[tax_mpa]=kraken_counts
			braken[tax_mpa]=braken_counts
			braken_norm[tax_mpa]=braken_norm_counts
#def taxa_iterations():

print header
#Run the previous function
tax_dict();
#for key , value in kraken.iteritems():
# 	print key, value
input_file= open (sys.argv[1],'r')
with open(sys.argv[1]) as input_file:
	next(input_file)
	for i, line in enumerate(input_file):
		my_line=line.split('\t')
		my_kraken_counts=int(my_line[3])
		my_braken_counts=int(my_line[4])
		my_braken_norm_counts=float(my_line[5])
		testing_tax=my_line[2]
		testing_tax_len=testing_tax.split('|')
		while len(testing_tax_len)>1:
			#del testing_tax_len [-1]
			#testing_tax='|'.join(testing_tax_len)
			tmp_len=testing_tax.split('|')
			last=tmp_len[-1]
			testing_tax="|".join(testing_tax.split('|')[:-1])
			#testing_tax=str(testing_tax)
			#print testing_tax
			if testing_tax in kraken:
				break
				#print "TRUE",testing_tax
			else:
				#print "FALSE", testing_tax
				if testing_tax in new_kraken:
					if pattern_unknown.match(last):
						#print "true"
					 	new_kraken[testing_tax]+=my_kraken_counts
					 	new_braken[testing_tax]+=my_braken_counts
					 	new_braken_norm[testing_tax]+=my_braken_norm_counts
					else:
						new_kraken[testing_tax] += my_kraken_counts
						new_braken[testing_tax] +=my_braken_counts
						new_braken_norm[testing_tax] +=my_braken_norm_counts
						#break
				else:
					new_kraken[testing_tax]=my_kraken_counts
					new_braken[testing_tax]=my_braken_counts
					new_braken_norm[testing_tax]=my_braken_norm_counts
			testing_tax_len=testing_tax.split('|')
		name=str(my_line[0])
		taxaid=str(my_line[1])
		name_mpa=str(my_line[2])
		kraken_c=str(my_line[3])
		braken=str(my_line[4])
		norm_braken=str(my_line[5])
		seq_len=re.sub('\n', '', my_line[6])
		print "%s\t%s\t%s\t%s\t%s\t%s\t%s"  % (name, taxaid, name_mpa, kraken_c, braken, norm_braken, seq_len)
for key, value in new_kraken.iteritems():
	name="added_level"
	taxaid="-"
	name_mpa=key
	kraken_c=value
	braken=new_braken[key]
	norm_braken=new_braken_norm[key]
	seq_len="-"
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s"  % (name, taxaid, name_mpa, kraken_c, braken, norm_braken, seq_len)
