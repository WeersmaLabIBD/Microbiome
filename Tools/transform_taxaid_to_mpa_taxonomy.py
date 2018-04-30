#!/bin/python

#Creator: Arnau Vich

#Year: 2016


##Create index file. From taxaid to mpa taxonomy like (metaphlan taxonomy)

##### USAGE: python report_statistics.py ../my_index_mpa.txt SizeOfAddedFastas.txt 9002000001352282 > result_corr.txt

import sys
import re
from collections import defaultdict

pattern_domain = re.compile("d__*")
pattern_kingdom = re.compile("k__*")
pattern_phylum = re.compile("p__*")
pattern_class = re.compile("c__*")
pattern_order = re.compile("o__*")
pattern_family= re.compile("f__*")
pattern_genus = re.compile("g__*")
pattern_specie = re.compile("s__*")

log_error = open('missing_taxids.txt', 'w')
mpa_map={}
size_map={}
sufix=sys.argv[3]
suma=0
head=1
levels=["Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Specie"]
header= "Name" "\t" "taxaid""\t""Taxa_mpa""\t""Kraken""\t""Braken""\t""Braken_norm""\t""Ref_length"

def load_mpa_taxonomy():
	input_file= open (sys.argv[1],'r')
	with open(sys.argv[1]) as input_file:
		for i, line in enumerate(input_file):
			each_line=line.split('\t')
			tax_mpa=each_line[0]
			tax_id=str(each_line[1])
			tax_id=re.sub('\n', '', tax_id)
			tax_mpa=chek_and_complete(tax_mpa)
			mpa_map[tax_id]=tax_mpa

def load_size():
	seq_file= open (sys.argv[2],'r')
	with open(sys.argv[2]) as seq_file:
		for i, line in enumerate(seq_file):
			each_line=line.split('\t')
			tax_id=each_line[0]
			genome_size=each_line[1]
			size_map[tax_id]=genome_size
def chek_and_complete(old_mpa):
	mylevels=old_mpa.split('|')
	#print mylevels
	#if len(mylevels)>=2:
#		del mylevels[-1]
	#for i in mylevels:
	#	print i
	# if len(mylevels) ==1:
	last=mylevels[-1]
	#print last
	# else:
	# 	last=mylevels[-1]
	# 	mylevels= mylevels.pop(-1)
	genus="g__unknwon"
	family="f__unknown"
	order="o__unknown"
	clas="c__unknown"
	phylum="p__unknown"
	domain="d__unknown"
	kingdom="k__unknown"
	if pattern_domain.match(last):
		return old_mpa
	elif pattern_kingdom.match(last):
		if len(mylevels) >=2:
			last=mylevels[-1]
			del mylevels [-1]
			for i in mylevels:
				if pattern_domain.match(i):
					domain=i
		else: 
			i=str(mylevels)
			if pattern_domain.match(i):
				domain=i
		old_mpa=domain+"|"+last
		return old_mpa
	elif pattern_phylum.match(last):
		if len(mylevels) >=2:
			last=mylevels[-1]
			del mylevels [-1]
			for i in mylevels:
				if pattern_domain.match(i):
					domain=i
					kingdom=re.sub("d__",'k__',i)
				elif pattern_kingdom.match(i):
					kingdom=i
					name=re.sub("k__",'',kingdom)
					if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
						name="Eukaryota"
					domain="d__"+name
		else:
			i=str(mylevels)
			if pattern_domain.match(i):
				domain=i
				kingdom=re.sub("d__",'k__',i)
			elif pattern_kingdom.match(i):
				kingdom=i
				name=re.sub("k__",'',kingdom)
				if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
					name="Eukaryota"
				domain="d__"+name	
		old_mpa=domain+"|"+kingdom+"|"+last
		return old_mpa
	elif pattern_class.match(last):
		if len(mylevels) >=2:
			last=mylevels[-1]
			del mylevels [-1]
			for i in mylevels:
				if pattern_domain.match(i):
					domain=i
					kingdom=re.sub("d__",'k__',i)
				elif pattern_kingdom.match(i):
					kingdom=i
					name=re.sub("k__",'',kingdom)
					if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
						name="Eukaryota"
					domain="d__"+name
				elif pattern_phylum.match(i):
					phylum=i
		else:
			i=str(mylevels)
			if pattern_domain.match(i):
				domain=i
				kingdom=re.sub("d__",'k__',i)
			elif pattern_kingdom.match(i):
				kingdom=i
				name=re.sub("k__",'',kingdom)
				if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
					name="Eukaryota"
				domain="d__"+name
			elif pattern_phylum.match(i):
				phylum=i
		old_mpa=domain+"|"+kingdom+"|"+phylum+"|"+last
		return old_mpa
	elif pattern_order.match(last):
		if len(mylevels) >=2:
			last=mylevels[-1]
			del mylevels [-1]
			for i in mylevels:
				if pattern_domain.match(i):
					domain=i
					kingdom=re.sub("d__",'k__',i)
				elif pattern_kingdom.match(i):
					kingdom=i
					name=re.sub("k__",'',kingdom)
					if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
						name="Eukaryota"
					domain="d__"+name
				elif pattern_phylum.match(i):
					phylum=i
				elif pattern_class.match(i):
					clas=i
		else:
			i=str(mylevels)
			if pattern_domain.match(i):
				domain=i
				kingdom=re.sub("d__",'k__',i)
			elif pattern_kingdom.match(i):
				kingdom=i
				name=re.sub("k__",'',kingdom)
				if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
					name="Eukaryota"
				domain="d__"+name
			elif pattern_phylum.match(i):
				phylum=i
			elif pattern_class.match(i):
				clas=i
		old_mpa=domain+"|"+kingdom+"|"+phylum+"|"+clas+"|"+last
		return old_mpa
	elif pattern_family.match(last):
		if len(mylevels) >=2:
			last=mylevels[-1]
			del mylevels [-1]
			for i in mylevels:
				if pattern_domain.match(i):
					domain=i
					kingdom=re.sub("d__",'k__',i)
				elif pattern_kingdom.match(i):
					kingdom=i
					name=re.sub("k__",'',kingdom)
					if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
						name="Eukaryota"
					domain="d__"+name
				elif pattern_phylum.match(i):
					phylum=i
				elif pattern_class.match(i):
					clas=i
				elif pattern_order.match(i):
					order=i
		else:
			i=str(mylevels)
			if pattern_domain.match(i):
				domain=i
				kingdom=re.sub("d__",'k__',i)
			elif pattern_kingdom.match(i):
				kingdom=i
				name=re.sub("k__",'',kingdom)
				if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
					name="Eukaryota"
				domain="d__"+name
			elif pattern_phylum.match(i):
				phylum=i
			elif pattern_class.match(i):
				clas=i
			elif pattern_order.match(i):
				order=i
		old_mpa=domain+"|"+kingdom+"|"+phylum+"|"+clas+"|"+order+"|"+last
		return old_mpa
	elif pattern_genus.match(last):
		if len(mylevels) >=2:
			last=mylevels[-1]
			del mylevels [-1]
			for i in mylevels:
				if pattern_domain.match(i):
					domain=i
					kingdom=re.sub("d__",'k__',i)
				elif pattern_kingdom.match(i):
					kingdom=i
					name=re.sub("k__",'',kingdom)
					if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
						name="Eukaryota"
					domain="d__"+name
				elif pattern_phylum.match(i):
					phylum=i
				elif pattern_class.match(i):
					clas=i
				elif pattern_order.match(i):
					order=i
				elif pattern_family.match(i):
					family=i
		else:
			i=str(mylevels)
			if pattern_domain.match(i):
				domain=i
				kingdom=re.sub("d__",'k__',i)
			elif pattern_kingdom.match(i):
				kingdom=i
				name=re.sub("k__",'',kingdom)
				if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
					name="Eukaryota"
				domain="d__"+name
			elif pattern_phylum.match(i):
				phylum=i
			elif pattern_class.match(i):
				clas=i
			elif pattern_order.match(i):
				order=i
			elif pattern_family.match(i):
				family=i
		old_mpa=domain+"|"+kingdom+"|"+phylum+"|"+clas+"|"+order+"|"+family+"|"+last
		return old_mpa
	elif pattern_specie.match(last):
		if len(mylevels) >=2:
			last=mylevels[-1]
			del mylevels [-1]
			for i in mylevels:
				if pattern_domain.match(i):
					domain=i
					kingdom=re.sub("d__",'k__',i)
				elif pattern_kingdom.match(i):
					kingdom=i
					name=re.sub("k__",'',kingdom)
					if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
						name="Eukaryota"
					domain="d__"+name
				elif pattern_phylum.match(i):
					phylum=i
				elif pattern_class.match(i):
					clas=i
				elif pattern_order.match(i):
					order=i
				elif pattern_family.match(i):
					family=i
				elif pattern_genus.match(i):
					genus=i
		else:
			i=str(mylevels)
			if pattern_domain.match(i):
				domain=i
				kingdom=re.sub("d__",'k__',i)
			elif pattern_kingdom.match(i):
				kingdom=i
				name=re.sub("k__",'',kingdom)
				if name=="Fungi" or name=="Metazoa" or name=="Viridiplantae":
					name="Eukaryota"
				domain="d__"+name
			elif pattern_phylum.match(i):
				phylum=i
			elif pattern_class.match(i):
				clas=i
			elif pattern_order.match(i):
				order=i
			elif pattern_family.match(i):
				family=i
			elif pattern_genus.match(i):
				genus=i
		old_mpa=domain+"|"+kingdom+"|"+phylum+"|"+clas+"|"+order+"|"+family+"|"+genus+"|"+last
		return old_mpa

print header
load_mpa_taxonomy();
load_size();

for i in levels:
	my_input= i + "_" + str(sufix) + ".txt"
	arxiu=open (my_input,'r')
	with open (my_input) as arxiu:
		next(arxiu)
		for i, line in enumerate(arxiu):
			my_cols=line.split('\t')
			name=my_cols[0]
			taxaid=my_cols[1]
			name_mpa=mpa_map.get(taxaid)
			kraken=my_cols[3]
			braken=my_cols[5]
			seq_len=size_map.get(taxaid)
			if seq_len is not None:
				norm_braken=int(braken)/float(seq_len)
				suma += float(norm_braken)
			else:
				message="Taxa id" + taxaid + "is missing"
				log_error.write(message)
	with open (my_input) as arxiu:
		next(arxiu)
		for i, line in enumerate(arxiu):
			my_cols=line.split('\t')
			name=my_cols[0]
			taxaid=my_cols[1]
			name_mpa=mpa_map.get(taxaid)
			kraken=my_cols[3]
			braken=my_cols[5]
			seq_len=size_map.get(taxaid)
			if seq_len is not None:
				norm_braken=int(braken)/float(seq_len)
				seq_len=size_map.get(taxaid)
				seq_len=re.sub('\n', '', seq_len)
				#rel_abun=(norm_braken/float(suma))*100
				print "%s\t%s\t%s\t%s\t%s\t%s\t%s"  % (name, taxaid, name_mpa, kraken, braken, norm_braken, seq_len)
			else:
				message="Taxa id" + taxaid + "is missing"
				log_error.write(message)

