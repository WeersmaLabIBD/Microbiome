###USAGE: python select_tax.py input output
#Creator: Arnau Vich
#Year:2016

#!/usr/bin/python
import sys
import re

input_file= open (sys.argv[1],'r')
output_file= open (sys.argv[2], 'w')
count_line= 1
#pattern="g__"
#pattern="_unclassified"
with open(sys.argv[1]) as input_file:
	for i, line in enumerate(input_file):
		if count_line==1:
			count_line=2
			header = line.split("\t")
			first=header[0]
			more=header[1:]
			#sep='.'
			output_file.write ('%s' % (first))
			for hd in more:
				#clean=hd.split(sep,1)[0]
				output_file.write ('\t' '%s' % (hd))
			output_file.write ('\n')
		else: 
			header = line.split("\t")
			first=header[0]
			#tax_levels= first.split(";")
			tax_levels= first.split("|")
			number_tax_levels= len(tax_levels)
			corrected=number_tax_levels -1 
			last=tax_levels[corrected]  ##ajust
			more=header[1:]
			output_file.write ('%s' % (last))
			for hd in more:
				output_file.write ('\t' '%s' % (hd))
				#output_file.write ('\n')
