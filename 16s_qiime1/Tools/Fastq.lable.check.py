#!/usr/bin/env python

# Modified from Greg Caporaso's code in qiime/split_libraries_fastq.py

# Usage:  python compare_fastq_labels.py fastq1_fp fastq2_fp

from itertools import izip
from sys import argv

from cogent.parse.fastq import MinimalFastqParser

header_index = 0
sequence_index = 1
quality_index = 2

read1 = open(argv[1], "U")
read2 = open(argv[2], "U")

def check_header_match_180_or_later(header1,header2):
    """ Confirm headers are compatible in CASAVA 1.8.0 or later
        
        These contain information on the read number, so can differ
    """
    header1 = header1.split(':')
    header2 = header2.split(':')
    for e1,e2 in zip(header1,header2):
        if e1.split(' ')[0] != e2.split(' ')[0]:
            return False

    return True

print "Printing labels that do not match before space character:"
mismatched_labels_found = False

for read1_data,read2_data in izip(MinimalFastqParser(read1,strict=False),
                                  MinimalFastqParser(read2,strict=False)):

    if not check_header_match_180_or_later(read1_data[header_index], read2_data[header_index]):
        print "Mismatched labels: %s, %s" % (read1_data[header_index], read2_data[header_index])
        mismatched_labels_found = True
    if not (len(read1_data[sequence_index]) == len(read1_data[quality_index])):
        print "Sequence and quality score lengths do not match for read 1 label %s " % read1_data[header_index]
    if not (len(read2_data[sequence_index]) == len(read2_data[quality_index])):
        print "Sequence and quality score lengths do not match for read 2 label %s " % read2_data[header_index]

if not mismatched_labels_found:
    print "No mismatched labels found."
