#! /usr/bin/env python

# Script to change the read labels from a fastq file

Usage = """
fastq_header_changer.py - version 3.0 created by Brett Younginger
Convert the headers of a fastq file to match other fastq headers
for usearch commands

Usage:
  fastq_header_changer.py infile.fastq > converted.fq
"""
import sys
import re

SearchStr = r'^(@M02149.+AAMFM:1:)(\d+):(\d+):(\d+);(sample=)(BY[\d-]+);' # Modify this search string if needed
ReplaceStr = r'@\6.\2\3\4'

if len(sys.argv)<2:
	print(Usage)
else:
	InFileName = sys.argv[1]
	
	InFile = open(InFileName, 'r') # Not sure what the 'U' does

	RegSub = re.compile(SearchStr)
	
	for Line in InFile:
		Line = Line.strip('\n')
		if Line.startswith('@M02149'):
			NewLine = RegSub.sub(ReplaceStr,Line)
			print(NewLine)
		else:
			print(Line)
	InFile.close()		
