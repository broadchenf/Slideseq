#!/usr/bin/python
import sys,os
#convert txt file to csv file
#usage 'python txt_to_csv.py directorypath'
def process(filename):
	with open(sys.argv[1]+'/'+filename,'r') as f:
		lines_after_3 = f.readlines()[3:]#get rid of the header
		f = []
		lines_after_3 = (line.replace("\t", ",")for line in lines_after_3 if line)
		with open(sys.argv[1]+'/'+filename + '.csv', 'w') as out_file:
			for line in lines_after_3:
				out_file.write(line)

for f in os.listdir(sys.argv[1]):
	process(f)
										     
