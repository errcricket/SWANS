'''
Author:	E. Reichenberger
Date:		10.17.2019		

Purpose: Re-usable script to determine which samples will be analyzed. If samples.sample_list is present, get the files from there, if not, do a greedy glob search. 

Required: 
'''

import glob
import os.path
import os
import sys
from pathlib import Path   

#files
samplefile_list = 'samples.sample_list'

project = sys.argv[1]
#print(project)

def get_samples(project):
	raw_matrix_path = ''

	#paths
	raw_matrix_path = 'data/endpoints/' + project + '/'

	SAMPLE_LIST = [] #list variables
	sample_path_dic = {}

	if os.path.isfile(samplefile_list) == True:
		with open(samplefile_list, 'r') as input_file:
			input_file.readline()
			for line in input_file.readlines():
					line = line.lstrip().rstrip()
					if len(line) > 1: #grrr....check for empty lines
						sLine = line.split('\t')
						sample = sLine[0]
						condition = sLine[1]
						path = sLine[2] + '/outs/'
						print(path)

						SAMPLE_LIST.append(sample)
						if sample not in sample_path_dic:
							sample_path_dic[sample] = path

		print(SAMPLE_LIST)
		with open('make_symbolic_links.sh', 'w') as output_file:
			for s in SAMPLE_LIST:
				path_list = 'data/endpoints/' + project + '/' +  s  + '/10X/'
				out_line = 'cd ' + path_list + '\n'
				output_file.write(out_line)
				out_line = 'ln -s ' + sample_path_dic[s] + ' . \n'
				output_file.write(out_line)
				out_line = 'cd - \n'
				output_file.write(out_line)
		
get_samples(project)
