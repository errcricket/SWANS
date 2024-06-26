import glob
import os
import sys
from sample_list import get_samples

project = sys.argv[1]

flag = ''
def is_empty():
	raw_matrix_path = 'data/endpoints/' + project + '/'
	SAMPLE_LIST = os.listdir(raw_matrix_path)

	if 'analysis' in SAMPLE_LIST:
		loc = SAMPLE_LIST.index('analysis')
		del SAMPLE_LIST[loc]

	for s in SAMPLE_LIST:
		sample_path = 'data'
		path_list = ['endpoints', project, s, '10X']
		
		for p in path_list:
			sample_path = sample_path + '/' + p
		if len(os.listdir(sample_path)) != 0:
			flag = 1
		if len(os.listdir(sample_path)) == 0:
			flag = 0
			print(flag)
	#return(flag)

is_empty()
