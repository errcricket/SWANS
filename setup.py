import os
import sys
from sample_list import get_samples

project = sys.argv[1]

SAMPLE_LIST = []
SAMPLE_LIST = get_samples(project)

for s in SAMPLE_LIST:
	sample_path = 'data'
	path_list = ['endpoints', project, s, '10X']
	
	for p in path_list:
		sample_path = sample_path + '/' + p
		if not os.path.exists(sample_path):
			os.mkdir(sample_path)
