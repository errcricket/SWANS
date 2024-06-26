import sys
from sample_list import get_samples

args = sys.argv
project = args[1]

samples, estimates = [], []

samples = get_samples(project)

'''
with open('samples.sample_list', 'r') as input_file:
	input_file.readline()

	for line in input_file.readlines():
		sample = line.split('\t')[0]
		
		if sample not in samples:
			samples.append(sample)
'''


def get_estimate(samples, estimates):
	estimate_file = 'html_report_' + project + '/estimates.txt'
	with open(estimate_file, 'w') as output_file:
		path = 'data/endpoints/' + project + '/'
		for s in samples:
			html_path = path + s + '/10X/outs/web_summary.html'

			with open(html_path, 'r') as input_file:
				#print(html_path)
				for line in input_file.readlines():
					line = line.lstrip()
					if line.startswith('const data'):
						sLine = line.split('annotated_cells')[1].split('Estimated number of cells')[0]
						#": {"threshold": "pass", "metric": "11,222", "name": "

						estimate_split = sLine.replace('"', '').replace(':', '').replace('{', '').replace(',', '').split(' ')
						estimate = ''
				
						for e in estimate_split:
							if e.isdigit():
								estimate = e
								if estimate not in estimates:
									estimates.append(estimate)
									output_file.write(s + '\t' + e + '\n')
					
								
	return(estimates)

estimate_list = get_estimate(samples, estimates)
