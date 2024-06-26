import sys

def get_sig_pc(filename):
	sigPC = 0

	with open(filename, 'r') as input_file:
		input_file.readline()
		for line in input_file.readlines():
			sigPC = int(line.lstrip().rstrip().replace('\n', ''))

	return(sigPC)
