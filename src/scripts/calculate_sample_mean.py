import sys

args = sys.argv
in_files = args[1:]

def calculate_averages(in_file, out_file):
	with open(out_file, 'w') as output_file:
		with open(in_file, 'r') as input_file:
			input_file.readline()
			header = 'gene\tsample_average\n' 
			output_file.write(header)

			for line in input_file.readlines():
				sLine = line.split('\t')
				gene = sLine.pop(0)
				avg = 0

				for s in sLine:
					avg = avg + float(s)

				sample_average = float(avg)/len(sLine)
				out_line = gene + '\t' + str(sample_average) + '\n'
				output_file.write(out_line)

for f in in_files:
	out_file = f.replace('clusters_', '')
	calculate_averages(f, out_file)
