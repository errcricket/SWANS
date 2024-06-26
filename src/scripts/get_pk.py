import sys

file_name = sys.argv[1]
bcms, pks = [], []
difference = [0]
abs_difference = [0]

with open(file_name, 'r') as input_file:
	input_file.readline()

	for line in input_file.readlines():
		sLine = line.replace('\n', '').replace('"', '').split('\t')
		pk, bcm = 0.0, 0.0
		if len(sLine) == 5:
			pk = float(sLine[1])
			bcm = float(sLine[4])
		if len(sLine) == 4:
			pk = float(sLine[0])
			bcm = float(sLine[3])
		pks.append(pk)
		bcms.append(bcm)

max_bcms_location = bcms.index(max(bcms))
pk = pks[max_bcms_location]

print(pk)


'''
for index,b in enumerate(bcms):
	if index+1 > 0 and index+1 < len(bcms):
		diff = bcms[index+1] - b
		a_diff = abs(bcms[index+1] - b)
		difference.append(diff)
		abs_difference.append(a_diff)

max_location = abs_difference.index(max(abs_difference))

pk = 0
if max_location < len(pks)-1:
	if difference[max_location] < 0: 
		max_location = max_location-1

	pk = pks[max_location]
'''

