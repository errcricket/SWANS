import sys
from sample_list import get_samples

project = sys.argv[1]


def test():
	SAMPLE_LIST = str(get_samples(project))
	SAMPLE_LIST = SAMPLE_LIST.replace(',', '').replace('[', '').replace(']', '').replace('\'', '')
	print(SAMPLE_LIST)
	
	return(SAMPLE_LIST)

test()
