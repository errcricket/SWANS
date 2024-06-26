'''
Author:	E. Reichenberger
Date:		8.18.2021

Purpose:	Score gene sets from CellMaker and PangaloDB to identify clusters.

Steps:

1. Get list of DGE genes from input file   
'''

import glob
import os
import os.path
import sys

args = sys.argv
path = args[1] 
panglaoDB_dir = args[2] #panglaoDB_dir location #panglaoDB_dir = 'reference_data/panglaoDB/'
cellmarker_dir = args[3] #	cellmarker_dir location #cellmarker_dir = 'reference_data/CellMarker/'
threshold = args[4]
str_threshold  = str(threshold)
threshold = float(threshold)
dge_files = list(args[5:])  #path + 'lung_markers_integrated_snn_res.0.8.txt'
#print(dge_files)

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# function to print out list/dictionary items + length
def parse_item(in_item):
	print(in_item)

	if isinstance(in_item, list):
		for i in in_item:
			print(i)

	if isinstance(in_item, dict):
		for i in in_item:
			print(i, in_item[i])

	print('length of item: ', len(in_item))
#---------------------------------------------------------------------------------------
	
#---------------------------------------------------------------------------------------
# function to get unique genes from DGE file into list
def make_DGE_list(in_file, gene_list = []):
	with open(in_file, 'r') as input_file:
		for line in input_file.readlines():
			sLine = line.split('\t')
			gene = sLine[6].replace('\n', '')
		
			if gene not in gene_list:
				gene_list.append(gene)

	return(gene_list)
#	#in the future I may want to consider making this into a dictionary #new_dic[cluster]=[gene_list]
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# function to create dictionary of genes from gene list (dge list) 
# of their mean value across the sample. need gene list and file of means
def make_mean_gene_dic(gene_list, f, dic={}):
	with open(f, 'r') as input_file:
		input_file.readline()
	
		for line in input_file.readlines():
			sLine = line.split('\t')
			gene = sLine[0]
			mean = float(sLine[1].replace('\n', ''))

			if gene in gene_list:
				if gene not in dic:
					dic[gene] = ''
				dic[gene] = mean
	
		return(dic)
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# need to create dictionary of clusters/genes/values
def make_cluster_mean_gene_dic(f, dic={}):
	with open(f, 'r') as input_file:
		header = input_file.readline()
		clusters = header.replace('\n', '').split('\t')

		for c in clusters: #clusters are strings
			if c not in dic:
				dic[c] ={}

		for line in input_file.readlines():
			sLine = line.replace('\n', '').split('\t')
			gene = sLine.pop(0)

			for index,c in enumerate(clusters):
				dic[c][gene] = sLine[index] 

	return{'dictionary':dic, 'clusters':clusters}
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# function to score marker file (similarity == #hits/total#genes)
def score_markers(dge_list, mean_dic, cluster_mean_dic, celltype_dir, out_file, clusters, path=path):
	cell_marker_files = os.listdir(celltype_dir)

	with open(out_file, 'w') as output_file:
		header = 'cell_type\t'
		for c in clusters:
			header = header + c + '\t'
		header = header + '\n'
		output_file.write(header)

		for cell in cell_marker_files:
			celltype = cell.replace('.txt', '')
			cell = celltype_dir + cell

			total_count, match_count = 0, 0

			out_line = celltype + '\t' 
			for c in clusters: #this works but is very time consuming
				with open(cell, 'r') as input_file:
					for line in input_file.readlines():
						total_count+=1
						gene = line.replace('\n', '')

						if ((gene in dge_list) and (gene in mean_dic) and (gene in cluster_mean_dic[c])  and (float(cluster_mean_dic[c][gene]) > float(mean_dic[gene]))):
							match_count+=1

				score = round(100*float(match_count)/total_count, 2)
				out_line = out_line + str(score) + '\t'
			output_file.write(out_line + '\n')

	return(out_file)
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
def print_possible_cell_types(file_list, dge_file, method_dir, threshold=threshold, str_threshold=str_threshold): #, dge_file=dge_file):
	possible_cell_types = dge_file.replace('_markers_', '_possible_cell_types_').replace('dge_plots', 'dataset_characterization').replace(method_dir, '')
	possible_cell_types = possible_cell_types.replace('.txt', '_') + str_threshold + '.txt'

	tool_name = ''
	with open(possible_cell_types, 'w') as output_file:
		pct = []
		marker_dir = ''
		for f in file_list:
			if 'panglaoDB' in f:
				print('pang')
				tool_name =  'panglaoDB'
				marker_dir = panglaoDB_dir
			if 'cellmarker' in f:
				marker_dir = cellmarker_dir 
				print('cm')
				tool_name =  'CellMarker'
				marker_dir = cellmarker_dir
			
			with open(f, 'r') as input_file:
				print(f)
				input_file.readline()
			
				for line in input_file.readlines():
					sLine = line.lstrip().rstrip().replace('\n', '').split('\t')
					celltype = sLine.pop(0)

					for score in sLine:
						try:
							score = float(score)
							if score > threshold:
								out_line = marker_dir + celltype + '.txt' + '\n'
								if out_line not in pct:
									pct.append(out_line)
						except:
							x=9
							#print(score)
		for p in pct:
			output_file.write(p)
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# multi-function calls
for dge_file in dge_files:
	#print(dge_file)
	dge_list = make_DGE_list(dge_file) #list of all differential genes regardless of cluster. 
	#parse_item(dge_list)

	met = ''

	if 'Single' in dge_file:
		met = 'Single'
	if 'Standard' in dge_file:
		met = 'Standard'
	if 'RPCA' in dge_file:
		met = 'RPCA'
	if 'SCT' in dge_file:
		met = 'SCT'

	method_dir = met + '/'
	rna_dataset_mean_file = dge_file.replace('_markers_', '_meanGE_RNA_').replace('dge_plots', 'dataset_characterization').replace(method_dir, '')
	mean_dic = make_mean_gene_dic(dge_list, rna_dataset_mean_file)
	#parse_item(mean_dic)

	rna_cluster_mean_file = dge_file.replace('_markers_', '_meanGE_RNA_clusters_').replace('dge_plots', 'dataset_characterization').replace(method_dir, '')
	make_cluster_mean_gene_dic_values = make_cluster_mean_gene_dic(rna_cluster_mean_file)

	dict_cluster_avg  = make_cluster_mean_gene_dic_values['dictionary']
	clusters  = list(make_cluster_mean_gene_dic_values['clusters'])
	#parse_item(dict_cluster_avg)
					
	score_file_c = dge_file.replace('_markers_', '_celltype_scores_cellmarker_').replace('dge_plots', 'dataset_characterization').replace(method_dir, '')
	score_file_p = dge_file.replace('_markers_', '_celltype_scores_panglaoDB_').replace('dge_plots', 'dataset_characterization').replace(method_dir, '')

	cluster_score_files_p = score_markers(dge_list, mean_dic, dict_cluster_avg, panglaoDB_dir, score_file_p, clusters)
	cluster_score_files_c = score_markers(dge_list, mean_dic, dict_cluster_avg, cellmarker_dir, score_file_c, clusters)
	cluster_list = [cluster_score_files_c, cluster_score_files_p]

	print_possible_cell_types(cluster_list, dge_file, method_dir)
