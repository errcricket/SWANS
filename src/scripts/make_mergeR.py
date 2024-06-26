#Notes: need to add method to read H5 file if a seurat object does not exist.

import sys
import os

from get_condition import get_condition

args = sys.argv
sample_file = args[1]
project = args[2]
project = project.lower()
organism = args[3]
lib_path = args[4] 
sequencing = args[5] 
mito_cutoff = args[6] #mito cut-off for scRNA
storage = args[7] 
#experiment_file = args[8]
min_feature_threshold = args[8]
max_feature_threshold = args[9]

input_data = 'doubletFinder'

#print(sample_file, project, organism, lib_path, project, sequencing, mito_cutoff, storage, experiment_file)

mito = '^MT-'

if organism.lower() == 'mouse':
	 mito = '^mt-'

flag = 0
#if os.path.isfile(experiment_file) == True:
if os.path.isfile(sample_file) == True:
	flag = 1

sample_list = []

with open(sample_file, 'r') as input_file:
	input_file.readline()
	for line in input_file.readlines():
		sLine = line.replace('\n', '').split('\t')
		sample = sLine[0]
		#sample = line.replace('\n', '')
		if sample not in sample_list:
			sample_list.append(sample)

def make_file(sample_file, project, organism):
	string1 = 'library(Seurat, lib.loc=\'' + lib_path + '\')\n'
	string2, string3 = '', ''
	string4 = 'library(patchwork, lib.loc=\'' + lib_path + '\')\n'
	string5 = 'library(dplyr)\n'
	string6 = '# CREATE DIRECTORIES\n#--------------------------------------------------------------------\n# Output that will be unchanged by component #s go here\n'
	string15 = 'base_directory=paste(\'data/endpoints/' + project + '/analysis/\')\n'
	string16 = 'dir.create(base_directory, showWarnings=FALSE)\n\n'
	string17 = 'subdirectory1=paste(base_directory, \'/figures\', sep=\'\')\n'
	string18 = ''
	if storage.lower() == 'h5':
		string2 = 'library(hdf5r, lib.loc=\'' + lib_path + '\')\n'
		string3 = 'library(RIOH5, lib.loc=\'' + lib_path + '\')\n'
		string18 = 'subdirectory2=paste(base_directory, \'/H5\', sep=\'\')\n'
	if storage.lower() == 'rds':
		string18 = 'subdirectory2=paste(base_directory, \'/RDS\', sep=\'\')\n'
	string19 = 'subdirectory3=paste(base_directory, \'/tables\', sep=\'\')\n'
	string20 = 'dir.create(subdirectory1, showWarnings=FALSE)\n'
	string21 = 'dir.create(subdirectory2, showWarnings=FALSE)\n'
	string22 = 'dir.create(subdirectory3, showWarnings=FALSE)\n#--------------------------------------------------------------------\n\n'
	string23 = ' # IMPORT DATA, CREATE SURATE OBJECT, ADD METADATA, CALCULATE MITO%\n\n#--------------------------------------------------------------------\n\n# Generic function to create seurat object\n'
	string24 = 'create_seurat_object <- function()\n{\n'

	combine_list = 'c('
	combine_list2 = 'c('
	for index,s in enumerate(sample_list):
		count = str(index+1)
		if s != sample_list[-1] and s != sample_list[0]:
			combine_list = combine_list + 'S' + count + ','
		if s != sample_list[-1]:
			combine_list2 = combine_list2 + '\'' + 'S' + count + '\','
		else:
			combine_list = combine_list + 'S' + count + ')'
			combine_list2 = combine_list2 + '\'' + 'S' + count + '\')'

	with open('src/scripts/create_merged_dataset.R', 'w') as output_file:
		output_file.write(string1 + string2 + string3 + string4 + string5 + string6)
		output_file.write(string15 + string16 + string17 + string18 + string19 + string20)
		output_file.write(string21 + string22 + string23 + string24)

		count = 0
		for s in sample_list:
			print(s)
			count+=1
			sample = 'S' + str(count)
			object_name = sample + '.data'
			seurat = 'data/endpoints/' + project + '/' + s + '/' + input_data + '/\'\n'
			seurat_string = '\t' + object_name + ' = \'' + seurat
			seurat_string1 = '\t' + sample + ' <- Read10X(data.dir=' + object_name + ')\n'
			seurat_string2 = '\t' + sample + '  <- CreateSeuratObject(counts=' + sample + ', project=\'' + sample + '\', min.cells=3, min.features=' + min_feature_threshold + ')\n'
			seurat_string3 = '\t' + sample + ' <- AddMetaData(' + sample + ', metadata=\'' + s + '\', col.name=\'Experiment\')\n'
			seurat_string3 = seurat_string3 + '\t' + sample + ' <- AddMetaData(' + sample + ', metadata=\'' + s + '\', col.name=\'Samples\')\n\n'
			if flag == 1:
				#seurat_string3 = '\t' + sample + ' <- AddMetaData(' + sample + ', metadata=\'' + get_condition(s, experiment_file) + '\', col.name=\'Experiment\')\n'
				seurat_string3 = '\t' + sample + ' <- AddMetaData(' + sample + ', metadata=\'' + get_condition(s, sample_file) + '\', col.name=\'Experiment\')\n'
				seurat_string3 = seurat_string3 + '\t' + sample + ' <- AddMetaData(' + sample + ', metadata=\'' + s + '\', col.name=\'Samples\')\n\n'
			output_file.write(seurat_string + seurat_string1 + seurat_string2 + seurat_string3)

			seurat_string4 = '\tS.merged <- merge(S1, y = ' + combine_list + ', add.cell.ids = ' + combine_list2 +  ', project = \'' + project + '\')\n'

			seurat_string5, seurat_string6 = '', ''
			s1, s2, s3, s4, s5, s6, s7, s8 = '', '', '', '', '', '', '', ''

			s1 = '\tf1 = paste(subdirectory1, \'/qc_1.pdf\', sep=\'\')\n'
			s2 = '\tf2 = paste(subdirectory1, \'/qc_2.pdf\', sep=\'\')\n'

		seurat_string5 = '\tS.merged[[\'percent.mito\']] <- PercentageFeatureSet(S.merged, pattern = \'' + mito + '\')\n\n'

		s3 = '\tpdf(file=f1, onefile=TRUE, width=11, height=8.5)\n'
		s4 = '\tprint(VlnPlot(S.merged, features = c(\'nFeature_RNA\', \'nCount_RNA\', \'percent.mito\'), ncol = 3))\n'
		s4a = '\tplot1 <- FeatureScatter(S.merged, feature1 = \'nCount_RNA\', feature2 = \'percent.mito\')\n'
		s4b = '\tplot2 <- FeatureScatter(S.merged, feature1 = \'nCount_RNA\', feature2 = \'nFeature_RNA\')\n'
		s4c = '\tprint(plot1 + plot2)\n'
		s5 = '\tdev.off()\n\n'

		seurat_string6 = '\tS <- subset(S.merged, subset = nFeature_RNA > ' + min_feature_threshold + ' & nFeature_RNA < ' + max_feature_threshold + ' & percent.mito < ' + mito_cutoff + ')\n\n'

		s6 = '\tpdf(file=f2, onefile=TRUE, width=11, height=8.5)\n'
		s7 = '\tprint(VlnPlot(S, features = c(\'nFeature_RNA\', \'nCount_RNA\', \'percent.mito\'), ncol = 3))\n'
		s8 = '\tdev.off()\n\n'

		seurat_string7 = ''
		if storage.lower() == 'h5':
			seurat_string7 = '\twrite_h5(data =S, object.type = \'seurat\', assay.name = \'RNA\', file = paste(subdirectory2, ' + '\'/' + project + '_merged_samples.h5\', sep=\'\'))\n'
		if storage.lower() == 'rds':
			seurat_string7 = '\tsaveRDS(S, file = paste(subdirectory2, ' + '\'/' + project + '_merged_samples.RDS\', sep=\'\'))\n'

		seurat_string8 = '}\n\ncreate_seurat_object()'

		output_file.write(seurat_string4 + seurat_string5 + s1 + s3 + s4 + s4a + s4b + s4c + s5 + seurat_string6 + s2 + s6 + s7 + s8 + seurat_string7 + seurat_string8)

make_file(sample_file, project, organism)
