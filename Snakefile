'''
Author:	E. Reichenberger
Date:		4.25.2022

Notes: Reference script to call all rules
'''

import sys
import glob
import os
import os.path
import subprocess
from sample_list import get_samples
from get_sigPC import get_sig_pc

include: 'src/scripts/Snakefile_config.py'
configfile: 'configs/local_configs.yaml'
configfile: 'configs/final_analysis_configs.yaml'

#local configs-------------------------
RPATH = config['RPATH']
PROJECT = config['PROJECT']
PROJECT = PROJECT.lower() #guarantees lower case
MATRIX = config['MATRIX']
COMPONENTS = config['ICOMPONENTS']
ORGANISM = config['ORGANISM']
SEQUENCING = config['SEQUENCING']
RESOLUTION = config['RESOLUTION']
RESOLUTIONS = []
if RESOLUTION != 'I':
	p_counts = str(RESOLUTION).count('.')
	if p_counts > 1:
		RESOLUTIONS = RESOLUTION.split(' ')

CONSERVED_GENES = config['CONSERVED_GENES']
MITO = config['MITO']
MITO_REGRESSION = config['MITO_REGRESSION']
STORAGE = config['STORAGE']
STORAGE = 'rds' #H5 does not work
METHOD = config['METHOD']
if METHOD.startswith('A') or METHOD.startswith('a'):
	METHOD = METHOD.upper() #guarantees 'ALL'
MIN_FEATURE_THRESHOLD = config['MIN_FEATURE_THRESHOLD']
MAX_FEATURE_THRESHOLD = config['MAX_FEATURE_THRESHOLD']
CCREGRESSION = config['CCREGRESSION']
CCREGRESSION_METHOD = config['CCREGRESSION_METHOD']
MARKER_THRESHOLD = config['MARKER_THRESHOLD']
USER_GENE_FILE = config['USER_GENE_FILE']
#--------------------------------------

#final configs-------------------------
FINAL_METHOD = config['FINAL_METHOD']
FINAL_RESOLUTION = config['FINAL_RESOLUTION']
CLUSTER_ANNOTATION_FILE = config['CLUSTER_ANNOTATION_FILE']
USER_SUPPLIED_SEURAT_OBJECT = config['USER_SUPPLIED_SEURAT_OBJECT']
PATH_USER_SUPPLIED_SEURAT_OBJECT = config['PATH_USER_SUPPLIED_SEURAT_OBJECT']
#--------------------------------------

s_file = 'samples.sample_list'
if os.path.exists(s_file) == False:
	print('You must supply a list of files')
	print('Create samples.sample_list first, then run this script again')
	sys.exit()
if os.path.exists(s_file) == True:
	SAMPLE_LIST = []
	SAMPLE_LIST = get_samples(PROJECT) #get sample list from samples.sample_list or greedy glob.
	#print(SAMPLE_LIST)
	#-------------------------------------------------------------------------------------

	#Storage function
	def storage(in_file, STORAGE):
		if STORAGE.lower() == 'h5':
			in_file = in_file + '.h5Seurat'
		if STORAGE.lower() == 'rds':
			in_file = in_file + '.RDS'

		return(in_file)
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	# SoupX (src/rules/soupX.rules) -----------------------------------------------
	soupX_path = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX', s=SAMPLE_LIST)

	# soupX output
	barcodes_path_sx = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX/barcodes.tsv.gz', s=SAMPLE_LIST)
	features_path_sx = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX/features.tsv.gz', s=SAMPLE_LIST)
	matrix_path_sx = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX/matrix.mtx.gz', s=SAMPLE_LIST)
	soupX_list = [barcodes_path_sx, features_path_sx, matrix_path_sx]
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	# DoubletFinder (src/rules/doubletFinder.rules) -----------------------------------------------
	doubletFinder_path = expand('data/endpoints/' + PROJECT + '/' + '{s}/doubletFinder', s=SAMPLE_LIST)

	# doubletFinder output
	barcodes_path_df = expand('data/endpoints/' + PROJECT + '/' + '{s}/doubletFinder/barcodes.tsv.gz', s=SAMPLE_LIST)
	features_path_df = expand('data/endpoints/' + PROJECT + '/' + '{s}/doubletFinder/features.tsv.gz', s=SAMPLE_LIST)
	matrix_path_df = expand('data/endpoints/' + PROJECT + '/' + '{s}/doubletFinder/matrix.mtx.gz', s=SAMPLE_LIST)
	doubletFinder_list = [barcodes_path_df, features_path_df, matrix_path_df]

	# MAY WANT TO ADD pk files
	cds_input = []
	for s in SAMPLE_LIST:
		bars = 'data/endpoints/' + PROJECT + '/' + s + '/doubletFinder/barcodes.tsv.gz'
		features = 'data/endpoints/' + PROJECT + '/' + s + '/doubletFinder/features.tsv.gz'
		matrix = 'data/endpoints/' + PROJECT + '/' + s + '/doubletFinder/matrix.mtx.gz'
		cds_input.append(bars)
		cds_input.append(features)
		cds_input.append(matrix)
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	# create_dynamic_script.rules output
	# src/scripts/create_merged_dataset.R
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	# Merge (src/rules/merge_samples.rules) -----------------------------------------------
	little_project = PROJECT.lower()
	out_parts = []
	if STORAGE.lower() == 'h5':
		out_parts = ['endpoints', PROJECT, 'analysis', 'H5']
	if STORAGE.lower() == 'rds':
		out_parts = ['endpoints', PROJECT, 'analysis', 'RDS']

	string_path_m = 'data'
	for o in out_parts:
		string_path_m = string_path_m + '/' + o
		if not os.path.exists(string_path_m):
			os.mkdir(string_path_m)
			
	storage_file = ''
	if STORAGE.lower() == 'h5':
		storage_file = string_path_m + '/' + little_project + '_merged_samples.h5Seurat'
	if STORAGE.lower() == 'rds':
		storage_file = string_path_m + '/' + little_project + '_merged_samples.RDS'
		
	figure_outs = 'data/endpoints/' + PROJECT + '/analysis/figures/'

	# merg_samples output
	f1 = figure_outs + 'qc_1.pdf'
	f2 = figure_outs + 'qc_2.pdf'
	merge_list = [f1, f2, storage_file]
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	# Get number sig. components ---------------------------------------------------------
	length = len(SAMPLE_LIST)
	r_file = ''
	if length == 1:
		r_file = 'src/scripts/calculate_sigPCs_1.R'
		METHOD = 'Single'
	if length > 1:
		r_file = 'src/scripts/calculate_sigPCs.R'
		
	out_parts = []
	if STORAGE.lower() == 'h5':
		out_parts = ['endpoints', PROJECT, 'analysis', 'H5']
	if STORAGE.lower() == 'rds': 
		out_parts = ['endpoints', PROJECT, 'analysis', 'RDS']
		
	string_path_s = 'data'
	for o in out_parts:
		string_path_s = string_path_s + '/' + o
		
	storage_file = ''
	if STORAGE.lower() == 'h5':
		storage_file = string_path_s + '/' + PROJECT.lower() + '_merged_samples.h5Seurat'
	if STORAGE.lower() == 'rds':
		storage_file = string_path_s + '/' + PROJECT.lower() + '_merged_samples.RDS'
		
	project_tables_path = 'data/endpoints/' + PROJECT + '/analysis/tables/'

	# get_sig_PC.rules output
	signicant_pc = project_tables_path + 'sigPC.txt'
	signicant_pcs = project_tables_path + 'sigPCs.txt'
	sigs = [signicant_pc, signicant_pcs]
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	# Define variables needed once signicant_pc exists.
	# GLOBAL
	PCA, COMPONENTS_S = '', ''
	#-------------------------------------------

	# analyze_sc_object.rules -------
	sc_r_file = ''
	if len(SAMPLE_LIST) > 1:
		sc_r_file = 'src/scripts/analyze_sc_object.R'
	if len(SAMPLE_LIST) == 1:
		sc_r_file = 'src/scripts/analyze_sc_object_1.R'
		
	sc_objects, integrated_seurat_objects = [], []
	integrated_seurat_object, meta_file = '', ''
	#-------------------------------------------

	# dge_plots.rules ------------
	resolution_file = ''
	dge_files, images_files, marker_files, marker_100files, proportion_files, proportion_images = [], [], [], [], [], []
	migsig_path = ''

	resolution_list = []
	resolution_file = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/' + PROJECT.lower() + '_resolution_list.txt' #PCA == '' this is ok
	#-------------------------------------------

	# calculate_means_plot_markers.rules ---------------------
	reference_panglao_dir = panglaoDB_config[ config['ORGANISM'] ]
	reference_cellmarker_dir = cellmarker_config[ config['ORGANISM'] ]
	characterization_files = []
	cluster_mean_files, dataset_mean_files, possible_markers, marker_files_a= [], [], [], []
	characterization_path, figures_path = '', ''
	#-------------------------------------------

	# final_analysis.rules --------------------------------------
	final_integrated_seurat_object, final_analysed_integrated_seurat_object, final_resolution, final_directory = '', '', '', ''
	finals = []
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	final_files = [] #files to be passed to biggie

	if os.path.isfile(signicant_pc) == False:
		final_files = [soupX_list, doubletFinder_list, 'src/scripts/create_merged_dataset.R', merge_list, sigs]
	#-------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------
	# analyze_sc_object.rules-------------------------------------------------------------------------------------
	# the parameters for this rule are located in src/rules/analyze_sc_object.rules.
	# integrated_seurat_objects is for ALL or the specific integrated_seurat_object unique to Standard, RPCA, or SCT is added to the integrated_seurat_objects list
	# will pass integrated_seurat_objects as list to next rule (DGE plots).

	pca_path = '' #future path to analysis/PCA_*
	if os.path.isfile(signicant_pc):
		COMPONENTS_S = get_sig_pc(signicant_pc)
		PCA = 'PCA_' + str(COMPONENTS_S) #This is used/re-used

		out_parts_sc = []
		if STORAGE.lower() == 'h5':
			out_parts_sc = ['endpoints', PROJECT, 'analysis', PCA, 'H5']
		if STORAGE.lower() == 'rds':
			out_parts_sc = ['endpoints', PROJECT, 'analysis', PCA, 'RDS']
		string_path = 'data'
		for o in out_parts_sc:
			string_path = string_path + '/' + o
		pca_path = string_path.replace('/H5', '').replace('/RDS', '')

		if METHOD != 'ALL':
			integrated_seurat_object = string_path + '/' + PROJECT + '_' + METHOD
			integrated_seurat_object = storage(integrated_seurat_object, STORAGE)
			integrated_seurat_objects.append(integrated_seurat_object)
			if integrated_seurat_object not in integrated_seurat_objects:
				sc_objects.append(integrated_seurat_object)
			meta_file = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/' + PROJECT.lower() + '_' + METHOD + '_metatags.txt'
			sc_objects.append(meta_file)

		if METHOD == 'ALL':
			methods = ['RPCA', 'SCT', 'Standard']
			for m in methods:
				integrated_seurat_object = string_path + '/' + PROJECT + '_' + m
				integrated_seurat_object = storage(integrated_seurat_object, STORAGE)

				if integrated_seurat_object not in integrated_seurat_objects:
					integrated_seurat_objects.append(integrated_seurat_object)

				meta_file = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/' + PROJECT.lower() + '_' + m + '_metatags.txt'
				sc_objects.append(integrated_seurat_object)
				sc_objects.append(meta_file)
	#-------------------------------------------------------------------------------------

		#-------------------------------------------------------------------------------------
		if (os.path.exists(pca_path) == False): #
			final_files = sc_objects
		#-------------------------------------------------------------------------------------

		#-------------------------------------------------------------------------------------
		# DGE PLOT RULES ---------------------------------------------------------------------
		if os.path.exists(pca_path): #sigPC exists and os.path.isfile(signicant_pc)
			with open(resolution_file, 'w') as res_file:
				for s in sc_objects:
					if 'metatags' in s:
						if os.path.isfile(s):
							with open(s, 'r') as input_file:
								for line in input_file.readlines():
									if (('integrate' in line or '_snn_res' in line) and (line not in resolution_list)): #snn_res guarantees RNA vs integration
										line = line.replace('\n', '').lstrip().rstrip()
										if line not in resolution_list:
											res_file.write(line + '\n')
											resolution_list.append(line)
						if os.path.isfile(s) == False:
							print('There are no metatages...something is amuck.')

			if METHOD != 'ALL':
				for re in resolution_list:
					fname = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/figures/dge_plots/' + METHOD + '/' + PROJECT.lower() + '_initial_cluster_plots_' + re + '_' + METHOD + '.pdf'
					if fname not in images_files:
						images_files.append(fname)

					proportion_file = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dge_plots/' + METHOD + '/' + PROJECT.lower() + '_clusterProportions_' + re + '_' + METHOD + '.txt'
					if proportion_file not in proportion_files:
						proportion_files.append(proportion_file)
						
					proportion_image = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/figures/dge_plots/' + METHOD + '/' + PROJECT.lower() + '_clusterProportions_' + re + '_' + METHOD + '.pdf'
					if proportion_image not in proportion_images:
						proportion_images.append(proportion_image)

					markers = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dge_plots/' + METHOD + '/' + PROJECT.lower() + '_markers_' + re + '_' + METHOD + '.txt'
					if markers not in marker_files:
						marker_files.append(markers)

					markers100 = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dge_plots/' + METHOD + '/' + PROJECT.lower() + '_top100_markers_' + re + '_' + METHOD + '.txt'
					if markers100 not in marker_100files:
						marker_100files.append(markers100)

			if METHOD == 'ALL':
				methods = ['Standard', 'RPCA', 'SCT']
				for m in methods:
					for re in resolution_list:
						fname = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/figures/dge_plots/' + m + '/' + PROJECT.lower() + '_initial_cluster_plots_' + re + '_' + m + '.pdf'
						if fname not in images_files:
							images_files.append(fname)

						proportion_file = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dge_plots/' + m + '/' + PROJECT.lower() + '_clusterProportions_' + re + '_' + m + '.txt'
						if proportion_file not in proportion_files:
							proportion_files.append(proportion_file)

						proportion_image = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/figures/dge_plots/' + m + '/' + PROJECT.lower() + '_clusterProportions_' + re + '_' + m + '.pdf'
						if proportion_image not in proportion_images:
							proportion_images.append(proportion_image)

						markers = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dge_plots/' + m + '/' + PROJECT.lower() + '_markers_' + re + '_' + m + '.txt'
						markers100 = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dge_plots/' + m + '/' + PROJECT.lower() + '_top100_markers_' + re + '_' + m + '.txt'
						if markers not in marker_files:
							marker_files.append(markers)
						if markers100 not in marker_100files:
							marker_100files.append(markers100)

			dge_files = [images_files, marker_files, proportion_files, proportion_images]
			migsig_path = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/msig/'
		#-------------------------------------------------------------------------------------

			# DATASET CHARACTERIZATION ----------------------------------------------------------------
			figures_path = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/figures/possible_markers/'

			characterization_path = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dataset_characterization/'

			if METHOD != 'ALL':
				for res in resolution_list:
					rna_cluster_means = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dataset_characterization/' + PROJECT.lower() + '_meanGE_RNA_clusters_' + res + '_' + METHOD + '.txt'
					cluster_mean_files.append(rna_cluster_means)

					dataset_rna_mean = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dataset_characterization/' + PROJECT.lower() + '_meanGE_RNA_' + res + '_' + METHOD + '.txt'
					dataset_mean_files.append(dataset_rna_mean)

					possible_marker = characterization_path + PROJECT.lower() + '_possible_cell_types_' + res + '_' + METHOD + '_' + str(MARKER_THRESHOLD) + '.txt' 
					possible_markers.append(possible_marker) 

			if METHOD == 'ALL':
				methods = ['Standard', 'RPCA', 'SCT']
				for m in methods:
					for res in resolution_list:
						rna_cluster_means = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dataset_characterization/' + PROJECT.lower() + '_meanGE_RNA_clusters_' + res + '_' + m + '.txt'
						cluster_mean_files.append(rna_cluster_means)

						dataset_rna_mean = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dataset_characterization/' + PROJECT.lower() + '_meanGE_RNA_' + res + '_' + m + '.txt'
						dataset_mean_files.append(dataset_rna_mean)

						possible_marker = characterization_path + PROJECT.lower() + '_possible_cell_types_' + res + '_' + m + '_' + str(MARKER_THRESHOLD) + '.txt' 
						possible_markers.append(possible_marker) 

			characterization_files = [cluster_mean_files, dataset_mean_files, possible_markers, PROJECT.lower() + '_celltype_annotation_plot_dummy.txt', PROJECT.lower() + '_msig_annotation_plot_dummy.txt']
		#-------------------------------------------------------------------------------------

		#-------------------------------------------------------------------------------------
			# REPORT 
			#html_files = []


			#-------------------------------------------------------------------------------------
			if METHOD == 'Single':
				final_files = [dge_files, characterization_files]

			if METHOD != 'Single':
				# FINAL ANALYSIS ----------------------------------------------------------------
				#if CLUSTER_ANNOTATION_FILE == '' or CLUSTER_ANNOTATION_FILE is None:
				if CLUSTER_ANNOTATION_FILE == '' or CLUSTER_ANNOTATION_FILE is None or os.path.exists(CLUSTER_ANNOTATION_FILE)==False:
					final_files = [dge_files, characterization_files]
					print('------------------------------------------------')
					print('AUCHTUNG AUCHTUNG AUCHTUNG!!!!')
					print('You must supply \'configs/final_analysis_configs.yaml\' and enter values into \'configs/final_analysis_configs.yaml to finish the pipeline\'')
					print('------------------------------------------------')
				else:
					if USER_SUPPLIED_SEURAT_OBJECT.lower() == 'n':
						for r in resolution_list:
								if str(FINAL_RESOLUTION) in r:
									final_resolution = r
							
						for i in integrated_seurat_objects:
							if (FINAL_METHOD in i):
								final_integrated_seurat_object = i
								print(final_integrated_seurat_object)

								#if final_integrated_seurat_object.is_file() == False:
								if os.path.isfile(final_integrated_seurat_object) == False:
									print('There is something wrong with your final selection, and the final integrated object cannot be found/determined')
									print('Check your final resolution and final method in the configs/final_analysis_configs.yaml.')
									break

					if USER_SUPPLIED_SEURAT_OBJECT.lower() == 'y':
						final_integrated_seurat_object = PATH_USER_SUPPLIED_SEURAT_OBJECT
						if final_integrated_seurat_object.is_file() == False:
							print('\n\n')
							print('AUCHTUNG AUCHTUNG AUCHTUNG!!!!')
							print('You have have not supplied a final integrated seurat object.')
							print('Provide a path and name to the final seurat object in the configs/final_analysis_configs.yaml.')
							sys.exit()

					final_figures = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/figures/final_analysis/' + PROJECT.lower() + '_final_cluster_plots_' + final_resolution + '_' + FINAL_METHOD + '.pdf'
					final_analysed_integrated_seurat_object = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/' + STORAGE.upper() + '/' + PROJECT.lower() + '_' + FINAL_METHOD + '_' + str(FINAL_RESOLUTION) + '_annotated' 
					final_analysed_integrated_seurat_object = storage(final_analysed_integrated_seurat_object, STORAGE)
					final_directory = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA
					finals = [final_figures, final_analysed_integrated_seurat_object]

					final_files = [dge_files, characterization_files, finals]
			#-------------------------------------------------------------------------------------

	#--------------------TARGET RULES-----------------------------------
	include:
		"src/rules/soupX.rules"
	include:
		"src/rules/doubletFinder.rules"
	include:
		"src/rules/create_dynamic_script.rules"
	include:
		"src/rules/merge_samples.rules"
	include:
		"src/rules/get_sig_PC.rules"

	if os.path.isfile(signicant_pc):
		include:
			"src/rules/analyze_sc_object.rules"

		if os.path.exists(pca_path):
			include:
				"src/rules/dge_plots.rules"

			if METHOD == 'ALL':
				include:
					"src/rules/calculate_means_plot_markers_ALL.rules"
			if METHOD != 'ALL':
				include:
					"src/rules/calculate_means_plot_markers.rules"

			include:
				"src/rules/final_analysis.rules"
	#-------------------------------------------------------------------------------------

	#--------------------MESSAGES-----------------------------------
	onsuccess:
		print("The main controller pipeline completed with no errors.")
		shell("mail -s 'The main controller pipeline completed with no errors.' "+ config['contact']+" < {log}")

	onerror:
		print("The main controller pipeline did not complete without errors."),
		shell("mail -s 'The main controller pipeline did not complete without errors, check the logs and try again.' "+ config['contact']+" < {log}")
	#-------------------------------------------------------------------------------------


	#--------------------OUTPUT--------------------------------------
	# listed by rule
	#soupX_list, #soupX.rules
	#doubletFinder_list, #doubletFinder.rules
	#'src/scripts/create_merged_dataset.R', #create_dynamic_script.rules
	#merge_list, #merge_samples.rules
	#sigs, #get_sig_PC.rules
	#sc_objects, #analyze_sc_object.rules
	#dge_files, #dge_plots.rules
	#characterization_files #calculate_means_plot_markers.rules
	#finals #final_analysis.rules

	#--------------------RULES---------------------------------------
	rule biggie:
		input:
			final_files
	#-------------------------------------------------------------------------------------
