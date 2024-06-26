#Author:	E. Reichenberger
# Date:		8.10.2021

# Purpose: 	Create cluster images and look at conserved and DE markers

#(r message=FALSE)

args = commandArgs(trailingOnly=TRUE)
compo <- ''
project <- ''
method <- ''
lib_path <- ''
storage <- ''
integrated_object <<- ''
res_file <- ''

if (length(args) < 7)
	{ 
	  stop('At least seven arguments must be supplied.', call.=FALSE)
	}
   
if (length(args)==7)
	{
		project = args[1]
		method = args[2]
		lib_path = args[3]
		storage = args[4]
		compo = args[5]
		res_file = args[6]
		integrated_object = args[7]
	}

if (length(args)>7)
	{
		project = args[1]
		method = args[2]
		lib_path = args[3]
		storage = args[4]
		compo = args[5]
		res_file = args[6]
		integrated_object = args[7:length(args)]
	}

res_list <- read.table(res_file, header=FALSE)
res_list <- as.vector(res_list$V1)

library(dplyr, lib.loc=lib_path)
library(Seurat, lib.loc=lib_path)
library(data.table, lib.loc=lib_path)

# CREATE DIRECTORIES
#--------------------------------------------------------------------
base_directory=paste('data/endpoints/', project, '/analysis', sep='')
dir.create(base_directory, showWarnings=FALSE)
base_directory=paste(base_directory, '/PCA_', compo, sep='')
dir.create(base_directory, showWarnings=FALSE)

subdirectory2 <- ''

if (tolower(storage) == 'h5')
{
   library(hdf5r, lib=lib_path)
   library(RIOH5, lib=lib_path)
   library(SingleCellExperiment)
	#library(SingleCellExperiment, lib=lib_path)
   subdirectory2=paste(base_directory, '/H5', sep='')
}

if (tolower(storage) == 'rds')
{
   subdirectory2=paste(base_directory, '/RDS', sep='')
}

subdirectory1=paste(base_directory, '/figures/dataset_characterization', sep='')
subdirectory3=paste(base_directory, '/tables/dataset_characterization', sep='')
dir.create(subdirectory1, showWarnings=FALSE)
dir.create(subdirectory2, showWarnings=FALSE)
dir.create(subdirectory3, showWarnings=FALSE)
#--------------------------------------------------------------------

# IMPORT DATA
#--------------------------------------------------------------------
storage_retrieval <- function(storage, integrated_seurat_object)
{
	S <- ''
	if (tolower(storage) == 'h5')
	{
		print('Reading H5 object')
		S = read_h5(file=integrated_seurat_object, target.object='seurat')
	}

	if (tolower(storage) == 'rds')
	{
		print('Reading RDS object')
		S = readRDS(integrated_seurat_object)
	}

	return(S)
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
calculate_averages  <- function(seurat_object, res_list, met)
{
	print('calculating average gene expression by cluster...')

	#appending to file, delete to avoid run-ons...
	#rna_dataset_characterization = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/dataset_characterization/' + PROJECT.lower() + '_meanGE_RNA_clusters_' + res + '_' + METHOD + '.txt'

	avg_files <- paste(subdirectory3, '/', tolower(project), '_meta_meanGE_filenames_', met, '.txt', sep='')
	#avg_files <- paste(subdirectory3, '/', tolower(project), '_meanGE_RNA_clusters_', res, '_', met, '.txt', sep='')
	if (file.exists(avg_files)) 
	{
	  file.remove(avg_files)
	}

	for (res in res_list)
	{
		print(res)
		ae <- AverageExpression(seurat_object, group.by=res) 

		f_name_rna = paste(subdirectory3, '/', tolower(project), '_meanGE_RNA_clusters_', res, '_', met, '.txt', sep='')
		#write.table(f_name_rna, file=avg_files, sep='\t', append=TRUE, quote=F, col.names=FALSE, row.names=FALSE)
		write.table(ae$RNA, file=f_name_rna, sep='\t', quote=F, col.names=TRUE, row.names=TRUE)
	}
}
#--------------------------------------------------------------------

if (method != 'ALL')
{
	S <- storage_retrieval(storage, integrated_object)
	calculate_averages(S, res_list, method)
}

if (method == 'ALL')
{
   Standard_file <- integrated_object[1]
   RPCA_file <- integrated_object[2]
   SCT_file <- integrated_object[3]

   trans_files <- c(Standard_file, RPCA_file, SCT_file)

   for (t in trans_files)
   {
		met = ''
      print(t)
      S <- '' #each loop needs to be fresh
		S <- storage_retrieval(storage, t)

		if ( grepl('Standard', t, fixed = TRUE) == TRUE)
		{
			met = 'Standard'
		}

		if ( grepl('RPCA', t, fixed = TRUE) == TRUE)
		{
			met = 'RPCA'
		}

		if ( grepl('SCT', t, fixed = TRUE) == TRUE)
		{
			met = 'SCT'
		}

		calculate_averages(S, res_list, met)
	}
}
