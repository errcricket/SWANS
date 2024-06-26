args = commandArgs(trailingOnly=TRUE)

project <- ''
method <- ''
lib_path <- ''
storage <- ''
compo <- ''
res_file <- ''
subdirectory1 <- '' #figures_path 
integrated_object <- ''
possible_file_directory <- ''
threshold <- 11
user_gene_file <- '' #may not exist

# test if there is at least 12 arguments: if not, return an error
if (length(args) < 11)
{ 
  stop('At least 11 arguments must be supplied.', call.=FALSE)
}
   
if (length(args)==11)
{  
   project = args[1]
   method = args[2]
   lib_path = args[3]
   storage = args[4]
   compo = args[5]
   res_file = args[6]
   subdirectory1 = args[7]
	integrated_object = args[8]
   possible_file_directory = args[9] 
	threshold = args[10]
   user_gene_file = args[11]
}

if (length(args)>11)
{  
   project = args[1]
   method = args[2]
   lib_path = args[3]
   storage = args[4]
   compo = args[5]
   res_file = args[6]
   subdirectory1 = args[7]
	integrated_object = args[8]
   possible_file_directory = args[9] 
	threshold = args[10]
   user_gene_file = args[11:length(args)]
}

print(res_file)
res_list <- read.table(res_file, header=FALSE)
res_list <- as.vector(res_list$V1)
print(res_list)

print(integrated_object)
print(possible_file_directory)
#--------------------------------------------------------------------
library(dplyr, lib.loc=lib_path)
library(Seurat, lib.loc=lib_path)
library(ggplot2, lib.loc=lib_path)
library(data.table, lib.loc=lib_path)
#--------------------------------------------------------------------

dir.create(subdirectory1, showWarnings=FALSE)

# IMPORT DATA
#--------------------------------------------------------------------
storage_retrieval <- function(storage, integrated_object)
{
   S <- ''
   if (tolower(storage) == 'h5')
   {
      print('Reading H5 object')
      S = read_h5(file=integrated_object, target.object='seurat')
   }

   if (tolower(storage) == 'rds')
   {
      print('Reading RDS object')
      S = readRDS(integrated_object)
   }

   return(S)
}
#--------------------------------------------------------------------


#--------------------------------------------------------------------
mark_cells <- function(seurat_object, possible_file_directory, res_list, met, directory, threshold)
{
	for (res in res_list)
	{
		print(res)
		Idents(object=seurat_object) <- res

		m_file <- paste(possible_file_directory, '/', project, '_possible_cell_types_', res, '_', met, '_', threshold, '.txt', sep='')
		if (file.exists(m_file) && file.size(m_file) != 0L)
		{
			print(m_file)
			marker_files <- read.table(m_file, header=FALSE) #this is a list of cell type file locations
			marker_file = as.list(marker_files$V1) 
				
			for (marker in marker_file)
			{
				base_name = basename(marker) #this gives me the name of the celltype file (Gastric_Cancer.txt)
				base_name = gsub('.txt', '', base_name)

				if (file.exists(marker) == TRUE)
				{
					print(marker)
					gene_set <- read.table(marker, header=FALSE)
					genes = as.vector(gene_set$V1) 

					print(genes)
	
					pdf(file=paste(directory, '/', tolower(project), '_pannotation_', res, '_', base_name, '_', met, '.pdf', sep=''), onefile=TRUE, width=11, height=8.5)
					try(print(DoHeatmap(subset(seurat_object), group.by=res, features=c(genes), size=3)))
					try(print(DotPlot(seurat_object, features=c(genes)) + RotatedAxis()))
					try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, features=c(genes))))
					try(print(RidgePlot(seurat_object, features = c(genes))))
	
					dev.off()
				}
			}
		}
	}
}
#--------------------------------------------------------------------

user_genes <- function(seurat_object, user_marker_file, res_list, met, directory)
{
   if (user_marker_file != '')
   {
		directory=paste(directory, '/user_defined_markers', sep='')
		dir.create(directory, showWarnings=FALSE)
	
		for (res in res_list)
		{
			print(res)
			Idents(object=seurat_object) <- res

			for (marker_file in user_marker_file)
			{
				if (file.exists(marker_file) == TRUE)
				{
					markers <- read.table(marker_file, header=FALSE) #this is a list of cell type file locations
					#genes = as.list(markers$V1) 
					genes = as.vector(markers$V1) 
				
					#gene_set <- read.table(marker, header=FALSE)

					print(genes)

					pdf(file=paste(directory, '/', tolower(project), '_user_genes_', res, '_', met, '.pdf', sep=''), onefile=TRUE, width=11, height=8.5)
					try(print(DoHeatmap(subset(seurat_object), group.by=res, features=c(genes), size=3)))
					try(print(DotPlot(seurat_object, features=c(genes)) + RotatedAxis()))
					try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, features=c(genes))))
				}
			}

		dev.off()

      }
   }
}
#--------------------------------------------------------------------

marker_files <- list.files(possible_file_directory)
marker_files <- subset(marker_files, subset = (grepl('possible', marker_files, fixed = TRUE) == TRUE))
print(marker_files)

# FUNCTION CALLS
#--------------------------------------------------------------------
if (method != 'ALL')
{
	subdirectory1=paste(subdirectory1, '/', method , sep='')
	dir.create(subdirectory1, showWarnings=FALSE)

	S <- storage_retrieval(storage, integrated_object)
	mark_cells(S, possible_file_directory, res_list, method, subdirectory1, threshold)
   if (file.exists(user_gene_file))
	{
		print('user has supplied a marker file')
		user_genes(S, user_gene_file, res_list, method, subdirectory1)
	}
}
