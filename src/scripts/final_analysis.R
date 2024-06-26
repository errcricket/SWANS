# Author:	E. Reichenberger
# Date:		8.10.2021

# Purpose: 	Create cluster images and look at conserved and DE markers

#(r message=FALSE)

args = commandArgs(trailingOnly=TRUE)
#print(args)
project <- ''
final_method <- ''
lib_path <- ''
storage <- ''
resolution <- ''
integrated_object <- ''   
celltype_assignment_file <- ''
base_directory <- ''
seurat_object_final <- ''

# test if there is at least 9 arguments: if not, return an error
if (length(args) < 9) 
{
  stop('At least nine arguments must be supplied.', call.=FALSE)
} 

if (length(args)==9) 
{
	project = args[1] 
	final_method = args[2] 
	lib_path = args[3]
	storage = args[4]
	resolution = args[5]
	integrated_object = args[6]
	celltype_assignment_file = args[7]
	base_directory = args[8]
	seurat_object_final = args[9]
}

library(RColorBrewer)
library(dplyr, lib.loc=lib_path)
library(ggrepel, lib.loc=lib_path)
library(stringr, lib.loc=lib_path)
library(Seurat, lib.loc=lib_path)
library(ggplot2, lib.loc=lib_path)
library(reshape2)
library(data.table, lib.loc=lib_path)

# CREATE DIRECTORIES
#--------------------------------------------------------------------
subdirectory2 <- ''

if (tolower(storage) == 'h5')
{
   library(hdf5r, lib=lib_path)
   library(RIOH5, lib=lib_path)
   library(SingleCellExperiment)
   subdirectory2=paste(base_directory, '/H5', sep='')
}

if (tolower(storage) == 'rds')
{
   subdirectory2=paste(base_directory, '/RDS', sep='')
}

subdirectory1=paste(base_directory, '/figures/final_analysis', sep='')
subdirectory3=paste(base_directory, '/tables/final_analysis', sep='')
dir.create(subdirectory1, showWarnings=FALSE)
dir.create(subdirectory3, showWarnings=FALSE)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# IMPORT DATA
#--------------------------------------------------------------------
storage_retrieval <- function(saved_object)
{
	S <- ''
	if (tolower(storage) == 'h5')
	{
		print('Reading H5 object')
		S = read_h5(file=saved_object, target.object='seurat')
	}

	if (tolower(storage) == 'rds')
	{
		print('Reading RDS object')
		S = readRDS(saved_object)
	}

	return(S)
}
#--------------------------------------------------------------------

save_object <- function(seurat_object, storage, file_name=seurat_object_final)
{
	print('saving ...')
	print(file_name)
   if (tolower(storage) == 'h5')
   {
      write_h5(data=seurat_object, object.type = 'seurat', file = file_name)
   }

   if (tolower(storage) == 'rds')
   {
      saveRDS(seurat_object, file = file_name)
   }
}

# Rename Clusters
#--------------------------------------------------------------------
rename_and_visualize <- function(seurat_object, celltype_file, met, s1=subdirectory1, res=resolution, storage=storage)
{
	print('renaming clusters...')
	Idents(object=seurat_object) <- res

	annotation <- read.table(celltype_file, header=TRUE, sep='\t')

	old.cluster.ids <- annotation$cluster
	new.cluster.ids <- annotation$celltype

	names(new.cluster.ids) <- levels(seurat_object)
	seurat_object <- RenameIdents(seurat_object, new.cluster.ids)

	seurat_object@meta.data$celltypes <- Idents(seurat_object)

	print('creating dimplot images...')

	f1 = paste(s1, '/', tolower(project), '_final_cluster_plots_', res, '_', met, '.pdf', sep='')
	pdf(file=f1, onefile=TRUE, width=11, height=8.5)

	print(DimPlot(seurat_object, reduction = 'umap', label=TRUE, group.by="celltypes")) #showing clusters with new names
	print(DimPlot(seurat_object, reduction = 'umap', label=TRUE, group.by="Phase")) #showing clusters with new names
	print(DimPlot(seurat_object, reduction = 'umap', label=TRUE)) #showing clusters with new names
	print(DimPlot(seurat_object, reduction = 'umap', label=TRUE, split.by='Experiment')) #showing clusters with new names by experiment
	print(DimPlot(seurat_object, reduction = 'umap', label=TRUE, split.by='Samples')) #showing clusters with new names by sample

	print(DimPlot(seurat_object, reduction = 'tsne', label=TRUE)) #showing clusters with new names
	print(DimPlot(seurat_object, reduction = 'tsne', label=TRUE, split.by='Experiment')) #showing clusters with new names by experiment
	print(DimPlot(seurat_object, reduction = 'tsne', label=TRUE, split.by='Samples')) #showing clusters with new names by sample

	dev.off()

	return(seurat_object)
}
#--------------------------------------------------------------------

# Find DE Genes in each cluster
#--------------------------------------------------------------------
dge <- function(seurat_object)
{
	samples <- length(unique(seurat_object@meta.data[["Samples"]]))

	if (samples > 1)
	{
		print('calculating difference in gene expression by cluster...')
		clusters <- unique(seurat_object@active.ident)
		
		for (c in clusters)
		{
			print(c)
			cluster <- subset(seurat_object, idents = c)
			cluster$celltype.exp <- paste(Idents(cluster), cluster$Experiment, sep = "_")
			cluster$celltype <- Idents(cluster)
			Idents(cluster) <- "celltype.exp"

			conditions <- unique(cluster$celltype.exp)

			#in the event there is more than 2 conditions
			for (i in conditions)
			{
				ci <- subset(x = cluster, subset = celltype.exp == i)
				cix <- nrow(ci@meta.data)

				for (j in conditions)
				{
				  cj <- subset(x = cluster, subset = celltype.exp == j)
				  cjx <- nrow(cj@meta.data)

				  if (cix > 3 && cjx > 3) #need at least 3 cells per cluster
				  {
						if (i != j)
						{
							table_name = paste(subdirectory3, '/', project, '_final_markers_', i, '_', j, '.txt', sep='')

							project.markers <- FindMarkers(cluster, ident.1 = i, ident.2 = j, verbose = TRUE, only.pos = FALSE)

							setDT(project.markers, keep.rownames = TRUE)[]
							names(project.markers)[1]<-paste('gene')
							write.table(project.markers, file=table_name, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
					 
							filteredSignificant <- subset(project.markers, p_val_adj <= 0.05)
							setDT(filteredSignificant, keep.rownames = TRUE)[]
							names(filteredSignificant)[1]<-paste('gene')

							sig_table_name = paste(subdirectory3, '/', project, '_significant_final_markers_', i, '_', j, '.txt', sep='')
							write.table(filteredSignificant, file=sig_table_name, sep='\t', quote=F, col.names=TRUE, row.names=FALSE, append=FALSE)
						}
					}
				}
			}
		}
	}
}
#--------------------------------------------------------------------

# FUNCTION CALLS ----------------------------------------------------

S <- storage_retrieval(integrated_object)
S <- rename_and_visualize(S, celltype_assignment_file, final_method, s1=subdirectory1, res=resolution, storage=storage)
save_object(S, storage)
dge(S)
