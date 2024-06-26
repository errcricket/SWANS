# Author:	E. Reichenberger
# Date:		8.10.2021

# Purpose: 	Create cluster images and look at conserved and DE markers

#(r message=FALSE)
#[1] "data/endpoints/schaefer/analysis/PCA_19/figures/dge_plots/Standard/schaefer_initial_cluster_plots_integrated_snn_res.0.5_Standardschaefer.pdf"


args = commandArgs(trailingOnly=TRUE)
compo <- ''
project <- ''
method <- ''
lib_path <- ''
storage <- ''
res_file <- ''
integrated_object <- ''   #this will now be a list of objects or can i call each one individually
conserved_genes <- ''   #this will now be a list of objects or can i call each one individually

# test if there is at least 7 arguments: if not, return an error
if (length(args) < 8) 
{
  stop('At least eight arguments must be supplied.', call.=FALSE)
} 

if (length(args)==8) 
{
	project = args[1] 
	method = args[2] 
	lib_path = args[3]
	storage = args[4]
	compo = args[5]
	res_file = args[6]
	conserved_genes = args[7]
	integrated_object = args[8]
}

if (length(args)>8) 
{
	project = args[1] 
	method = args[2] 
	lib_path = args[3]
	storage = args[4]
	compo = args[5]
	res_file = args[6]
	conserved_genes = args[7]
	integrated_object = args[8:length(args)]
}

res_list <- read.table(res_file, header=FALSE)
res_list <- as.vector(res_list$V1)
#print(integrated_object)
#print(res_list)

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

subdirectory1=paste(base_directory, '/figures/dge_plots', sep='')
subdirectory3=paste(base_directory, '/tables/dge_plots', sep='')
dir.create(subdirectory1, showWarnings=FALSE)
dir.create(subdirectory2, showWarnings=FALSE)
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

# Visualize clusters 
#--------------------------------------------------------------------
integrated_visualization <- function(seurat_object, res_list, sample, met, s1=subdirectory1)
{
	m_subdirectory1=paste(s1, '/', met, sep='')
	dir.create(m_subdirectory1, showWarnings=FALSE)

	for (res in res_list)
	{
		print(res)
		Idents(object=seurat_object) <- res
		
		f1 = paste(m_subdirectory1, '/', tolower(project), '_initial_cluster_plots_', res, '_', met, sample, '.pdf', sep='')
		print(f1)
		pdf(file=f1, onefile=TRUE, width=11, height=8.5)
		p1 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'Experiment')
		p2 <- DimPlot(seurat_object, reduction = 'umap', label = TRUE, repel = TRUE)
		print(p1 + p2)

		print(DimPlot(seurat_object, reduction = 'umap', group.by = 'orig.ident'))
		print(DimPlot(seurat_object, reduction = 'umap', split.by = 'orig.ident'))

		x1 <- DimPlot(seurat_object, reduction = 'tsne', group.by = 'Experiment')
		x2 <- DimPlot(seurat_object, reduction = 'tsne', label = TRUE, repel = TRUE)
		print(x1 + x2)

		print(DimPlot(seurat_object, reduction = 'umap', split.by = 'Experiment'))
		print(DimPlot(seurat_object, reduction = 'tsne', split.by = 'Experiment'))

		print(DimPlot(seurat_object, reduction = 'umap', label=TRUE, repel=TRUE))
		print(DimPlot(seurat_object, reduction = 'tsne', label=TRUE, repel=TRUE))

		print(DimPlot(seurat_object, reduction = 'umap', group.by = 'Phase'))
		print(DimPlot(seurat_object, reduction = 'tsne', group.by = 'Phase'))
		dev.off()
	}
}
#--------------------------------------------------------------------

#Visualize proportions
#--------------------------------------------------------------------
proportions <- function(seurat_object, res_list, sample, met, s1=subdirectory1, s3=subdirectory3)
{
	m_subdirectory1=paste(s1, '/', met, sep='')
	m_subdirectory3=paste(s3, '/', met, sep='')
	dir.create(m_subdirectory3, showWarnings=FALSE)

	for (res in res_list)
	{
		cluster_count <- unique(seurat_object@meta.data[[res]])
		print(cluster_count)

		print(res)
		Idents(object=seurat_object) <- res

		# use double brackets to handle res not being a colunmn name (e.g., $res)
		number_perCluster <- table(seurat_object@meta.data$Experiment, seurat_object@meta.data[[res]])

		table_name = paste(m_subdirectory3, '/', tolower(project), '_clusterProportions_', res, '_', met, sample, '.txt', sep='')
		write.table(number_perCluster, file=table_name, sep='\t', quote=F, col.names=TRUE, row.names=FALSE, append=FALSE)

		#will need to melt table
		df <- melt(number_perCluster)
		names(df)[1] <- 'Sample'
		names(df)[2] <- 'Cluster'

		df$Cluster <- as.factor(df$Cluster) 

		file_name = paste(m_subdirectory1, '/', tolower(project), '_clusterProportions_', res, '_', met, sample, '.pdf', sep='')
		pdf(file=file_name, onefile=TRUE, width=11, height=8.5)

		print(ggplot(df, aes(x = Sample, y=value, fill = Cluster)) +
		geom_bar(stat = 'identity', position='fill') + 
		ylab('Cluster Proportion') + xlab('Sample') + guides(fill=guide_legend(title='Cluster')))
		dev.off()
	}
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
get_markers  <- function(seurat_object, res_list, met, s3=subdirectory3)
{
	print('finding differentially expressed genes...')

	m_subdirectory3=paste(s3, '/', met, sep='')

	for (res in res_list)
	{
		print(res)
		Idents(object=seurat_object) <- res

		project.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

		f_name = paste(m_subdirectory3, '/', tolower(project), '_markers_', res, '_', met, '.txt', sep='')
		project.markers %>% group_by(cluster) 
		write.table(project.markers, file=f_name, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

		project.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
		top100 <- project.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

		f_name1 = paste(m_subdirectory3, '/', tolower(project), '_top100_markers_', res, '_', met, '.txt', sep='')
		write.table(top100, file=f_name1, sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
	}
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
get_conserved_markers  <- function(seurat_object, res_list, met, s3=subdirectory3)
{
	samples <- unique(seurat_object@meta.data[['orig.ident']])

	if (length(samples) > 1)
	{
		m_subdirectory3=paste(s3, '/', met, sep='')
		m_subdirectory3=paste(m_subdirectory3, '/conserved_markers', sep='')
		dir.create(m_subdirectory3, showWarnings=FALSE)

		print('finding differentially (conserved) expressed genes...')

		for (res in res_list)
		{
			print(res)
			Idents(object=seurat_object) <- res

			cluster_count <- unique(seurat_object@meta.data[[res]])

			for (count in cluster_count)
			{
				count <- as.numeric(count)
				print(count)
				project.conserved.markers <- FindConservedMarkers(seurat_object, ident.1=count, ident.2=NULL, grouping.var='Experiment', verbose=TRUE)

				file_name <- paste(m_subdirectory3, '/', tolower(project), '_conservedMarkers_cluster_', count, '_', res, '_', met, '.txt', sep='')
				write.table(project.conserved.markers, file=file_name, sep='\t', quote=F, col.names=TRUE, row.names=TRUE, append=FALSE)
			}
		}
	}
}
#--------------------------------------------------------------------

# FUNCTION CALLS ----------------------------------------------------

if (method != 'ALL')
{
	S <- storage_retrieval(integrated_object)
	samples <- unique(S@meta.data[['orig.ident']])
	if (samples == project)
	{
		samples <- unique(S@meta.data[['Samples']])
	}
	print(samples)
	integrated_visualization(S, res_list, '', method) #on combined dataset

	for (s in samples)
	{
		print(s)
		x <- subset(S, subset = orig.ident == s)
		integrated_visualization(x, res_list, s, method) #individual samples
		proportions(x, res_list, s, method)
	}

	proportions(S, res_list, '', method)

	get_markers(S, res_list, method)

	if (conserved_genes == 'y')
	{
		get_conserved_markers(S, res_list, met)
	}
}

if (method == 'ALL')
{
	Standard_file <- integrated_object[1]
	RPCA_file <- integrated_object[2]
	SCT_file <- integrated_object[3]

	trans_files <- c(Standard_file, RPCA_file, SCT_file)
	
	for (t in trans_files)
	{
		print(t)
		S <- '' #each loop needs to be fresh
		S <- storage_retrieval(t)
		print(S)
		samples <- unique(S@meta.data[['orig.ident']])
		if (samples == project)
		{
			samples <- unique(S@meta.data[['Samples']])
		}
		print(samples)

		met = ''
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

		integrated_visualization(S, res_list, '', met)

		for (s in samples)
		{
			print(s)
			x <- subset(S, subset = orig.ident == s)
			integrated_visualization(x, res_list, s, met) #individual samples
			proportions(x, res_list, s, met)
		}

		proportions(S, res_list, '', met)

		get_markers(S, res_list, met)

		if (conserved_genes == 'y')
		{
			get_conserved_markers(S, res_list, met)
		}
	}
}
