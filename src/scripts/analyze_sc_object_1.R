# Author:	E. Reichenberger
# Date:		4.18.2021

# Purpose: 	Using the calculated sigPC to normalize (standard, Single, SCT) data, regress out %mito/cell cycling, and find clusters

#(r message=FALSE)

args = commandArgs(trailingOnly=TRUE)
compo <- ''
project <<- ''
method <- ''
lib_path <- ''
storage <- ''
merged_file <- ''
res <<- ''
cell_cycling <<- ''
cc_method <<- ''
organism <<- ''
mito_regression <<- ''

# test if there is at least 11 arguments: if not, return an error
if (length(args) < 11) 
{
  stop('At least 11 arguments must be supplied.', call.=FALSE)
} 

if (length(args)==11) 
{
	compo = args[1]
	project = args[2] 
	method = args[3] 
	lib_path = args[4]
	mito_regression = args[5]
	storage = args[6]
	cell_cycling = args[7]
	cc_method = args[8]
	organism = args[9]
	merged_file = args[10]
	res = args[11]
}

if (length(args)>11) 
{
	compo = args[1]
	project = args[2] 
	method = args[3] 
	lib_path = args[4]
	mito_regression = args[5]
	storage = args[6]
	cell_cycling = args[7]
	cc_method = args[8]
	organism = args[9]
	merged_file = args[10]
	res = args[11:length(args)]
}

compo = as.numeric(compo)
resc <- ''


library(dplyr, lib.loc=lib_path)
library(ggrepel, lib.loc=lib_path)
library(Seurat, lib.loc=lib_path)
library(patchwork, lib.loc=lib_path)
library(sctransform, lib.loc=lib_path)
library(data.table, lib.loc=lib_path)
#library(dplyr)

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
   subdirectory2=paste(base_directory, '/H5', sep='')
}

if (tolower(storage) == 'rds')
{
   subdirectory2=paste(base_directory, '/RDS', sep='')
}

subdirectory1=paste(base_directory, '/figures', sep='')
subdirectory3=paste(base_directory, '/tables', sep='')
dir.create(subdirectory1, showWarnings=FALSE)
dir.create(subdirectory2, showWarnings=FALSE)
dir.create(subdirectory3, showWarnings=FALSE)
#--------------------------------------------------------------------

# CYCLING GENE --------------------------------------------------------------------------
#Cell cycle genes for human included with seurat but must convert from seurat's human markers to mouse IDs.
  #https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

convertHumanGeneList <- function(x){
	suppressMessages(library('biomaRt', lib.loc=lib_path))
	human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl', host = 'https://dec2021.archive.ensembl.org/')
	mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl', host = 'https://dec2021.archive.ensembl.org/')
	genesV2 = getLDS(attributes = c('hgnc_symbol'), filters = 'hgnc_symbol', values = x , mart = human, attributesL = c('mgi_symbol'), martL = mouse, uniqueRows=T)
	humanx <- unique(genesV2[, 2])

	return(humanx)
}

s.genes <<- cc.genes$s.genes
g2m.genes <<- cc.genes$g2m.genes

if (tolower(organism) == 'mouse')
{
	s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
	g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
}
#--------------------------------------------------------------------------------------

# IMPORT DATA
#--------------------------------------------------------------------
storage_function <- function(storage)
{
	S <- ''

	if (tolower(storage) == 'h5')
	{
		print('H5')
		S = read_h5(file=merged_file, target.object='seurat')
	}

	if (tolower(storage) == 'rds')
	{
		print('RDS')
		S = readRDS(merged_file)
	}

	return(S)
}
#--------------------------------------------------------------------

# SAVE FILE
save_function <- function(seurat_object, storage, filename, met)
{
	if (tolower(storage) == 'h5')
	{
		print('saving object as H5')
		filename <- paste(filename, '.H5', sep='')
		write_h5(data = seurat_object, object.type = 'seurat', file = filename)
	}

	if (tolower(storage) == 'rds')
	{
		print('saving object as RDS')
		filename <- paste(filename, '.RDS', sep='')
		saveRDS(seurat_object, file=filename)
	}

	object_names <- names(seurat_object@meta.data)
	base_filename <- paste(subdirectory3, '/', tolower(project), '_', met, '_metatags.txt', sep='')
	write.table(object_names, file=base_filename, sep='\t', quote=F, col.names=FALSE, row.names=FALSE)
}

variable_features <- function(seurat_object) #could add if statement to calling this (e.g., if sample# == 1)
{
	top20 <- head(FindVariableFeatures(seurat_object), 20)

	plot1 <- VariableFeaturePlot(seurat_object)
	plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

	f2 = paste(subdirectory1, '/', project, '_top20_variable_features_', method, '.pdf', sep='')
	pdf(file=f2, onefile=TRUE, width=11, height=8.5)
	print(plot1 + plot2)
	dev.off()
}

# PLOT CELL CYCLING SCORES
#--------------------------------------------------------------------
plot_cell_cycle <- function(seurat_object, met, res, state, features, project_name=project, directory=subdirectory1, compos=compo)
{
	print('making cell cycling images...')

	pdf(file=paste(directory, '/', project_name, '_cell_cycling_', met, '_', state, '.pdf', sep=''), onefile=TRUE, width=11, height=8.5)

	if (state == 'pre')
	{
		print('scaling...')
		seurat_object <- ScaleData(seurat_object, verbose = FALSE)

		print('running PCA')
		seurat_object <- RunPCA(seurat_object, features = c(s.genes, g2m.genes), npcs = compos, verbose = TRUE)
		print('making violin plots')
		print(VlnPlot(seurat_object, features = c('S.Score', 'G2M.Score', 'CC.Difference'), group.by = 'Phase', ncol = 3, pt.size = 0.1))
		print('making dimplots plots')
		#print(ggplot(seurat_object, aes(Phase)) + geom_bar())
		print(DimPlot(seurat_object)) 
		
		dev.off()
	}
      
	if (state == 'post')
	{  
		if (cell_cycling == 'y')
		{  
			print('accounting for cell cycles')
			
			if (cc_method == 'standard')
			{  
				if (mito_regression == 'n')
				{  
					seurat_object <- ScaleData(seurat_object, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(seurat_object))
				}
				
				if (mito_regression == 'y')
				{  
					seurat_object <- ScaleData(seurat_object, vars.to.regress = c('percent.mito', 'S.Score', 'G2M.Score'), features = rownames(seurat_object))
				}
			}
			
			if (cc_method == 'alternative')
			{  
				if (mito_regression == 'n')
				{  
					seurat_object <- ScaleData(seurat_object, vars.to.regress = 'CC.Difference', features = rownames(seurat_object))
				}
				
				if (mito_regression == 'y')
				{  
					seurat_object <- ScaleData(seurat_object, vars.to.regress = c('percent.mito', 'CC.Difference'), features = rownames(seurat_object))
				}
			}
		}

		else
		{  
			if (mito_regression == 'y')
			{  
				seurat_object <- ScaleData(seurat_object, vars.to.regress='percent.mito', verbose = FALSE)
			}

			if (mito_regression == 'n')
			{
				seurat_object <- ScaleData(seurat_object, verbose = FALSE)
			}
		}

		seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
		print('About to return seurat object')
		return(seurat_object)
	}
}
#--------------------------------------------------------------------

# Get clustering info (same for all approaches)
#--------------------------------------------------------------------
#get_clusters <- function(seurat_object, compos, resc, res, met)
get_clusters <- function(seurat_object, compos, res, met)
{
	print('Running UMAP, TSNE and finding Neighbors and Clusters')
	print('Running UMAP')
	seurat_object <- RunUMAP(seurat_object, reduction = 'pca', dims = 1:compos)
	print('Running TSNE')
	seurat_object <- RunTSNE(seurat_object, reduction = 'pca', dims = 1:compos)
	print('Running FindNeighbors')
	seurat_object <- FindNeighbors(seurat_object, reduction = 'pca', dims = 1:compos)

	print('Running FindClusters')
	if (res == '')
	{  
		print('res = 0.8')
		seurat_object <- FindClusters(seurat_object, resolution=0.8)
	}

	else if (res == 'I')
	{  
		print('res = 0.5, 0.8, 1.2')
		seurat_object <- FindClusters(seurat_object, resolution = c(0.5, 0.8, 1.2))
	}

	# single value
	else if ( (sum(lengths(regmatches(res, gregexpr("\\.", res)))) == 1) || (length(res) == 1) )
	{
		print('single user-supplied resolution... ')
		print('res = ')
		res <- as.numeric(res)
		print(res)
		seurat_object <- FindClusters(seurat_object, resolution = res)
	}

	# multiple values -- how is this different
	#else
	else if (sum(lengths(regmatches(res, gregexpr("\\.", res)))) > 1)
	{  
		res_list <- list()
		print('multiple user-supplied resolutions... ')
		for (r in res)
		{
			r <- as.numeric(r)
			print(r)
			res_list = c(res_list, r)
		}
		seurat_object <- FindClusters(seurat_object, resolution = res_list)
	}

	else
	{
		print('There is something wrong with the resolutions')
	}

	return(seurat_object)
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
transform_object <- function(seurat_object, directory, res, compos)
{
	seurat_object <- NormalizeData(seurat_object)
	seurat_object <- FindVariableFeatures(seurat_object, selection.method = 'vst', nfeatures = 2000)
	seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	seurat_object@meta.data[['CC.Difference']] <- seurat_object@meta.data[['S.Score']] - seurat_object@meta.data[['G2M.Score']]

	print(names(seurat_object@meta.data))

	#variable_features(seurat_object)

	print('plot')#
	plot_cell_cycle(seurat_object, 'Single', res, 'pre', features)
	print('finished plotting')

	print('plotting 2')#
	seurat_object <- plot_cell_cycle(seurat_object, 'Single', res, 'post', features)
	print('finished plotting round 2')
	seurat_object <- get_clusters(seurat_object, compos, res, 'Single') #calling RunUMAP, RunTSNE, FindNeighbors, FindClusters

	filename = paste(subdirectory2, '/', project, '_Single', sep='')
	save_function(seurat_object, storage, filename, 'Single')
}

#--------------------------------------------------------------------

# FUNCTION CALLS
#--------------------------------------------------------------------
S <- storage_function(storage)
transform_object(S, subdirectory2, res, compo)
