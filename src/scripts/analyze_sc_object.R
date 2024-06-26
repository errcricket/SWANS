# Author:	E. Reichenberger
# Date:		4.18.2021

# Purpose: 	Using the calculated sigPC to normalize (standard, RPCA, SCT) data, regress out %mito/cell cycling, and find clusters

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
	#top30 <- head(VariableFeatures(seurat_object), 30)
	top30 <- head(FindVariableFeatures(seurat_object), 30)

	plot1 <- VariableFeaturePlot(seurat_object)
	plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)

	f2 = paste(subdirectory1, '/', project, '_top30_variable_features_', method, '.pdf', sep='')
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

	if (met != 'SCT' && state == 'pre')
	{
		seurat_object <- ScaleData(seurat_object, verbose = FALSE)
	}

	print('running PCA')
	seurat_object <- RunPCA(seurat_object, features = c(s.genes, g2m.genes), npcs = compos, verbose = TRUE)
	print('making violin plots')
	print(VlnPlot(seurat_object, features = c('S.Score', 'G2M.Score', 'CC.Difference'), group.by = 'Phase', ncol = 3, pt.size = 0.1))
	print('making dimplots plots')
	#print(ggplot(seurat_object, aes(Phase)) + geom_bar())
	print(DimPlot(seurat_object)) 
	print('running PCA again')
	seurat_object <- RunPCA(seurat_object, features = features, npcs = compos, verbose = TRUE)
	print(DimPlot(seurat_object))
	
	dev.off()
      
	if (state == 'post')
	{  
		print('About to return seurat object')
		return(seurat_object)
	}
}
#--------------------------------------------------------------------

# Regression Scaling
#--------------------------------------------------------------------
scale_regression <- function(seurat_object)
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

	return(seurat_object)
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
		res_list <- as.numeric(res_list)
		print(res_list)
		seurat_object <- FindClusters(seurat_object, resolution = res_list)
	}

	else
	{
		print('There is something wrong with the resolutions')
	}

	return(seurat_object)
}
#--------------------------------------------------------------------

# RUN PCA-TSNE
#--------------------------------------------------------------------
transform_object_standard <- function(S, directory, res, compos)
{
	#SO.list <- SplitObject(S, split.by = 'Experiment')
	SO.list <- SplitObject(S, split.by = 'Samples')

	print('using standard analysis approach...')

	SO.list <- lapply(X = SO.list, FUN = function(x) {
		x <- NormalizeData(x)
		x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
		x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	})

	# for alternative cell cycling 
	for (i in 1:length(SO.list))
	{
		SO.list[[i]]@meta.data[['CC.Difference']] <- SO.list[[i]]@meta.data[['S.Score']] - SO.list[[i]]@meta.data[['G2M.Score']]
	}

	features <- SelectIntegrationFeatures(object.list = SO.list)

	SO.anchors <- FindIntegrationAnchors(object.list = SO.list, anchor.features = features)
	seurat_object <- IntegrateData(anchorset = SO.anchors)

	#variable_features(seurat_object)

	DefaultAssay(seurat_object) <- 'integrated'

	plot_cell_cycle(seurat_object, 'Standard', res, 'pre', features)
	seurat_object <- scale_regression(seurat_object) #regressions based on user input

	print('PCA')
	seurat_object <- plot_cell_cycle(seurat_object, 'Standard', res, 'post', features)
	#seurat_object <- get_clusters(seurat_object, compos, resc, res, 'Standard') #calling RunUMAP, RunTSNE, FindNeighbors, FindClusters
	seurat_object <- get_clusters(seurat_object, compos, res, 'Standard') #calling RunUMAP, RunTSNE, FindNeighbors, FindClusters

	filename = paste(subdirectory2, '/', project, '_Standard', sep='')
	save_function(seurat_object, storage, filename, 'Standard')
}

#--------------------------------------------------------------------
transform_object_rpca <- function(S, directory, res, compos)
{
	SO.list <- SplitObject(S, split.by = 'Samples')

	print('using RPCA analysis approach...')

	# normalize and identify variable features for each dataset independently
	SO.list <- lapply(X = SO.list, FUN = function(x) {
		x <- NormalizeData(x)
		x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
		x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	})

	for (i in 1:length(SO.list))
	{
		SO.list[[i]]@meta.data[['CC.Difference']] <- SO.list[[i]]@meta.data[['S.Score']] - SO.list[[i]]@meta.data[['G2M.Score']]
	}

	# select repeatedly variable features across datasets for integration run PCA on each dataset using these features
	features <- SelectIntegrationFeatures(object.list = SO.list)

	#------------------------------------------------------------------
	SO.list <- lapply(X = SO.list, FUN = function(x) {
		x <- ScaleData(x, features = features, verbose = FALSE)
		x <- RunPCA(x, features = features, verbose = FALSE)
	})

	#k.anchor can be increased to improve alignment strength (5 is default)
	SO.anchors <- FindIntegrationAnchors(object.list = SO.list, anchor.features = features, reduction = 'rpca', k.anchor = 15)
	seurat_object <- IntegrateData(anchorset = SO.anchors)

	# unmodified data still resides in the 'RNA' assay
	DefaultAssay(seurat_object) <- 'integrated'

	print('plot')#
	plot_cell_cycle(seurat_object, 'RPCA', res, 'pre', features)
	seurat_object <- scale_regression(seurat_object) #regressions based on user input

	seurat_object <- plot_cell_cycle(seurat_object, 'RPCA', res, 'post', features)
	seurat_object <- get_clusters(seurat_object, compos, res, 'RPCA') #calling RunUMAP, RunTSNE, FindNeighbors, FindClusters

	filename = paste(subdirectory2, '/', project, '_RPCA', sep='')
	save_function(seurat_object, storage, filename, 'RPCA')
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
transform_object_sct <- function(S, directory, res, compos)
{
	print('using SCT analysis approach...')
	SO.list <- SplitObject(S, split.by = 'Samples')

	library(glmGamPoi, lib.loc=lib_path)
	comp <- compos + 10 #number of sig PC slightly higher for SCT than Standard approach. Adding 10....

	SO.list <- lapply(X = SO.list, FUN = function(x) {
		x <- NormalizeData(x)
		x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
		x <- SCTransform(x, verbose=TRUE, method = 'glmGamPoi', assay='RNA')
	})

	for (i in 1:length(SO.list))
	{
		SO.list[[i]]@meta.data[['CC.Difference']] <- SO.list[[i]]@meta.data[['S.Score']] - SO.list[[i]]@meta.data[['G2M.Score']]
	}

	# this section plots cc scores for preservation not keep seurat object
	#------------------------------------------------------------------
	features <- SelectIntegrationFeatures(object.list = SO.list, nfeatures = 3000)
	SO.list <- PrepSCTIntegration(object.list = SO.list, anchor.features = features)
	SO.list <- lapply(X = SO.list, FUN = RunPCA, features = features)

	SO.anchors <- FindIntegrationAnchors(object.list = SO.list, normalization.method = 'SCT', reduction='rpca', anchor.features = features, dims=1:compos, verbose=TRUE, k.anchor=15)
	seurat_object <- IntegrateData(anchorset = SO.anchors, normalization.method = 'SCT', dims=1:comp)

	DefaultAssay(seurat_object) <- 'integrated'
	#variable_features(SO.integrated)

	plot_cell_cycle(seurat_object, 'SCT', res, 'pre', features)
	seurat_object <- scale_regression(seurat_object) #regressions based on user input

	seurat_object <- plot_cell_cycle(seurat_object, 'SCT', res, 'post', features)
	seurat_object <- get_clusters(seurat_object, compos, res, 'SCT') #calling RunUMAP, RunTSNE, FindNeighbors, FindClusters

	filename = paste(subdirectory2, '/', project, '_SCT', sep='')
	save_function(seurat_object, storage, filename, 'SCT')
}
#--------------------------------------------------------------------

# FUNCTION CALLS
#--------------------------------------------------------------------
S <- storage_function(storage)

if (method == 'Standard')
{
	transform_object_standard(S, subdirectory2, res, compos=compo)
}

if (method == 'RPCA')
{
	transform_object_rpca(S, subdirectory2, res, compos=compo)
}

if (method == 'SCT')
{
	transform_object_sct(S, subdirectory2, res, compos=compo)
}

if (method == 'ALL')
{
	transform_object_rpca(S, subdirectory2, res, compos=compo)
	transform_object_standard(S, subdirectory2, res, compos=compo)
	transform_object_sct(S, subdirectory2, res, compos=compo)
}
