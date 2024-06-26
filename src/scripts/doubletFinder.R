# Author:		E. Reichenberger
# Date:			10.31.2020
# Modifified: 	3.11.2022

# Purpose: Remove doublets from samples. The script should be called on each individual sample after the samples have been processed with soupX. 

# Notes: Doublets are identified in the 10X data and the doublet barcodes are removed from the soupX data. 

lib_path <- ''
sample <- ''
project <- ''
components <- ''
soup_path <- ''
doubletFinder_path <- ''
organism <- ''
sequencing <- ''
data_type <- ''

args = commandArgs(trailingOnly=TRUE)
#print(length(args))
#print(args)

# test if there is at least 9 arguments: if not, return an error
if (length(args) < 9) {
  stop('At least 9 arguments must be supplied.', call.=FALSE)
} else if (length(args)==9) {
   lib_path = args[1]
   sample = args[2] 
   project = args[3] #project name
   components = args[4]
   soup_path = args[5]
   doubletFinder_path = args[6]
   organism = args[7]
   sequencing = args[8]
	data_type = args[9]
}

suppressMessages(library('ggplot2', lib.loc=lib_path))
suppressMessages(library('DoubletFinder', lib.loc=lib_path))
suppressMessages(library('Seurat', lib.loc=lib_path))
suppressMessages(library('DropletUtils', lib.loc=lib_path))

# ASSIGN VARIABLES
pN <- 0.25
doublet_predictor <- 0 #NOTE this number will be based on # of recovered cells (cellranger, filtered)

# CREATE DIRECTORIES
#--------------------------------------------------------------------
# Create figure and tables folders under sample folder
base_directory=gsub('soupX', '', soup_path)

dir.create(base_directory, showWarnings=FALSE)
subdirectory1=paste(base_directory, '/figures', sep='')
subdirectory2=paste(base_directory, '/tables', sep='')

dir.create(subdirectory1, showWarnings=FALSE)
dir.create(subdirectory2, showWarnings=FALSE)
soup_path = paste(soup_path, '/', sep='')
#print(soup_path)
#--------------------------------------------------------------------

# READ IN SOUPX MATRIX FILES & CREATE SEURAT OBJECT
#--------------------------------------------------------------------
tenX_to_seuratOB <- function(soup_path, organism, components)
{
	print('listing files in 10X directory as a sanity check...')
	#this will need to change if h5 is used.
	X10_path=gsub('soupX/', '10X/outs/filtered_feature_bc_matrix/', soup_path)
	print(list.files(path=X10_path))
	sample.data <- Read10X(data.dir = X10_path)
	so <- CreateSeuratObject(counts = sample.data, project = project, min.cells = 0, min.features = 0)

	return(so)
}
#--------------------------------------------------------------------

#GET MULTIPLET RATE FROM NUMBER RECOVERED CELLS
#--------------------------------------------------------------------
determine_multiplet_rate <- function(seurat_object)
{
	print('determing multiplet rate')
	cell_count <- length(seurat_object@meta.data[["orig.ident"]])
	print('Number of recovered cells...')
	print(cell_count)

	if (cell_count > 0 && cell_count <= 740)
		{doublet_predictor <- 0.004}
	if (cell_count > 740 && cell_count <= 1400)
		{doublet_predictor <- 0.008}
	if (cell_count > 1500 && cell_count <= 2400)
		{doublet_predictor <- 0.016}
	if (cell_count > 2500 && cell_count <= 3400)
		{doublet_predictor <- 0.023}
	if (cell_count > 3400 && cell_count <= 4400)
		{doublet_predictor <- 0.031}
	if (cell_count > 4400 && cell_count <= 5400)
		{doublet_predictor <- 0.039}
	if (cell_count > 5400 && cell_count <= 6400)
		{doublet_predictor <- 0.046}
	if (cell_count > 6400 && cell_count <= 7400)
		{doublet_predictor <- 0.054}
	if (cell_count > 7400 && cell_count <= 8400)
		{doublet_predictor <- 0.061}
	if (cell_count > 8400 && cell_count <= 9400)
		{doublet_predictor <- 0.069}
	if (cell_count > 9400 && cell_count <= 10400)
		{doublet_predictor <- 0.076}

	return(doublet_predictor)
}
#--------------------------------------------------------------------

# FILTER SEURAT OBJECT
#--------------------------------------------------------------------
filter_seuratOB <- function(seurat_object)
{
	#filter data
	#if (tolower(sequencing) == 'cell')
	#{
		if (tolower(organism) == 'human')
		{
			seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
		}

		if (tolower(organism) == 'mouse')
		{
			seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
		}
	#}

	seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)
	seurat_object <- NormalizeData(seurat_object)
	seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(seurat_object)
	seurat_object <- ScaleData(seurat_object, features = all.genes)
	seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
	seurat_object <- RunUMAP(seurat_object, dims = 1:components)
	seurat_object <- FindNeighbors(seurat_object, dims = 1:components)
	seurat_object <- FindClusters(seurat_object, resolution=0.5)

	return(seurat_object)
}
#--------------------------------------------------------------------

# GET PK VALUE
#--------------------------------------------------------------------
get_pk <- function(seurat_object, sample_name, project)
{
	print('Calculating pk')
	## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
	sweep.res.list_bm <- paramSweep_v3(seurat_object, PCs = 1:components, sct = FALSE)
	sweep.stats_bm <- summarizeSweep(sweep.res.list_bm, GT = FALSE)
	bcmvn_bm <- find.pK(sweep.stats_bm)

	write.table(bcmvn_bm, paste(subdirectory2, '/', sample_name, '_', project, '_pk_values.txt', sep=''), append = FALSE, sep = '\t', dec = '.', quote=FALSE, row.names = FALSE, col.names = TRUE)

	pdf(file=paste(subdirectory1, '/pk_sweep_plot_', sample_name, '_', project, '.pdf', sep=''), onefile=TRUE, width=11, height=8.5)
	print(ggplot(bcmvn_bm, aes(x = pK, y = BCmetric)) + geom_point())
	dev.off()

	file_name = paste(subdirectory2, '/', sample_name, '_', project, '_pk_values.txt', sep='')
	pk <- as.numeric(system2('python', args = c('src/scripts/get_pk.py', file_name), stdout=TRUE) )
	write.table(pk, paste(subdirectory2, '/', sample_name, '_', project, '_pk_value.txt', sep=''), append = FALSE, sep = '\t', dec = '.', quote=FALSE, row.names = FALSE, col.names = FALSE)

	return(pk)
}
#--------------------------------------------------------------------

# RUN DATA THROUGH DOUBLETFINDER
#--------------------------------------------------------------------
run_doubletfinder <- function(seurat_object, doublet_rate, pk_value)
{
	print('running doubletFinder')
	## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
	pann_data = colnames(seurat_object@meta.data)[grepl("pANN", colnames(seurat_object@meta.data))]
	annotations <- seurat_object@meta.data$seurat_clusters
	homotypic.prop <- modelHomotypic(annotations)
	nExp_poi <- round(doublet_rate*nrow(seurat_object@meta.data))  
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	seurat_object <- doubletFinder_v3(seurat_object, PCs = 1:components, pN = pN, pK = pk_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	seurat_object <- doubletFinder_v3(seurat_object, PCs = 1:components, pN = pN, pK = pk_value, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

	return(seurat_object)
}
#--------------------------------------------------------------------

# Save object
#--------------------------------------------------------------------
save_object <- function(seurat_object, sample)
{
	print('saving output')
	file_name = paste(base_directory, '/', sample, '_DF.RDS', sep='')
	saveRDS(seurat_object, file = file_name)
}
#--------------------------------------------------------------------

# REMOVE DOUBLETS FROM SEURAT OBJECT
#--------------------------------------------------------------------
remove_doublets <- function(seurat_object, out_path, sample, project, soup_path)
{
	print('plotting  doublets')
	meta_data = colnames(seurat_object@meta.data)[grepl("DF.classification", colnames(seurat_object@meta.data))] #there could be more than 1
	f_name = paste(subdirectory1, '/', sample, '_', project, '_doublet.pdf', sep='')
	d_name <- paste(subdirectory2, '/', sample, '_', project, '_doublet_ids.txt', sep='')
	multi_doublet_ids <- ''

	#I need to append to file if there is more than one meta_data slot.
	if (file.exists(d_name))
  	{
		file.remove(d_name)
	}
	
	pdf(file=f_name, onefile=TRUE, width=11, height=8.5)
	print(DimPlot(seurat_object, group.by = meta_data) + NoAxes())
	print(VlnPlot(seurat_object, features = "nFeature_RNA", group.by = meta_data, pt.size = 0.1))
	dev.off()

	if (length(meta_data) == 1)
	{
		#need cell ids for doublets
		doublets = seurat_object[, seurat_object@meta.data[, meta_data] == "Doublet"]
		doublet_ids <- colnames(doublets)
		write.table(doublet_ids, file=d_name, quote=FALSE, sep='', row.names=FALSE)
	}

	if (length(meta_data) > 1)
	{
		for (meta in meta_data)
		{
			doublets = seurat_object[, seurat_object@meta.data[, meta] == "Doublet"]
			doublet_ids <- colnames(doublets)
			write.table(doublet_ids, file=d_name, quote=FALSE, sep='', append=TRUE, row.names=FALSE)
		}
	}

	# Create seurate object from soupX data
	print(list.files(path=soup_path))

	#so = ''
	#if (data_type == 'outs')
	#{
	sample.data <- Read10X(data.dir = soup_path)
	so <- CreateSeuratObject(counts = sample.data, project = project, min.cells = 3, min.features = 200)
	#}

#	if (data_type == 'no_clusters')
#	{
#		f=paste(in_path, 'filtered_feature_bc_matrix/', sep='')
#		filt <- Read10X(data.dir=f)
#		so <- CreateSeuratObject(counts=filt , project=project, min.cells=0, min.features=0)
#	}
#
#	if (data_type == 'h5')
#	{
#		f=paste(in_path, 'filtered_feature_bc_matrix.h5', sep='')
#		filt <-Read10X_h5(f, use.names=TRUE)
#		so <- CreateSeuratObject(counts=filt , project=project, min.cells=0, min.features=0)
#	}

	# Remove doublets from soupX Seurat object
	print('subsetting DF for singlets')

	ids <- read.csv(file=d_name, header=TRUE, sep='') #there may be multiple meta data
	doublet_ids <- unique(ids$x)
	seurat_object <- so[,!colnames(so) %in% doublet_ids]

	print('writing output in 10x format...')
	#note, if this already exists, it will fail sans overwrite=TRUE
	write10xCounts(x = seurat_object@assays$RNA@counts, path = out_path, overwrite = TRUE, version='3')
}
#--------------------------------------------------------------------

## FUNCTION CALLS
#--------------------------------------------------------------------
S <- tenX_to_seuratOB(soup_path, organism, components)
doublet_predictor <- determine_multiplet_rate(S)
S <- filter_seuratOB(S)
pk <- get_pk(S, sample, project) 
S_df <- run_doubletfinder(S, doublet_predictor, pk)
save_object(S_df, sample)
remove_doublets(S_df, doubletFinder_path, sample, project, soup_path)
#--------------------------------------------------------------------
