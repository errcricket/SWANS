args = commandArgs(trailingOnly=TRUE)
project <- ''
lib_path <- ''
compo <<- ''
merged_file <- ''
storage <- ''
cell_cycle <<- ''
organism <- ''

subdirectory2 <<- ''
S <<- ''

# test if there is at least one argument: if not, return an error
if (length(args) < 7) {
  stop('At least seven arguments must be supplied.', call.=FALSE)
} else if (length(args)==7) {
   project = args[1] # starting data_type (out, filtered, merged_file)
   lib_path = args[2]
   compo = args[3]
   merged_file = args[4] 
   storage = args[5] 
	cell_cycle = args[6]
	organism = args[7]
}

compo <- as.numeric(compo)

library(dplyr, lib.loc=lib_path)
library(Seurat, lib.loc=lib_path)
library(patchwork, lib.loc=lib_path)
library(data.table, lib.loc=lib_path)

# CREATE DIRECTORIES
#--------------------------------------------------------------------
# Output that will be unchanged by component #s go here
base_directory=paste('data/endpoints/', project, '/analysis/', sep='')
dir.create(base_directory, showWarnings=FALSE)

subdirectory1=paste(base_directory, '/figures', sep='')

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
subdirectory3=paste(base_directory, '/tables', sep='')
dir.create(subdirectory1, showWarnings=FALSE)
dir.create(subdirectory2, showWarnings=FALSE)
dir.create(subdirectory3, showWarnings=FALSE)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# IMPORT DATA

import_data <- function(seurat_object, store, filename)
{
	print('READING DATA')

	if (tolower(store) == 'h5')
	{
		print('H5')
		seurat_object = read_h5(file=filename, target.object='seurat')
		#seurat_object = read_h5(file=filename, target.object='seurat', assay_name='RNA')
	}

	if (tolower(store) == 'rds')
	{
		print('RDS')
		seurat_object = readRDS(filename)
	}

	return(seurat_object)
}

S <- import_data(S, storage, merged_file)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Save object

save_object <- function(seurat_object, storage)
{	
	if (tolower(storage) == 'h5')
	{
		write_h5(data=seurat_object, object.type = 'seurat', file = paste(subdirectory2, '/', tolower(project), '_transformed.h5', sep=''))
	}

	if (tolower(storage) == 'rds')
	{
		saveRDS(seurat_object, file = paste(subdirectory2, '/', tolower(project), '_transformed.RDS', sep=''))
	}
}

#--------------------------------------------------------------------

# remove cell cycing genes
remove_cc_genes <- function(seurat_object)
{
	if (cell_cycle == 'y')
	{
		print('removing cell cycling genes')
		convertHumanGeneList <- function(x)
		{
			suppressMessages(library('biomaRt', lib.loc=lib_path))
			human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl', host = 'https://dec2021.archive.ensembl.org/')
			mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl', host = 'https://dec2021.archive.ensembl.org/')
			genesV2 = getLDS(attributes = c('hgnc_symbol'), filters = 'hgnc_symbol', values = x , mart = human, attributesL = c('mgi_symbol'), martL = mouse, uniqueRows=T)
			humanx <- unique(genesV2[, 2])

			return(humanx)
		}

		s.genes <- cc.genes$s.genes
		g2m.genes <- cc.genes$g2m.genes

		if (tolower(organism) == 'mouse')
		{  
			s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
			g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
		}

		toRemove <- c(s.genes, g2m.genes)
		counts <- GetAssayData(seurat_object, assay = 'RNA')
		counts <- counts[-(which(rownames(counts) %in% toRemove)),] #c('HBB','HBA1','HBA2','HBG2','HBG1','HBD'))),]
		seurat_object <- subset(seurat_object, features = rownames(counts))
	}

	return(seurat_object)
}

S <- remove_cc_genes(S)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Generic function to run PCA-UMAP

transform_object <- function(seurat_object)
{
	#print(length(seurat_object@meta.data[["nCount_RNA"]]))
	print('splitting samples by experiment...')

	SO.list <- SplitObject(seurat_object, split.by = 'Experiment')

	print('normalizing')
	SO.list <- lapply(X = SO.list, FUN = function(x) {
			 x <- NormalizeData(x)
			 x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
	})

	print('features')
	features <- SelectIntegrationFeatures(object.list = SO.list)
	print(features)

	print('finding anchors...')
	SO.anchors <- FindIntegrationAnchors(object.list = SO.list, anchor.features = features)
	print(SO.anchors)

	print('integrating...')
	SO.integrated <- IntegrateData(anchorset = SO.anchors)
	print(SO.integrated)

	DefaultAssay(SO.integrated) <- 'integrated'

	SO.integrated <- ScaleData(SO.integrated, verbose = TRUE)

	print('running PCA')
	SO.integrated <- RunPCA(SO.integrated, npcs = compo, verbose = TRUE)

	elbow_file <- paste(subdirectory1, '/', project, '_elbow_plot.pdf', sep='')
	pdf(file=elbow_file, onefile=TRUE, width=11, height=8.5)
	print(ElbowPlot(SO.integrated, ndims = compo, reduction = 'pca'))
	dev.off()

	save_object(SO.integrated, storage)

	return(SO.integrated)
}

S <- transform_object(S)
#S <- readRDS('data/endpoints/chaefer/analysis/RDS/chaefer_transformed.RDS')
#--------------------------------------------------------------------

#--------------------------------------------------------------------
get_sigPCs <- function(seurat_object)
{
	print('sigPCs')
	print(seurat_object)

	# Determine percent of variation associated with each PC
	pct <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100

	# Calculate cumulative percents for each PC
	cumu <- cumsum(pct)
	df <- cbind(pct, cumu)
	print('writing PC values')
	write.csv(df, file=paste(subdirectory3, '/sigPCs.txt', sep=''), quote=FALSE, row.names=FALSE)

	# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
	co1 <- which(cumu > 90 & pct < 5)[1]

	co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

	pcs <- min(co1, co2)

	print('writing out final PC value')
	write.csv(pcs, file=paste(subdirectory3, '/sigPC.txt', sep=''), quote=FALSE, row.names=FALSE)

	print(pcs)
}

get_sigPCs(S)
