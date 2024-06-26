args = commandArgs(trailingOnly=TRUE)

project <- ''
method <- ''
lib_path <- ''
storage <- ''
res_file <- ''
organism <- ''
subdirectory3 <- '' #figures_path 
integrated_object <- ''
top100marker_files <- ''

# test if there is at least one argument: if not, return an error
if (length(args) < 9) 
{
  stop('At least 9 arguments must be supplied.', call.=FALSE)
} 

if (length(args)==9) 
{
	project = args[1]
	method = args[2]
	lib_path = args[3]
	storage = args[4]
	res_file = args[5]
	organism = args[6]
	subdirectory3 = args[7]
	integrated_object = args[8]
	top100marker_files = args[9] 
}

if (length(args)>9) 
{
	project = args[1]
	method = args[2]
	lib_path = args[3]
	storage = args[4]
	res_file = args[5]
	organism = args[6]
	subdirectory3 = args[7]
	integrated_object = args[8]
	top100marker_files = args[9:length(args)]
}

res_list <- read.table(res_file, header=FALSE)
res_list <- as.vector(res_list$V1)
print(res_list)

#--------------------------------------------------------------------
library(dplyr, lib.loc=lib_path)
library(Seurat, lib.loc=lib_path)
library(ggplot2, lib.loc=lib_path)
library(data.table, lib.loc=lib_path)
library(fgsea, lib.loc=lib_path)
library(msigdbr, lib.loc=lib_path)
#--------------------------------------------------------------------

# CREATE DIRECTORY
#--------------------------------------------------------------------
dir.create(subdirectory3, showWarnings=FALSE)

# CREATE MSIG LISTS
#--------------------------------------------------------------------
celltypes = msigdbr(species=organism, category = 'C8')
curated_gene_sets <- msigdbr(species=organism, category = 'C2', subcategory = 'CP')
hallmarks = msigdbr(species=organism, category = 'H')
ontology_gene_sets <- msigdbr(species=organism, category = 'C5')
#print(head(celltypes))
#print(head(curated_gene_sets))
#print(head(hallmarks))
#print(head(ontology_gene_sets))
#--------------------------------------------------------------------

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

# GSEA w/ MSIGDBR
#--------------------------------------------------------------------
f_gsea <- function(seurat_object, marker_file, res, gmt_df, gmt_name, met, directory=subdirectory3)
{
	directory=paste(directory, '/', met, sep='')
	dir.create(directory, showWarnings=FALSE)

	Idents(object=seurat_object) <- res 
	cluster_count <- unique(seurat_object@meta.data[[res]])
	print(cluster_count)
	dge_markers <-read.table(marker_file, header=T, sep='\t')

	for (c in cluster_count)
	{
		dge_cluster_markers <- subset(dge_markers, cluster == c)
		rankData <- dge_cluster_markers$avg_log2FC
		names(rankData) <- dge_cluster_markers$gene
		print(dge_cluster_markers$gene)
		  
		msigdbr_list <- split(x = gmt_df$gene_symbol, f = gmt_df$gs_name)
		  
		gsea_output <- fgsea(pathways = msigdbr_list, rankData, minSize=3)
		gsea_output <- subset(gsea_output, padj < 0.05)
		x1 <- length(gsea_output$padj)
		print('')
		print(x1)
		  
		if (x1 >= 1)
		{   
			print('making files')
			fname = paste(directory, '/', project, '_', res, '_', c, '_', gmt_name, '_', met, '.csv', sep='')
			print(fname)
			fwrite(gsea_output, file =fname)
		}   
	}
}
#--------------------------------------------------------------------

# FUNCTION CALLS
#--------------------------------------------------------------------
marker_files = top100marker_files

if (method != 'ALL')
{
	S <- storage_retrieval(storage, integrated_object)

	for (i in 1:length(res_list))
	{
		marker_file = marker_files[i]
		f_gsea(S, marker_file, res_list[i], celltypes, 'C8', method)
		f_gsea(S, marker_file, res_list[i], curated_gene_sets, 'C2_CP', method)
		f_gsea(S, marker_file, res_list[i], hallmarks, 'hallmarks', method)
		f_gsea(S, marker_file, res_list[i], ontology_gene_sets, 'C5', method)
	}
}
