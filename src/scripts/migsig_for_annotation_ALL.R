args = commandArgs(trailingOnly=TRUE)

project <- ''
method <- ''
lib_path <- ''
storage <- ''
res_file <- ''
organism <- ''
subdirectory3 <- '' #figures_path 
integrated_object1 <- ''
integrated_object2 <- ''
integrated_object3 <- ''
#top100marker_files_files1 <- ''
#top100marker_files_files2 <- ''
#top100marker_files_files3 <- ''

#"Rscript {input.script} {params.project} {params.method} {params.rpath} {params.storage} {params.reso_file} {params.organism} {params.path} {params.sample_file1} {params.sample_file2} {params.sample_file3} 2> {log.log_output}     && touch {output.out_files}"

#"Rscript {input.script} {params.project} {params.method} {params.rpath} {params.storage} {params.reso_file} {params.organism} {params.path} {params.sample_file1} {params.sample_file2} {params.sample_file3} 2> {log.log_output}     && touch {output.out_files}"

#migsig_path = 'data/endpoints/' + PROJECT.lower() + '/analysis/' + PCA + '/tables/msig/'

# test if there is at least one argument: if not, return an error
#if (length(args) < 14) {
#  stop('At least 14 arguments must be supplied.', call.=FALSE)
#} else if (length(args)==14) {
if (length(args) < 10) {
  stop('At least 10 arguments must be supplied.', call.=FALSE)
} else if (length(args)==10) {
	project = args[1]
	method = args[2]
	lib_path = args[3]
	storage = args[4]
	res_file = args[5]
	organism = args[6]
	subdirectory3 = args[7]
	integrated_object1 = args[8]
	integrated_object2 = args[9]
	integrated_object3 = args[10]
	#top100marker_files_files1 = args[12] 
	#top100marker_files_files2 = args[13] 
	#top100marker_files_files3 = args[14] 
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
integration_methods <- c('RPCA', 'SCT', 'Standard')
integrated_objects <- c(integrated_object1, integrated_object2, integrated_object3)
integrated_objects <- sort(integrated_objects) #this should put everything in the correct order (RPCA, SCT, Standard)

top100_path <- gsub('msig', 'dge_plots', subdirectory3)
if (method == 'ALL')
{
	for (i in (1:length(integration_methods)))
	{
		S <- storage_retrieval(storage, integrated_objects[i])
		base_directory = paste0(subdirectory3, integration_methods[i], '/')

		for (r in res_list)
		{
			f_name = paste0(top100_path, '/', integration_methods[i], '/', project, '_top100_markers_', r, '_', integration_methods[i], '.txt')
			if (file.exists(f_name))
			{
				f_gsea(S, f_name, r, celltypes, 'C8', integration_methods[i])
				f_gsea(S, f_name, r, curated_gene_sets, 'C2_CP', integration_methods[i])
				f_gsea(S, f_name, r, hallmarks, 'hallmarks', integration_methods[i])
				f_gsea(S, f_name, r, ontology_gene_sets, 'C5', integration_methods[i])
			}

			if (file.exists(f_name) == FALSE)
			{
				print(paste0(f_name, ' does not exist! This migsig analysis cannot be completed'))
			}

		}
	}
}
#--------------------------------------------------------------------
