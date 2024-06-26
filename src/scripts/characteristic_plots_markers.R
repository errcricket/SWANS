# this script must be called from the base directory.
#(r message=FALSE)
# 80 
# 81 #define method
# 82 if (grepl('SCT', merged_file, fixed = TRUE) == TRUE)
# 83 {met = 'SCT'}
# 84 if (grepl('Standard', merged_file, fixed = TRUE) == TRUE)
# 85 {met = 'Standard'}
# 86 if (grepl('RPCA', merged_file, fixed = TRUE) == TRUE)
# 87 {met = 'RPCA'}`

args = commandArgs(trailingOnly=TRUE)
project <<- ''
method <<- ''
lib_path <<- ''
storage <<- ''
components <<- ''
organism <<- ''

# test if there is at least one argument: if not, return an error
if (length(args) < 6) {
  stop('At least six arguments must be supplied.', call.=FALSE)
} else if (length(args)==6) {
	project = args[1] 
	method = args[2] 
	lib_path = args[3]
	storage = args[4]
	components = args[5] 
	organism = args[6] 
}

library(Seurat, lib.loc=lib_path)
library(DropletUtils, lib.loc=lib_path)
library(patchwork)
library(magrittr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(data.table)
library(ggrepel)
library(metap)
library(multtest)
library(stringr)

# CREATE DIRECTORIES
#--------------------------------------------------------------------
base_directory=paste('data/endpoints/', project, '/analysis', sep='')
dir.create(base_directory, showWarnings=FALSE)
base_directory=paste(base_directory, '/PCA/', components, sep='')
dir.create(base_directory, showWarnings=FALSE)

subdirectory2 <<- ''

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

subdirectory1=paste(base_directory, '/figures', sep='')
subdirectory3=paste(base_directory, '/tables', sep='')
dir.create(subdirectory1, showWarnings=FALSE)
dir.create(subdirectory2, showWarnings=FALSE)
dir.create(subdirectory3, showWarnings=FALSE)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# IMPORT DATA

import_data <- function(method, filename)
{
	filename <<- ''
	S <<- ''

	if (method == 'Standard')
	{
		filename = paste(subdirectory2, '/', project, '_Standard', sep='')
	}
	if (method == 'RPCA')
	{
		filename = paste(subdirectory2, '/', project, '_RPCA', sep='')
	}
	if (method == 'SCT')
	{
		filename = paste(subdirectory2, '/L_Final_SCT', sep='')
	}


	if (tolower(storage) == 'h5')
	{
		filename <- paste(filename, '.H5', sep='')
		S = read_h5(file=filename, target.object='seurat')
	}

	if (tolower(storage) == 'rds')
	{
		filename <- paste(filename, '.RDS', sep='')
		S = readRDS(filename)
	}

	return(S)
}

S <- import_data(method, '')
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Generic function to plot object features

plot_characteristics <- function(seurat_object, directory, output_name) 
{
	print('Creating characteristic plots...')
	#Idents(object=seurat_object) <- ident #need to set this....add as argument

	pdf(file=paste(subdirectory1, '/dims_unnamed.pdf', sep=''), onefile=TRUE, width=11, height=8.5)
	plots <- DimPlot(seurat_object, group.by = 'Experiment') + ggtitle('Unnamed UMAP Clusters')
	print(plots & theme(legend.position = 'top') & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
	
	#file_name <- paste(subdirectory1, '/dimplot_by_experiment.pdf', sep='')
	#pdf(file=file_name, onefile=TRUE, width=11, height=8.5)
	print(DimPlot(seurat_object, reduction = 'umap', split.by = 'Experiment', label=TRUE, label.size=3, repel=TRUE) + NoLegend() + ggtitle('UMAP Split by Experiment')

	for (red in c('tsne', 'umap'))
	{   
		print(red)
		title1=paste('dimplot_unnamed_', red, sep='')
		#pdf(file=file_name1, onefile=TRUE, width=11, height=8.5)
		print(DimPlot(seurat_object, reduction=red, label=TRUE) + ggtitle(title1))  #+ plot_annotation(title = output_name))
	}   

	dev.off()
}

plot_characteristics(S, subdirectory1)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Get Cell Count Proportions

find_cell_count <- function(seurat_object, directory)
{
	c=table(Idents(seurat_object))
	y <- as.data.frame(t(c))
	names(y)[1] <- 'Proportion'
	names(y)[2] <- 'Cluster'
	names(y)[3] <- 'NumCell'

	p=round(prop.table(table(Idents(seurat_object))), 3)
	x <- as.data.frame(t(p))
	names(x)[2] <- 'Cluster'
	names(x)[3] <- 'Proportion'

	y['Proportion'] <- x$Proportion
	y=y[c(2,3,1)]
	
	file=paste(directory, '/', tolower(project), '_cluster_proportions.txt', sep='')
	write.table(y, file, append=FALSE, sep='\t', dec='.', row.names=FALSE, col.names=TRUE)
}

find_cell_count(S, subdirectory3)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Find All Markers by Cluster

find_all_markers <- function(seurat_object, directory, project)
{
	file=paste(directory, '/', tolower(project), '_markers.txt', sep='')
	markers <- FindAllMarkers(seurat_object, min.pct=0.01, logfc.threshold=0.25)
	markers %>% group_by(cluster) %>% top_n(n=75, wt=avg_logFC)
	write.table(markers, file, append=FALSE, sep='\t', dec='.', row.names=FALSE, col.names=TRUE)
}

find_all_markers(S, subdirectory3, project)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Find Conserved Markers by Cluster

find_conserved_markers <- function(seurat_object, directory)
{
	for (i in c(0:(length(levels(seurat_object@meta.data[['seurat_clusters']]))-1)))
	{
		print(i)
		file=paste(directory, '/', project, '_markers.txt', sep='')
		markers <- FindConservedMarkers(seurat_object, ident.1=i, grouping.var='Experiment')
		write.table(markers, file, append=FALSE, sep='\t', dec='.', row.names=TRUE, col.names=TRUE)
	}
}

find_conserved_markers(S, subdirectory3)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Create Heatmaps to show panglaodb markers
mark_cells <- function(seurat_object, markers, directory)
{  
	#Idents(object=S) <- 'SCT_snn_res.0.8' #need to set this by approach (SCT vs RPCA vs Standard)  and resolution

	f_name <- paste('data/genes/', markers, '.txt', sep='')
	markerS <- read.table(f_name, header=FALSE)
   gene_list=markerS$V1
	marker_list = gene_list

	if tolower(organism) == 'human' 
	{
		genes <- lapply(gene_list, toupper)
		marker_list = genes
	}

	if tolower(organism) == 'mouse'
	{
		genes <- lapply(gene_list, tolower())
		marker_list <- lapply(genes, str_to_title(genes)) #this may not work...
	}

   pdf(file=paste(directory, '/', markers, '_heatmap.pdf', sep=''), onefile=TRUE, width=11, height=8.5)
	try(print(DoHeatmap(subset(seurat_object, downsample=200), features=c(marker_list), size=3)))
   dev.off()

   #pdf(file=paste(directory, '/', markers, '_dotplot.pdf', sep=''), onefile=TRUE, width=11, height=8.5)
	#try(print(DotPlot(seurat_object, features=c(gene_list)) + RotatedAxis()))
   #dev.off()

   #pdf(file=paste(directory, '/', markers, '_featureplot.pdf', sep=''), onefile=TRUE, width=11, height=8.5)
	#try(print(FeaturePlot(seurat_object, label=TRUE, repel=TRUE, label.size=2.5, features=c(gene_list))))
   #dev.off()
}


mark_cells(S, 'Mammary_epithelial_cells', subdirectory1)
mark_cells(S, 'Pancreatic_stellate_cells', subdirectory1)
mark_cells(S, 'Undefined_placental_cells', subdirectory1)
mark_cells(S, 'Goblet_cells', subdirectory1)
mark_cells(S, 'T_helper_cells', subdirectory1)
mark_cells(S, 'Macrophages', subdirectory1)
mark_cells(S, 'Astrocytes', subdirectory1)
mark_cells(S, 'Anterior_pituitary_gland_cells', subdirectory1)
mark_cells(S, 'Tanycytes', subdirectory1)
mark_cells(S, 'Purkinje_fiber_cells', subdirectory1)
mark_cells(S, 'Erythroblasts', subdirectory1)
mark_cells(S, 'Pyramidal_cells', subdirectory1)
mark_cells(S, 'T_regulatory_cells', subdirectory1)
mark_cells(S, 'Decidual_cells', subdirectory1)
mark_cells(S, 'Oligodendrocytes', subdirectory1)
mark_cells(S, 'Thymocytes', subdirectory1)
mark_cells(S, 'Epithelial_cells', subdirectory1)
mark_cells(S, 'Erythroid-like_and_erythroid_precursor_cells', subdirectory1)
mark_cells(S, 'Spermatozoa', subdirectory1)
mark_cells(S, 'Ductal_cells', subdirectory1)
mark_cells(S, 'Stromal_cells', subdirectory1)
mark_cells(S, 'Luteal_cells', subdirectory1)
mark_cells(S, 'Peri-islet_Schwann_cells', subdirectory1)
mark_cells(S, 'Eosinophils', subdirectory1)
mark_cells(S, 'Neurons', subdirectory1)
mark_cells(S, 'Podocytes', subdirectory1)
mark_cells(S, 'Cholangiocytes', subdirectory1)
mark_cells(S, 'Cardiac_stem_and_precursor_cells', subdirectory1)
mark_cells(S, 'Trophoblast_cells', subdirectory1)
mark_cells(S, 'Granulosa_cells', subdirectory1)
mark_cells(S, 'Plasma_cells', subdirectory1)
mark_cells(S, 'Natural_killer_T_cells', subdirectory1)
mark_cells(S, 'Motor_neurons', subdirectory1)
mark_cells(S, 'Spermatocytes', subdirectory1)
mark_cells(S, 'Myofibroblasts', subdirectory1)
mark_cells(S, 'Meningeal_cells', subdirectory1)
mark_cells(S, 'NK_cells', subdirectory1)
mark_cells(S, 'Pericytes', subdirectory1)
mark_cells(S, 'Sebocytes', subdirectory1)
mark_cells(S, 'Hematopoietic_stem_cells', subdirectory1)
mark_cells(S, 'Basal_cells', subdirectory1)
mark_cells(S, 'Trichocytes', subdirectory1)
mark_cells(S, 'Fibroblasts', subdirectory1)
mark_cells(S, 'Mast_cells', subdirectory1)
mark_cells(S, 'Trophoblast_stem_cells', subdirectory1)
mark_cells(S, 'Delta_cells', subdirectory1)
mark_cells(S, 'Neuroblasts', subdirectory1)
mark_cells(S, 'Pluripotent_stem_cells', subdirectory1)
mark_cells(S, 'Epiblast_cells', subdirectory1)
mark_cells(S, 'Interneurons', subdirectory1)
mark_cells(S, 'Gamma_(PP)_cells', subdirectory1)
mark_cells(S, 'Acinar_cells', subdirectory1)
mark_cells(S, 'Endothelial_cells', subdirectory1)
mark_cells(S, 'Alpha_cells', subdirectory1)
mark_cells(S, 'Gamma_delta_T_cells', subdirectory1)
mark_cells(S, 'Hepatic_stellate_cells', subdirectory1)
mark_cells(S, 'Noradrenergic_neurons', subdirectory1)
mark_cells(S, 'His_bundle_cells', subdirectory1)
mark_cells(S, 'Purkinje_neurons', subdirectory1)
mark_cells(S, 'T_cells', subdirectory1)
mark_cells(S, 'Taste_receptor_cells', subdirectory1)
mark_cells(S, 'Beta_cells', subdirectory1)
mark_cells(S, 'Dopaminergic_neurons', subdirectory1)
mark_cells(S, 'Chondrocytes', subdirectory1)
mark_cells(S, 'Sertoli_cells', subdirectory1)
mark_cells(S, 'Pulmonary_alveolar_type_II_cells', subdirectory1)
mark_cells(S, 'Satellite_glial_cells', subdirectory1)
mark_cells(S, 'Juxtaglomerular_cells', subdirectory1)
mark_cells(S, 'Schwann_cells', subdirectory1)
mark_cells(S, 'Myoblasts', subdirectory1)
mark_cells(S, 'B_cells_memory', subdirectory1)
mark_cells(S, 'T_follicular_helper_cells', subdirectory1)
mark_cells(S, 'Tuft_cells', subdirectory1)
mark_cells(S, 'Distal_tubule_cells', subdirectory1)
mark_cells(S, 'Basophils', subdirectory1)
mark_cells(S, 'Trophoblast_progenitor_cells', subdirectory1)
mark_cells(S, 'Radial_glia_cells', subdirectory1)
mark_cells(S, 'Airway_goblet_cells', subdirectory1)
mark_cells(S, 'Ciliated_cells', subdirectory1)
mark_cells(S, 'Pancreatic_progenitor_cells', subdirectory1)
mark_cells(S, 'Glycinergic_neurons', subdirectory1)
mark_cells(S, 'Embryonic_stem_cells', subdirectory1)
mark_cells(S, 'B_cells', subdirectory1)
mark_cells(S, 'Clara_cells', subdirectory1)
mark_cells(S, 'Proximal_tubule_cells', subdirectory1)
mark_cells(S, 'T_memory_cells', subdirectory1)
mark_cells(S, 'Retinal_ganglion_cells', subdirectory1)
mark_cells(S, 'Parathyroid_chief_cells', subdirectory1)
mark_cells(S, 'Ependymal_cells', subdirectory1)
mark_cells(S, 'Oxyphil_cells', subdirectory1)
mark_cells(S, 'Hepatoblasts', subdirectory1)
mark_cells(S, 'Immature_neurons', subdirectory1)
mark_cells(S, 'Müller_cells', subdirectory1)
mark_cells(S, 'B_cells_naive', subdirectory1)
mark_cells(S, 'Trigeminal_neurons', subdirectory1)
mark_cells(S, 'Adipocytes', subdirectory1)
mark_cells(S, 'Microglia', subdirectory1)
mark_cells(S, 'Olfactory_epithelial_cells', subdirectory1)
mark_cells(S, 'Airway_epithelial_cells', subdirectory1)
mark_cells(S, 'Glutaminergic_neurons', subdirectory1)
mark_cells(S, 'Alveolar_macrophages', subdirectory1)
mark_cells(S, 'Glomus_cells', subdirectory1)
mark_cells(S, 'Gastric_chief_cells', subdirectory1)
mark_cells(S, 'Myocytes', subdirectory1)
mark_cells(S, 'Smooth_muscle_cells', subdirectory1)
mark_cells(S, 'Serotonergic_neurons', subdirectory1)
mark_cells(S, 'Foveolar_cells', subdirectory1)
mark_cells(S, 'Follicular_cells', subdirectory1)
mark_cells(S, 'Enterochromaffin_cells', subdirectory1)
mark_cells(S, 'Airway_smooth_muscle_cells', subdirectory1)
mark_cells(S, 'Pulmonary_alveolar_type_I_cells', subdirectory1)
mark_cells(S, 'Urothelial_cells', subdirectory1)
mark_cells(S, 'Dendritic_cells', subdirectory1)
mark_cells(S, 'Plasmacytoid_dendritic_cells', subdirectory1)
mark_cells(S, 'Satellite_cells', subdirectory1)
mark_cells(S, 'Adrenergic_neurons', subdirectory1)
mark_cells(S, 'Cholinergic_neurons', subdirectory1)
mark_cells(S, 'Enteric_glia_cells', subdirectory1)
mark_cells(S, 'Parietal_cells', subdirectory1)
mark_cells(S, 'Oligodendrocyte_progenitor_cells', subdirectory1)
mark_cells(S, 'Hepatocytes', subdirectory1)
mark_cells(S, 'Pulmonary_vascular_smooth_muscle_cells', subdirectory1)
mark_cells(S, 'Keratinocytes', subdirectory1)
mark_cells(S, 'Osteoclast_precursor_cells', subdirectory1)
mark_cells(S, 'Paneth_cells', subdirectory1)
mark_cells(S, 'T_cells_naive', subdirectory1)
mark_cells(S, 'T_cytotoxic_cells', subdirectory1)
mark_cells(S, 'Neural_stem_precursor_cells', subdirectory1)
mark_cells(S, 'Adipocyte_progenitor_cells', subdirectory1)
mark_cells(S, 'Platelets', subdirectory1)
mark_cells(S, 'Myoepithelial_cells', subdirectory1)
mark_cells(S, 'Endothelial_cells_(blood_brain_barrier)', subdirectory1)
mark_cells(S, 'Choroid_plexus_cells', subdirectory1)
mark_cells(S, 'Kupffer_cells', subdirectory1)
mark_cells(S, 'Crypt_cells', subdirectory1)
mark_cells(S, 'Reticulocytes', subdirectory1)
mark_cells(S, 'Neuroendocrine_cells', subdirectory1)
mark_cells(S, 'Transient_cells', subdirectory1)
mark_cells(S, 'Kidney_progenitor_cells', subdirectory1)
mark_cells(S, 'Enteroendocrine_cells', subdirectory1)
mark_cells(S, 'Retinal_progenitor_cells', subdirectory1)
mark_cells(S, 'Principal_cells', subdirectory1)
mark_cells(S, 'Osteoblasts', subdirectory1)
mark_cells(S, 'Epsilon_cells', subdirectory1)
mark_cells(S, 'Melanocytes', subdirectory1)
mark_cells(S, 'Chromaffin_cells', subdirectory1)
mark_cells(S, 'Photoreceptor_cells', subdirectory1)
mark_cells(S, 'Enterocytes', subdirectory1)
mark_cells(S, 'Merkel_cells', subdirectory1)
mark_cells(S, 'Mesothelial_cells', subdirectory1)
mark_cells(S, 'Endothelial_cells_(aorta)', subdirectory1)
mark_cells(S, 'GABAergic_neurons', subdirectory1)
mark_cells(S, 'Peritubular_myoid_cells', subdirectory1)
mark_cells(S, 'Nuocytes', subdirectory1)
mark_cells(S, 'Myeloid-derived_suppressor_cells', subdirectory1)
mark_cells(S, 'Luminal_epithelial_cells', subdirectory1)
mark_cells(S, 'Mesangial_cells', subdirectory1)
mark_cells(S, 'Megakaryocytes', subdirectory1)
mark_cells(S, 'Hemangioblasts', subdirectory1)
mark_cells(S, 'Germ_cells', subdirectory1)
mark_cells(S, 'Cajal-Retzius_cells', subdirectory1)
mark_cells(S, 'Ionocytes', subdirectory1)
mark_cells(S, 'Osteoclasts', subdirectory1)
mark_cells(S, 'Pinealocytes', subdirectory1)
mark_cells(S, 'Monocytes', subdirectory1)
mark_cells(S, 'Intercalated_cells', subdirectory1)
mark_cells(S, 'Salivary_mucous_cells', subdirectory1)
mark_cells(S, 'Neutrophils', subdirectory1)
mark_cells(S, 'Bergmann_glia', subdirectory1)
mark_cells(S, 'Loop_of_Henle_cells', subdirectory1)
mark_cells(S, 'Osteocytes', subdirectory1)
mark_cells(S, 'Vascular_smooth_muscle_cells', subdirectory1)
mark_cells(S, 'Microfold_cells', subdirectory1)
mark_cells(S, 'Cardiomyocytes', subdirectory1)
mark_cells(S, 'Leydig_cells', subdirectory1)
mark_cells(S, 'Langerhans_cells', subdirectory1)
mark_cells(S, 'Enteric_neurons', subdirectory1)
mark_cells(S, 'Red_pulp_macrophages', subdirectory1)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Calculate gene expression
calc_helper <- function(seurat_object, genes)
{
	counts=seurat_object[['RNA']]@counts
	ncells=ncol(counts)

	if(genes %in% row.names(counts)){
	sum(counts[genes,]>0)
	}else{return(NA)}
}

#genes <- c('Cd34', 'Kit')
#cd34 <- calc_helper(BM0H, 'Cd34')
#kit <- calc_helper(BM0H, 'Kit')
#--------------------------------------------------------------------

### PRINT QC METRICS (PRE/POST-FILTERING)
###--------------------------------------------------------------------
##dim_pre <- dim(BM_analysis)
#dim_pre0h <- dim(x=BM0H)
#dim_pre2h <- dim(x=BM2H)
#dim_pre4h <- dim(x=BM4H)
#dim_pre8h <- dim(x=BM8H)
#dim_pre24h <- dim(x=BM24H)
#dim_pre48h <- dim(x=BM48H)
#summary(BM0H)
#
#for (dim in c(dim_pre0h, dim_pre2h, dim_pre4h, dim_pre8h, dim_pre24h, dim_pre48h))
#{
#	print(dim)
#	write.table(dim, 'data/endpoints/tables/pre_dimensions.txt', append=TRUE, sep='\t', dec='.', row.names=TRUE, col.names=TRUE)
#}
##write.table(dim_pre, 'data/endpoints/tables/pre_dimensions.txt', append=FALSE, sep='\t', dec='.', row.names=TRUE, col.names=TRUE)

### View all samples @ once
### ------------------------------------------------------------------------------------------
##print('Creating comprehensive data plots')
##figure_directory=paste('data/endpoints/PCA', components, '/figures/', sep='')
##
##fname=paste(figure_directory, 'Features.pdf', sep='')
##pdf(file=fname, onefile=TRUE, width=11, height=8.5)
##print(VlnPlot(BM_analysis, features='nFeature_RNA', ncol=1))
##dev.off()
##
##fname=paste(figure_directory, 'Counts.pdf', sep='')
##pdf(file=fname, onefile=TRUE, width=11, height=8.5)
##print(VlnPlot(BM_analysis, features='nCount_RNA', ncol=1))
##dev.off()
##
##fname=paste(figure_directory, 'percent.mito.pdf', sep='')
##pdf(file=fname, onefile=TRUE, width=11, height=8.5)
##print(VlnPlot(BM_analysis, features='percent.mito', ncol=1))
##dev.off()
### ------------------------------------------------------------------------------------------
##
##f1=paste('data/endpoints/figures/individual_comprehensive_data_plots.pdf', sep='')
##pdf(file=f1, width=11, height=8.5) 
##
##counts_per_cell <- Matrix::colSums(BM_analysis)
##cat('counts per cell: ', counts_per_cell[1:7], '\n') ## counts for first 5 cells
##hist(log10(counts_per_cell+1), main='counts per cell', col='wheat')
##
##for (item in c('nFeature_RNA', 'nCount_RNA', 'percent.mito'))
##{
##	print(VlnPlot(BM_analysis, features=item, ncol=1))
##}
##
##print(FeatureScatter(BM_analysis, feature1='nCount_RNA', feature2='percent.mito'))
##print(FeatureScatter(BM_analysis, feature1='nCount_RNA', feature2='nFeature_RNA'))
##print(FeatureScatter(BM_analysis, feature1='nFeature_RNA', feature2='percent.mito'))
##
##for (ds in c(BM0H, BM2H, BM4H, BM8H, BM24H, BM48H))
##{
##
##	counts_per_cell <- Matrix::colSums(ds)
##	cat('counts per cell: ', counts_per_cell[1:7], '\n') ## counts for first 5 cells
##	hist(log10(counts_per_cell+1), main='counts per cell', col='wheat')
##
##	for (item in c('nFeature_RNA', 'nCount_RNA', 'percent.mito'))
##	{
##		print(VlnPlot(ds, features=item, ncol=1))
##	}
##
##	print(FeatureScatter(ds, feature1='nCount_RNA', feature2='percent.mito'))
##	print(FeatureScatter(ds, feature1='nCount_RNA', feature2='nFeature_RNA'))
##	print(FeatureScatter(ds, feature1='nFeature_RNA', feature2='percent.mito'))
##}
##dev.off()
#
### View all samples @ once
### ------------------------------------------------------------------------------------------
##figure_directory=paste('data/endpoints/PCA', components, '/figures/', sep='')
##
##fname=paste(figure_directory, 'Features_filtered.pdf', sep='')
##pdf(file=fname, onefile=TRUE, width=11, height=8.5)
##print(VlnPlot(BM_filtered, features='nFeature_RNA', ncol=1))
##dev.off()
##
##fname=paste(figure_directory, 'Counts_filtered.pdf', sep='')
##pdf(file=fname, onefile=TRUE, width=11, height=8.5)
##print(VlnPlot(BM_filtered, features='nCount_RNA', ncol=1))
##dev.off()
##
##fname=paste(figure_directory, 'percent.mito_filtered.pdf', sep='')
##pdf(file=fname, onefile=TRUE, width=11, height=8.5)
##print(VlnPlot(BM_filtered, features='percent.mito', ncol=1))
##dev.off()
### ------------------------------------------------------------------------------------------
##
##for (item in c('nFeature_RNA', 'nCount_RNA', 'percent.mito'))
##{
##	print(VlnPlot(BM_filtered, features=item, ncol=1))
##}
##
##f1=paste('data/endpoints/figures/individual_comprehensive_filtered_data_plots.pdf', sep='')
##pdf(file=f1, width=11, height=8.5) 
##
##counts_per_cell <- Matrix::colSums(BM_filtered)
##cat('counts per cell: ', counts_per_cell[1:7], '\n') ## counts for first 5 cells
##hist(log10(counts_per_cell+1), main='counts per cell', col='wheat')
##
##for (item in c('nFeature_RNA', 'nCount_RNA', 'percent.mito'))
##{
##	print(VlnPlot(BM_filtered, features=item, ncol=1))
##}
