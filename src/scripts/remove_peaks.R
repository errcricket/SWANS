library(Seurat, lib.loc=lib_path)
library(DropletUtils, lib.loc=lib_path)

samples = c('A1', 'A2', 'A3', 'B1', 'B2')

for (s in samples)
{
	out_path_r <- paste(s, '/outs/raw_feature_bc_matrix/', sep='')
	out_path_f <- paste(s, '/outs/filtered_feature_bc_matrix/', sep='')

	print(out_path_r)
	print(out_path_f)

	# raw
	seurat_data <- Read10X(data.dir = paste0(s, '/raw_feature_bc_matrix'))
	counts <- seurat_data$`Gene Expression`
	seurat.object <- CreateSeuratObject(counts = counts, project = s)
	write10xCounts(x = seurat.object@assays$RNA@counts, path = out_path_r, version='3')
	
	# filtered 
	seurat_data <- Read10X(data.dir = paste0(s, '/filtered_feature_bc_matrix'))
	counts <- seurat_data$`Gene Expression`
	seurat.object <- CreateSeuratObject(counts = counts, project = s)
	write10xCounts(x = seurat.object@assays$RNA@counts, path = out_path_f, version='3')
}
