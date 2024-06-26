# SWANS
* *******************************************************************************
## NOTES: 
	* You will need to unzip cellmarker and pangelaoDB files.
* *******************************************************************************
## Table of Contents
* [Description](#description)
* [Requirements](#requirements)
* [Quick Start Guide](#quickstart)
    * Files  
		 * [configuration files](#configs) | [example_sample.sample_list](#sample_list)   
		 * [set-up](#setup)
    * [Modules/Targets and Output](#output)
		 * SoupX
		 * DoubletFinder
		 * Create Dynamic Script
		 * Merge Samples
		 * Get Significant PC
		 * Analyze Seurat Object
		 * DGE Round 1
* [Troubleshooting](#troubleshooting)
* [Tailoring Analysis Parameters](#tailor)
* [Reading Pause](#pause)
* [Post Annotation Analysis](#post_annotation)
* [HTML + Shiny Report](#HTML)
* Modules/Targets and Output II    
    * [Renaming Clusters and DGEA by Experimental Condition](#annotation)
* [Authors](#authors)
* [License](#license)

## <a name="description">Description: SWANS  Pipeline</a>
The pipeline/workflow or targets in this repo are available for your analyzing pleasure. These pipeline relies heavily on <a href="https://satijalab.org/seurat/index.html" target="_blank">Seurat 4.1.0</a>, and the analysis can be configured to work with single cell or single nucleus, human or mouse organisms, one or multiple samples, and a myriad of parameters that can be tailored to user needs. The pipeline must start with output data from <a href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger" target="_blank">CellRanger</a> (or have matrix, barcodes, and features files) and will finish with clustered data and DGE to help identify cell types. 

<!--- The targets can be operated individually (provided the necessary input is available) and are controlled overwhelmingly by Snakemake. See [Snakefile](#snakefile) for more information. -->

#### Analysis Overview
Contaminant RNA is removed from each sample -- making new matrix/features/barcodes files before identifying and removing doublets from the data (Note: 10X filtered and unfiltered outputs are required). Samples are then filtered to retain quality cells before being merged into one single cell object. A cursory analysis is done to calculate the percent of variation in each principle component, and the number of components carrying the majority of variation is identified and used to re-run the clustering analysis. Differential gene expression analysis is performed on each cluster to help identify cell types, and images used to visualize <a href="https://panglaodb.se/index.html" target="_blank">panglaoDB</a>,  <a href="http://bio-bigdata.hrbmu.edu.cn/CellMarker/" target="_blank">CellMarker</a>, and if available, user-supplied markers. 

## <a name="requirements">Requirements</a>
The pipelines require Python, Snakemake, bash, and R (plus several R packages like Seurat).

A list of all R tools and versions used in the workflow can be found in `date_Rsession_info.txt` however a brief overview of installed tools can be found below.

### R Software Tools & Installation
<details>
  <summary>Click to expand!</summary>

```
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux
```

attached base and other packages:

```
parallel  stats4    stats     graphics  grDevices utils     datasets methods   base     
glmGamPoi_1.2.0             reshape2_1.4.4    
biomaRt_2.46.3              RIOH5_0.1.3            
hdf5r_1.3.5                 data.table_1.14.2    
sctransform_0.3.3           patchwork_1.1.1           
DropletUtils_1.10.3         SingleCellExperiment_1.12.0
SummarizedExperiment_1.20.0 Biobase_2.50.0                 
GenomicRanges_1.42.0        GenomeInfoDb_1.26.2          
IRanges_2.24.1              S4Vectors_0.28.1            
BiocGenerics_0.36.0         MatrixGenerics_1.2.0     
matrixStats_0.58.0          sp_1.4-4                     
SeuratObject_4.1.0          Seurat_4.1.1     
stringr_1.4.0               DoubletFinder_2.0.3    
SoupX_1.5.2                 ggrepel_0.9.1               
ggplot2_3.3.3               dplyr_1.0.9           
msigdbr_7.5.1               fgsea_1.16.0
RColorBrewer_1.1-2              
```
#### Installation Tips

	* install.packages('withr', lib='path_to_installed_R_packages', dependencies=TRUE, repos='http://cran.rstudio.com/')   
<!--- 	* with_libpaths(new = 'path_to_installed_R_packages', install_github('JiekaiLab/RIOH5@HEAD', force=TRUE))    -->

	* Sys.unsetenv("GITHUB_PAT") #this is necessary to install anything with devtools.
	* devtools::install_github("hhoeflin/hdf5r")   
</details>
    
### Python Software Tools 
<details>
  <summary>Click to expand!</summary>

<a href="https://www.python.org/" target="_blank">Python (3.6.8)</a>,
<a href="https://snakemake.readthedocs.io/en/stable/#" target="_blank">Snakemake 7.9.0</a>

</details>

## <a name="quickstart">Quick Start Guide</a>

**This pipeline has been designed to work in a Linux Environment. This pipeline will not work in a Windows environment and has not been tested on the MacOS platform.** As a bare minumum, you must have some CellRanger outputs for each sample you would like to analyze, and have Python, snakemake, and R installed (and github to easily clone this repo). All R libraries have been made available. 


* Clone this repo into your working directory (instructions <a href="https://help.github.com/en/articles/cloning-a-repository" target="_blank">here</a>).

```
git clone ADD REPO ADDRESS HERE 
```

#### <a name="sample_list">samples.sample_list</a>
<details>
  <summary>Click to expand!</summary>
* Create a `samples.sample_list` file that will be placed at the top of your working directory (an example is included in this repo) and serves a few purposes: a) it defines which samples will be analyzed, b) it links each sample to the experimental condition of the samples, and c) it provides the location of the initial starting data. If there are more than one sample in the dataset, the `samples.sample_list` file should have 1 sample name per line along with the experimental condition and path to 10X data.   

**Tabs separate `samples`, `condition`, and `path_to_10X` with no extra spaces or empty lines**. 

Note that the pipeline expects an `outs` directory to exist in the folder containing the input 10X data; if none exists, you will need to create an outs directory for your input data and place the raw/filtered files therin. You do not need to include `/outs` in the 10X path in the `samples.sample_list` file.

```
samples  condition   path_to_10X
A1 control  full_path_to_A1_10X_data
A2 control  full_path_to_A2_10X_data
A3 control  full_path_to_A3_10X_data
B1 non_controls full_path_to_B1_10X_data
B2 non_controls full_path_to_B2_10X_data
```
</details>


#### <a name="configs">configs/example_local_configs.yaml</a>
<details>
  <summary>Click to expand!</summary>

* Customize the `configs/local_configs.yaml` file (an example file has been provided [configs/example_local_configs.yaml](configs/example_local_configs.yaml)) by supplying your email and other options that best correspond to the analysis you would like to perform (*e.g.,* single cell vs. single nucleus...).

```
# contact (email will be sent when jobs complete or fail) 
contact: user_name@foobar.foo

# R library location 
RPATH: /path_to_your_installed_R_packages/

# project name (IN LOWER CASE) (e.g., name of output directory under data/endpoints/project_name)
PROJECT: project_name

# starting files for soupX: (outs, no_clusters, h5)
MATRIX: outs

#organism (e.g., mouse, human)
ORGANISM: mouse

# mito cutoff (e.g., numeric value)
MITO: 15

# mito regression (y/n)
MITO_REGRESSION: y

# sequencing (e.g., cell, nucleus)
SEQUENCING: nucleus

# cell cycling regression (y/n)
CCREGRESSION: n

# cell regression method (standard/alternative)
CCREGRESSION_METHOD: standard

# feature thresholds 
# e.g., S <- subset(S.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
MIN_FEATURE_THRESHOLD: 200
MAX_FEATURE_THRESHOLD: 3000

# number of components for initial exploration
ICOMPONENTS: 50

# normalization method (e.g., RPCA, SCT, Standard, ALL)
METHOD: ALL

#resolution value(s) examples: single value: 0.5; if I, resolution = c(0.5, 0.8, 1.2), if blank, res = 0.8)
RESOLUTION: I

# find conserved genes across conditions (y/n)
CONSERVED_GENES: n

#cell marker gene set identity matching threshold
MARKER_THRESHOLD: 30

# full path to file(s) containing genes of interest
USER_GENE_FILE: genes_of_interest.txt
```
</details>

#### <a name="setup">Set-up</a>
<details>
  <summary>Click to expand!</summary>
* Call the `run_snakemake.sh` script for from the command line (`sh run_snakemake.sh`) and the following struture will be created for you. 

```
├── data 
     ├── endpoints
     │   └── project_name
     │       └── sample_name1 
     │           ├── 10X/outs 
     │       └── sample_name2 
     │           ├── 10X/outs
     │       └── sample_namen 
     │           ├── 10X/outs 
```

The `run_snakemake.sh` script will check to see that data has been added to the 10X input folders. If the path to the 10X data exists in the `samples.sample_list`, a symbolic link is created to the 10X data. If the folders are still empty (e.g., the paths were written incorrectly), you will receive the message below until there is 10X data for every sample listed in the `samples.sample_list` file. 
```
 'You have not placed any data in your 10X folders'
 'Add data to your 10X folders and run this again'
```
If you receive the message above, this is the moment to either fix the `samples.sample_list` file (preferrred) or manually add the 10X data to the appropriate folder (not ideal but doable). **This pipeline expects to see an `outs` directory under the input `10X` directory** and will automatically create an `outs` directory for the output data. If copying the 10X data into the output `outs` directory is cumbersome, make a symbolic link to the data (e.g., `ln -s`). Once 10X data are present, call the script again from the commandline (`sh run_snakemake.sh`).

NOTE The following files [configs/example_local_configs.yaml](configs/example_local_configs.yaml) and [example_samples.sample_list](example_samples.sample_list) are included in this repo. All files must be tailored for the end user, and once altered, **the names must be changed to `config/local_configs.yaml` and  `samples.sample_list` respectively**.
</details>

## <a name="output">Modules/Targets and Output </a>
<details>
  <summary>Click to expand!</summary>

#### <a name="soupX">SoupX</a>  
**Rule**: `src/rules/soupX.rules`  
**Scripts**: `src/scripts/soupX.R`, `src/rules/sample_list.py`   
**Input**: cellranger output in `path_to_10X/outs` folder   
**Output**: barcodes,genes,matrix files  under `data/endpoints/project/sample/soupX` 
    
This script will remove contaminant RNA using the <a href="https://github.com/constantAmateur/SoupX" target="_blank">SoupX</a> tool. By using the config file (`configs/local_configs.yaml`), you can run the 3 following options: `outs, no_clusters, h5`. If you have CellRanger output and you have the raw, filtered, clustering data, and other items in a folder called `outs`, use `outs` in the config file. If you have the raw and filtered CellRanger data, but none of the clustering information, choose `no_clusters`. If you only have h5 files -- put `h5` in the config **(be certain to place the raw/filtered h5 files in `path_to_10X_data/outs`**, you may need to create the `outs` directory).   

#### <a name="doubletfinder">DoubletFinder</a>
**Rule**: `src/rules/doubletFinder.rules`  
**Scripts**: `src/scripts/doubletFinder.R`  
**Input**: cellranger output in `path_to_10X/outs` folder   
**Output**: barcodes,features,matrix files under `data/endpoints/project/sample/doubletFinder` folder  
    
This script will remove cells identified as being a doublet using the <a href="https://github.com/chris-mcginnis-ucsf/DoubletFinder" target="_blank">DoubletFinder</a> tool. This tools needs to estimate the predicted doublet rate from the original CellRanger output so the 10X folder is used as input and the ids of cells identified as doublets are recorded and saved under `endpoints/project_name/sample/tables/sample_project_doublet_ids.txt`. The identified doublet cells are then removed from the soupX output and saved under `data/endpoints/project/sample/doubletFinder/`.

#### <a name="dynamic">Create Dynamic Script</a>
**Rule**: `src/rules/create_dynamic_script.rules`  
**Scripts**: `src/scripts/make_mergeR.py` or `src/scripts/make_mergeR_1.py`  
**Input**: `samples.sample_list`   
**Output**: `src/scripts/create_merged_dataset.R`  

This script will call one of two scripts (`make_mergeR_1.py`, `make_mergeR.py`) depending on the number of samples in `samples.sample_list`. The python scripts will dynamically create a new script called `src/scripts/create_merged_dataset.R` based on the names and the experimental conditions found in `samples.samples_list` (example below). Note that the `min.features` is defined in the `local_configs.yaml` file (MIN_FEATURE_THRESHOLD).

```
create_seurat_object <- function()
{
   S1.data = 'data/endpoints/project_name/A1/doubletFinder/'
   S1 <- Read10X(data.dir=S1.data)
   S1  <- CreateSeuratObject(counts=S1, project=project_name, min.cells=3, min.features=200)
   S1 <- AddMetaData(S1, metadata='control', col.name='Experiment')
   S1 <- AddMetaData(S1, metadata='A1', col.name='Samples')

   S2.data = 'data/endpoints/project_name/A2/doubletFinder/'
   S2 <- Read10X(data.dir=S2.data)
   S2  <- CreateSeuratObject(counts=S2, project=project_name, min.cells=3, min.features=200)
   S2 <- AddMetaData(S2, metadata='control', col.name='Experiment')
   S2 <- AddMetaData(S2, metadata='A2', col.name='Samples')

   S3.data = 'data/endpoints/project_name/A3/doubletFinder/'
   S3 <- Read10X(data.dir=S3.data)
   S3  <- CreateSeuratObject(counts=S3, project=project_name, min.cells=3, min.features=200)
   S3 <- AddMetaData(S3, metadata='control', col.name='Experiment')
   S3 <- AddMetaData(S3, metadata='A3', col.name='Samples')
   ...
   S.merged <- merge(S1, y = c(S2,S3,S4,S5), add.cell.ids = c('S1','S2','S3','S4','S5'), project = project_name)
   S.merged[['percent.mito']] <- PercentageFeatureSet(S.merged, pattern = '^mt-')
   ...
}
```

#### <a name="merge">Merge Samples</a>
**Rule**: `src/rules/merge_samples.rules`   
**Scripts**: `src/scripts/create_merged_dataset.R`   
**Input**:     
**Output**: `data/endpoints/project_name/analysis/figures/qc_1.pdf` `data/endpoints/project_name/analysis/figures/qc_2.pdf`, `data/endpoints/project_name/RDS/project_name_merged_samples.RDS` 

When the `src/rules/create_dynamic_script.rules` is called (`merge_samples.rules`), two files (qc_1.pdf, qc_2.pdf) containing information about nFeature_RNA, nCount_RNA, and percent.mito are created under `data/endpoints/project_name/analysis/figures/`. qc_1.pdf is the unfiltered data and qc_2.pdf is post filtering and the filtering values are user-defined in the `local_configs.yaml` file.
```
subset = nFeature_RNA > MIN_FEATURE_THRESHOLD & nFeature_RNA < MAX_FEATURE_THRESHOLD & percent.mito < MITO) 
```
This merged Seurat object is also saved as an RDS file in `data/endpoints/project_name/RDS/project_name_merged_samples.RDS`.

#### <a name="PC">Calculate Sig PC</a>
**Rule**: `src/rules/get_sig_PC.rules`   
**Scripts**: `src/scripts/calculate_sigPCs.R` or `src/scripts/calculate_sigPCs_1.R` depending on # of samples  
**Input**: `data/endpoints/project_name/analysis/RDS/project_name_merged_samples.RDS`.   
**Output**: `data/endpoints/project_name/analysis/sigPC.txt` and `data/endpoints/project_name/analysis/sigPCs.txt`  
    
This rule essentially does a quick analysis (normalizes, variable features, scaling, integrating) using a `Standard` approach with the sole purpose of determining how much variation is found in each principal component (sigPCs.txt) and returning the number of components containing the majority of variation (sigPC.txt). If you are interested in cell cycling, if `CCREGRESSION` is set to `y`, genes affiliated with cell cycling will be removed to ensure the PCs are not influenced by cell cycling genes. 

#### <a name="analyze_sc">Analyze Seurat Object</a>
**Rule**: `src/rules/analyze_sc_object.rules`  
**Scripts**: `src/scripts/analyze_SCdata.R`  
**Input**: `data/endpoints/project_name/analysis/sigPC.txt`, `data/endpoints/project_name/analysis/RDS/project_name_merged_samples./RDS`.   

**Output**: 
   * `data/endpoints/project_name/analysis/PCA_numeric_value/tables/project_name_METHOD_metatags.txt` (depending on the `METHOD` parameter in the `local_configs.yaml` file. )
   *  `data/endpoints/project_name/analysis/PCA_numeric_value/RDS/project_name_METHOD_RESOLUTION.RDS` (if ALL is selected, there will be a saved object for each of the three integration methods (`Standard`, `RPCA`, and `SCT`)).

This rule will create an integrated Seurat Object comprised of all the samples in the dataset. Depending on what integration `METHOD` is used, (Standard, RPCA, SCT, or ALL), the object will be normalized using the integration method found in the `local_configs.yaml` file. If 'ALL' is listed, 3 Seurat objects will be created, one for each integration method. If the data is unknown to the user, the flexibility of this pipeline will allow the user to compare different integration methods and determine which approach makes the most biological sense.   

If the user has `y` for `MITO_REGRESSION` or `CCREGRESSION`, mitochondria or cell cycling genes (resp.) will be regressed out. There are two approaches to looking at cell cycling (`standard`, `alternative`). Regardless of what is in the config file, this script will create images for the associated `METHOD`(s)  (e.g., `data/endpoints/project_name/analysis/PCA_X/figures/project_name_cell_cycling_METHOD_RESOLUTION_pre.pdf`, `data/endpoints/project_name/analysis/PCA_X/figures/project_name_cell_cycling_METHOD_RESOLUTION_post.pdf`) that shows how and whether the cells cluster by cell cycling pre and post regression. This provides the end user with an opportunity to understand their data and will be performed regardless of whether they requested cell cycling regression or not.
    
**NOTE**: The Snakefile should be updated to expect cell cycling images for this rule. The files are created but are not listed as a (input/output) requirement. 

#### <a name="dge">DGE Round 1: Characteristic Plots and Differential Gene Expression Analysis</a>
**Snakefile**: `src/rules/dge_plots.rules`  
**Scripts**: `src/scripts/create_images_DGE.R`   
**Input**: integrated Seurat Object(s) located `data/endpoints/project_name/analysis/PCA_numeric_value/RDS/project_name_METHOD_RESOLUTION.RDS` (if ALL is selected, there will be a saved object for each of the three integration methods (`Standard`, `RPCA`, and `SCT`)).  
**Output**: 

  * `data/endpoints/project_name/analysis/PCA_numeric/figures/dge_plots/METHOD/project_name_initial_cluster_plots_integrated_snn_res.RESOLUTION_METHOD.pdf` #UMAP/TSNE clustering images (1)
  * `data/endpoints/project_name/analysis/PCA_numeric/figures/dge_plots/METHOD/project_name_clusterProportions_integrated_snn_res.RESOLUTION_METHOD.pdf` # barplot image of cell counts per cluster (2)
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dge_plots/METHOD/project_name_clusterProportions_integrated_snn_res.RESOLUTION_METHOD.txt`  # numbers of cell counts per cluster (3)
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dge_plots/METHOD/project_name_markers_integrated_snn_res.RESOLUTION_METHOD.txt` #markers for each cluster (FindAllMarkers) 
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dge_plots/METHOD/project_name_top100_markers_integrated_snn_res.RES_METHOD.txt` #top 100 markers per cluster
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dge_plots/METHOD/conserved_markers/project_name_conservedMarkers_cluster_X_integrated_snn_res.RES_METHOD.txt`  #conserved markers by cluster (4)

This module will create cluster plots (UMAP, TSNE) and cluster proportions bar plots (.txt too) as well as finding DGE (`FindAllGenes`) and conserved genes (`FindConservedGenes`) for each cluster that will be used to help identify cell types. If there is more than one resolution (e.g., I (0.5, 0.8, 1.2)), there will be output for each resolution. 

**NOTES** 
	* If multiple integration methods and resolutions are selected, finding conserved markers (regardless of experimental condition) for each cluster is very time consuming. To reduce the analysis time, set `CONSERVED_GENES` to `n` in the `configs/example_local_configs.yaml` file. If the user chooses `y`, there is a conserved marker file for each cluster (item 4 above).   
	* Only upregulated genes are preserved during this step (trying to annotate cells, not full pathway analysis)
	* For items 1-3, there are corresponding files for each individual sample in the dataset where the sample name is at the end of the file name before the extension, (e.g., ...StandardA1.pdf). 
	* Also note that if `ALL` is selected as the `METHOD`, the above description applies to and will have `Standard`, `RPCA`, and `SCT` output.

#### <a name="markers">Plot Possible Markers</a>
**Snakefile**: `src/rules/calculate_means_plot_markers.rules`  
**Scripts**: `src/scripts/dataset_characterization.R`, `src/scripts/calculate_sample_mean.py`, `src/scripts/score_gene_sets.py`, `src/scripts/plot_cell_markers.R`   
**Input**: integrated Seurat Object(s) located `data/endpoints/project_name/analysis/PCA_numeric_value/RDS(H5)/project_name_METHOD_RESOLUTION.RDS` (if ALL is selected, there will be a saved object for each of the three integration methods (`Standard`, `RPCA`, and `SCT`)).  
**Output**: 

  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dataset_characterization/project_name_meanGE_RNA_clusters_integrated_snn_res.RESOLUTION_METHOD.txt` 
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dataset_characterization/project_name_meanGE_RNA_integrated_snn_res.RESOLUTION_METHOD.txt` 
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dataset_characterization/project_name_celltype_scores_cellmarker_integrated_snn_res.RESOLUTION_METHOD.txt`   
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dataset_characterization/project_name_celltype_scores_panglaoDBr_integrated_snn_res.RESOLUTION_METHOD.txt`   
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/dataset_characterization/project_name_possible_cell_types_integrated_snn_res.RESOLUTION_METHOD_MARKER_THRESHOLD.txt`  
  * `data/endpoints/project_name/analysis/PCA_numeric/figures/possible_markers/METHOD/project_name_pannotation_integrated_snn_res.RESOLUTION_celltype_name_METHOD.pdf`  

This pipeline uses gene sets from <a href="https://panglaodb.se/" target="_blank">PangeloDB</a> and <a href="http://bio-bigdata.hrbmu.edu.cn/CellMarker/" target="_blank">CellMarker</a> to identify potential cell types. Each cell type has its own file with affiliated genes therein. In the author's experiences, with rare exception, these databases have not yielded clear results.   

This rule...   
  * Runs scripts that calculate the mean gene expression for each cluster and calculates the mean for the entire dataset.   
  * Runs scripts that opens each celltype marker file (above) and fetches the gene names and **for each cluster** (and for every integration method and every resolution) and a) the expression of each gene is checked to see if the average value of the gene in the cluster is greater than the gene's average expression for the entire dataset; b) the gene set (genes in celltype marker file) is scored by counting how many genes meet criteria **a** divided by the total number of genes in the set; and c) if any cluster score for that cell type is greater than the `MARKER_THRESHOLD` parameter in the config file, the cell type is recorded and
  * The saved genes/cell types are plotted in heatmaps, dotplots, and featureplots to show the gene expression by cluster. Note that the rule is not expecting anything output except  `PROJECT.lower() + '_celltype_annotation_plot_dummy.txt'`; this is because it is not known in advance which files will be created. The created files are found here: `data/endpoints/project_name/analysis/PCA_numeric_value/figures/possible_markers/METHOD/`. 
  * If the user has provided a `USER_GENE_FILE`, the genes are plotted in heatmaps, dotplots, and featureplots to show the gene expression by cluster, the files are found here: `data/endpoints/project_name/analysis/PCA_numeric_value/figures/possible_markers/METHOD/user_defined_markers/`.
  * The top100 most differentially expressed genes (for each cluster) are run through msigdb (Hallmarks, Curated, Oncology, Celltypes) and if there are any significant hits, it is recorded and saved under `data/endpoints/project_name/analysis/PCA_numeric_value/tables/msig/`.
</details>

## 
### <a name="troubleshooting">Troubleshooting Pipeline</a>
<details>
  <summary>Click to expand!</summary>

The `Snakefile` contains a %$#@-ton of information and towards the bottom, there is something that looks like this:

```
include:
   "src/rules/soupX.rules"
include:
   "src/rules/doubletFinder.rules"
include:
   "src/rules/create_dynamic_script.rules"
include:
   "src/rules/merge_samples.rules"
include:
   "src/rules/get_sig_PC.rules"
include:
   "src/rules/analyze_sc_object.rules"
include:
   "src/rules/dge_plots.rules"
include:
   "src/rules/calculate_means_plot_markers.rules"

#--------------------MESSAGES-----------------------------------
onsuccess:
   print("The main controller pipeline completed with no errors.")
   shell("mail -s 'The main controller pipeline completed with no errors.' "+ config['contact']+" < {log}")

onerror:
   print("The main controller pipeline did not complete without errors."),
   shell("mail -s 'The main controller pipeline did not complete without errors, check the logs and try again.' "+ config['contact']+" < {log}")

#--------------------RULES---------------------------------------
rule biggie:
   input:
      final_files
#-------------------------------------------------------------------------------------
```

* ##### Swapping out Targets
Below is a list of targets (and the associated rule) that could be used instead of `final_files`. 
```
soupX_list, #soupX.rules step
doubletFinder_list, #doubletFinder.rules step
'src/scripts/create_merged_dataset.R', #create_dynamic_script.rules step
merge_list, #merge_samples.rules
sigs, #get_sig_PC.rules
sc_objects, #analyze_sc_object.rules
dge_files, #dge_plots.rules
characterization_files #calculate_means_plot_markers.rules
```

Each input item (e.g., target) listed above is defined earlier in the Snakefile and defines the output for each rule. The way Snakemake works is the output for one rule (e.g., `src/rules/create_dynamic_script.rules`) will be in input for the following rule (e.g., `src/rules/merge_samples.rules`). It is not necessary to list all the inputs/targets for the `biggie` rule  except `final_files` because the input for `src/rules/final_analysis.rules` is the output for `src/rules/calculate_means_plot_markers*.rules` and so on (kinda like a daisy chain); in short, calling `src/rules/final_analysis.rules` will call all the other rules until all the required inputs/targets have been created . However, if there is a need to troubleshoot a problem, the target can be changed. For example, if there is a problem with the DoubletFinder steps, below would only run the soupX and doubletFinder steps. Again, there is no need to explicitly add the soupx_list to the `biggie` rule as the files associated with this variable are required as input for the rule associated with doubletfinder_list.   
```
rule biggie:
   input:
      #soupX_list, #soupX.rules step
      doubletFinder_list, #doubletFinder.rules step
      #final_files
```

* ##### Redo sample merging
* If steps prior to `src/rules/merge_samples.rules` need to be re-run, **be certain to delete `src/scripts/create_merged_dataset.R`**.   

* ##### Checking Log Files
Each rule will create a log file specific to the rule. Regardless if the pipeline completes without an error, an email will be sent to the address in the config file, and will be one of the two below:
```
onsuccess:
   print("The main controller pipeline completed with no errors.")
   shell("mail -s 'The main controller pipeline completed with no errors.' "+ config['contact']+" < {log}")

onerror:
   print("The main controller pipeline did not complete without errors."),
   shell("mail -s 'The main controller pipeline did not complete without errors, check the logs and try again.' "+ config['contact']+" < {log}")
```

If there is a failure, there will be additional information in the email (as well as printing the error to the screen where the pipeline is being run). For example...  

```
rule biggie: (Hive)
include: "src/rules/calculate_means_plot_markers.rules": (Queen Bee)
rule msigdbr: (worker bee (located in src/rules/calculate_means_plot_markers.rules file)) 
```

There is one Hive and there are several Queen Bees, and each Queen Bee may have multiple worker bees/rules in the Queen Bee file*. Each of these sub-rules will create its own specific log file if there is an error. If there is a failure, the email will list something like this:

###### (\*These Queens play well together in one hive.)

```
Error in rule plot_cell_markers:
    jobid: 2
    output: project_name_celltype_annotation_plot_dummy.txt
    log: logs/celltype_annotation/project_name_plot_cell_markers.log (check log file(s) for error message)
    shell:
        Rscript src/scripts/plot_cell_markers.R project_name .....
```

Be sure to open up the log file and read about what the error is.   

* ##### Inspect Actual Command
The other item that will be in the email is the actual command that was used, you can see a bit of it above under `shell` (full command parameters not listed). If you are having issues running the pipeline to its conclusion, look at the actual command and see if any parameters are missing or considering just running the command from the commandline -- this may be helpful if you are trying to determine if the error is Snakemake or if it is a scripting issue.  
</details>

### <a name="tailor">Tailoring Analysis Parameters and Double Checking Analysis</a>
<details>
  <summary>Click to expand!</summary>
The most obvious starting place is in the configuration file (`configs/example_local_configs.yaml`). For example, if you are interested in changing the mitochondria filtering threshold, change the `MITO` value. In an earlier draft of this pipeline, if single nucleus data was placed in `SEQUENCING`, mitochondria thresholding was skipped. However, experience has shown that nuclear-encoded mitochondrial proteins will show up. IF exclusion of all mito genes are requested, a) the `MITO` value can be changed to 0 and/or choose `y` for `MITO_REGRESSION` and mitochondria-associated genes will be regressed out during analysis.  

   * If multiple integration methods and resolutions are selected, finding conserved markers (regardless of experimental condition) for each cluster is very time consuming. To reduce the analysis time, set `CONSERVED_GENES` to `n` in the `configs/example_local_configs.yaml` file.   

  * The DoubletFinder pk value for each sample is stored under `data/endpoints/sample_name/tables/sample_name_project_pk_value.txt`. This value should correspond to the highest peak on the `data/endpoints/sample_name/figures/pk_sweep_plot_sample_name_project_name.pdf`. If you find that this is not the case, change the value in the `sample_name_project_pk_value.txt` and re-run the analysis.

  * The number of significant components is located here: `data/endpoints/project_name/analysis/sigPC.txt`. If there is concern about loss of signal from excluding principal components for example, open this file and increase the value.  
   * Along these lines, Seurat authors suggested increasing the number of components for SCT integration. Instructions in ascertaining "how much should it be increased" is rather vague, so 10 is added to the number of significant components when using SCTransform. If you want to change this value, you can change it here: `src/scripts/analyze_SCdata.R`   

```
transform_object_sct <- function(S, directory, res, compos=compo)
{
	print('using SCT analysis approach...')
	SO.list <- SplitObject(S, split.by = 'Samples')

	library(glmGamPoi, lib.loc=lib_path)
	compo <- compo + 10 #number of sig PC slightly higher for SCT than Standard approach. Adding 10....
```
</details>

## 
### <a name="pause">Reading Pause</a>
It is not lost on the author(s) the amount of reading on this page, so take a quick break. <a href="https://www.youtube.com/watch?v=lm6IU6V-dE8&ab_channel=BurrnsLuciano" target="_blank">Reading Break</a>   
## 

### <a name="HTML">HTML + Shiny Report</a>
<details>
  <summary>Click to expand!</summary>
Interactive files are located under `html_report_project_name/` and will allow the user to compare multiple clustering schemas simultaneously in a R-studio environment. The interactive shiny app will provide users with an overview of YAML parameters, samples (names, conditions, recovered cells, DoubletFinder characteristics), QC metrics, PC selection, and an html file with hyperlinks to all output data. Additionally, if for example, the user employs all three Seurat's integration methods and chooses two resolutions, they will have 6 different clustering schemas. Using the previous example, the shiny app will give an overview in table format for the integration method, resolution, and number of clusters, number of cells per cluster for all 6 clustering arrangements. Lastly, the user can choose up to three schemas to compare simultaneously. The selected arrangement is rendered along with a table showing the differentially expressed genes by cluster.   

Individual .Rmd files for each sample will be created (under `html_report_project_name/`), however the user must use R-studio to convert (`knit`) them into html files.  

</details>


### <a name="post_annotation">Post Annotation Analysis</a>
<details>
  <summary>Click to expand!</summary>

At this point, the user has perused all the data that has been created in order to help identify cell types for each cluster. Moving forward, the user must select a final resolution, a final integration method (Standard, RPCA, SCT), (or provide a path to a different RDS file) and supply a text file that attaches an original cluster number to a cell type. This is defined in a second configuration file `configs/final_analysis_configs.yaml` (`configs/example_final_analysis_configs.yaml` has been included in this repo and requires the user to input final choices and must be renamed to `configs/final_analysis_configs.yaml`).   

```
# Final Normalization Method (e.g., RPCA, SCT, Standard)
FINAL_METHOD: Standard

#Final Resolution Value (e.g., 0.5)
FINAL_RESOLUTION: 0.5

#Path to Cluster Annotation File
CLUSTER_ANNOTATION_FILE: final_cluster_assignment.txt

#If user has performed additional work on seurat object
# Supply final Seurat Object (y/n)
USER_SUPPLIED_SEURAT_OBJECT: n

# full path & file name of user-supplied seurat object
PATH_USER_SUPPLIED_SEURAT_OBJECT:
```

The `CLUSTER_ANNOTATION_FILE` can have just about any name (the first character of the file name must be a letter (e.g., no numbers)). The header must be as below (cluster^Icelltype) where ^I represents a tab, and the cluster number and its celltype are also separated by a tab.

```
cluster  celltype
0  Fish
1  Fry 
2  Erin
3  Myke
4  Sound
5  Depeche
6  Cave
7  Hydra
8  Bicycle
9  Camping
10 Bianchi
11 Pear
```
</details>

## <a name="outputII">Modules/Targets and Output II</a>
#### <a name="annotation">Renaming Clusters and DGEA by Experimental Condition</a>
<details>
  <summary>Click to expand!</summary>

The information below is already inside the Snakefile.
```
include:
   "src/rules/final_anlysis.rules" 

rule biggie:
   input:
	data/endpoints/project_name/analysis/PCA_numeric/figures/final_analysis/project_name_final_cluster_plots_FINAL_RESOLUTION_FINAL_METHOD.pdf
	
```

**Snakefile**: `src/rules/final_analysis.rules`  
**Scripts**: `src/scripts/final_analysis.R`   
**Input**: characterization_files, integrated Seurate Object(s) located `data/endpoints/project_name/analysis/PCA_numeric_value/RDS(H5)/project_name_METHOD_RESOLUTION.RDS`,`FINAL_METHOD`, `FINAL_RESOLUTION`, `CLUSTER_ANNOTATION_FILE`   
**Output**: 
  * `data/endpoints/project_name/analysis/PCA_numeric/RDS/project_name_FINAL_METHOD_final.RDS`
  * `data/endpoints/project_name/analysis/PCA_numeric/figures/final_analysis/project_name_final_cluster_plots_FINAL_RESOLUTION_FINAL_METHOD.pdf` #UMAP/TSNE clustering images (1)
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/final_analysis/project_name_final_markers_CLUSTER_CELLTYPE.txt` #differentially expressed genes between conditions for each cluster (FindAllMarkers)
  * `data/endpoints/project_name/analysis/PCA_numeric/tables/final_analysis/project_name_significant_final_markers_CLUSTER_CELLTYPE.txt` #significantly differentially expressed genes between conditions for each cluster
</details>

## 
The author has done their best to provide out of the box code that _should_ work provided you have followed the directions. However, there are always items that were not considered. Try to solve the problem and if you cannot resolve it, feel free to contact me.

## <a name="authors">Authors-Contributors</a>
<a href="https://github.com/errcricket" target="_blank">Erin R. Reichenberger</a>
## <a name="license">License</a>
MIT License

Copyright (c) 2024

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# SWANS
