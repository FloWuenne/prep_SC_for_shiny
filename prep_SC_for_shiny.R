## How to run tool
# Rscript /cvmfs/soft.galaxy/v2.1/server/tools/ --input1 $seurat_object --output $shiny_data

# Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  
  library("getopt")
  library("Seurat")
  library("feather")
  library("data.table")
  library("dplyr")
})

# Import required libraries
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_specification = matrix(c(
  'input1', 'i1', 2, 'character',
  'input2', 'i2', 2, 'character',
  'output1', 'o1', 2, 'character',
  'output2', 'o2', 2, 'character',
  'output3', 'o3', 2, 'character',
  'output4', 'o4', 2, 'character'
), byrow=TRUE, ncol=4)

# Parse options
options = getopt(option_specification)

## Read in alevin sparse Matrix as RDS file
#seurat_file <- "/Users/florian_wuennemann/Postdoc/Genap/prep_SC_for_shiny/test_data/seurat_tsne_test.Rds"
# seurat_file <- "/Users/florian_wuennemann/Postdoc/Genap/data/large_seurat_test.Rds"
# seurat_object <- readRDS(seurat_file)
#marker_list <- fread("/Users/florian_wuennemann/Postdoc/Genap/data/test_marker_list.txt")

seurat_object <- readRDS(options$input1)
marker_list <- fread(options$input2)

## Extract mapping from Seurat object
cell_embeddings <- as.data.frame(seurat_object@dr$tsne@cell.embeddings)
metadata <- seurat_object@meta.data
metadata$cell_classification <- seurat_object@ident
norm_data <- t(as.data.frame(as.matrix(seurat_object@data)))

cell_embeddings_with_expression <- merge(cell_embeddings,metadata,by=0)
rownames(cell_embeddings_with_expression) <- cell_embeddings_with_expression$Row.names
cell_embeddings_with_expression <- cell_embeddings_with_expression[2:ncol(cell_embeddings_with_expression)]
cell_embeddings_with_expression <- merge(cell_embeddings_with_expression,norm_data,by=0)

## calculate cluster centers and other metadata
cluster_centers <- cell_embeddings_with_expression %>%
  dplyr::group_by(cell_classification) %>%
  summarise_at(vars(tSNE_1,tSNE_2),funs(mean(., na.rm=TRUE)))
cluster_centers <- as.data.frame(cluster_centers)

## get gene names
gene_names <- rownames(seurat_object@data)
gene_names_df <- data.frame("genes" = gene_names)

## Format Marker list
marker_list_formatted <- marker_list %>%
  select(-V1) %>%
  mutate("pct.diff" = pct.2 - pct.1)

# Write feather file
write_feather(cell_embeddings_with_expression,
path = options$output1)

write_feather(cluster_centers,
path = options$output2)

fwrite(gene_names_df,
       file = options$output3)

fwrite(marker_list_formatted,
       file = options$output4)

# ## Save embeddings with expression data
# write_feather(cell_embeddings_with_expression,
#               path = "/Users/florian_wuennemann/Postdoc/Genap/data/clustering_shiny.feather")
# 
# ## Save centroids and cluster information
# write_feather(cluster_centers,
#               path = "/Users/florian_wuennemann/Postdoc/Genap/data/cluster_info_shiny.feather")
# 
# fwrite(gene_names_df,
#        file = "/Users/florian_wuennemann/Postdoc/Genap/data/gene_names.tsv")
# 
# fwrite(marker_list_formatted,
#        file = "/Users/florian_wuennemann/Postdoc/Genap/data/marker_table.tsv")

cat("\n Successfully transformed data! \n")
