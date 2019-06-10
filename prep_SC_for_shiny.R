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
  'input3', 'i3', 2,'character',
  'input4', 'i4', 2,'character',
  'output1', 'o1', 2, 'character',
  'output2', 'o2', 2, 'character',
  'output3', 'o3', 2, 'character'
), byrow=TRUE, ncol=4)

# Parse options
options = getopt(option_specification)

## Read in alevin sparse Matrix as RDS file
#seurat_file <- "/Users/florian_wuennemann/Postdoc/Genap/prep_SC_for_shiny/test_data/seurat_tsne_test.Rds"
# seurat_file <- "/Users/florian_wuennemann/Postdoc/Genap/data/large_seurat_test.Rds"
# seurat_object <- readRDS(seurat_file)
#marker_list <- fread("/Users/florian_wuennemann/Postdoc/Genap/data/test_marker_list.txt")

## Check whether user used Seurat or Scanpy
marker_list <- fread(options$input2)

filetype <- options$input3
if(filetype == "scanpy"){
  pbmc3k <- ReadH5AD(file = options$input1)
  marker_list_formatted <- marker_list
} else if (filetype == "seurat"){
  seurat_object <- readRDS(options$input1)
  seurat_object <- UpdateSeuratObject(seurat_object)
  ## Format Marker list
  marker_list_formatted <- marker_list %>%
    select(-V1)
}

## Extract mapping from Seurat object

## If Scanpy get UMAP, if Seurat get tSNE
## Check whether user used Seurat or Scanpy
embeddings <- options$input4
if(embeddings == "umap"){
  cell_embeddings <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
} else if (embeddings == "tsne"){
  cell_embeddings <- as.data.frame(seurat_object@reductions$tsne@cell.embeddings)
}

metadata <- seurat_object@meta.data
metadata$cell_classification <- seurat_object@active.ident
norm_data <- t(as.data.frame(as.matrix(seurat_object@assays$RNA@data)))

cell_embeddings_with_expression <- merge(cell_embeddings,metadata,by=0)
rownames(cell_embeddings_with_expression) <- cell_embeddings_with_expression$Row.names
cell_embeddings_with_expression <- cell_embeddings_with_expression[2:ncol(cell_embeddings_with_expression)]
cell_embeddings_with_expression <- merge(cell_embeddings_with_expression,norm_data,by=0)

## calculate cluster centers and other metadata

## get gene names
gene_names <- rownames(seurat_object@data)
gene_names_df <- data.frame("genes" = gene_names)

# Write feather file
write_feather(cell_embeddings_with_expression,
path = options$output1)

fwrite(gene_names_df,
       file = options$output2)

fwrite(marker_list_formatted,
       file = options$output3)

cat("\n Successfully transformed data! \n")
