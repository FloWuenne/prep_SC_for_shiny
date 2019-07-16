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
  'input4', 'i3', 2,'character',
  'output1', 'o1', 2, 'character',
  'output2', 'o2', 2, 'character',
  'output3', 'o3', 2, 'character',
  'output4', 'o4', 2, 'character'
), byrow=TRUE, ncol=4)

# Parse options
options = getopt(option_specification)

## Check whether user used Seurat or Scanpy
## If Seurat, update Object to v3
## If Scanpy, convert to seurat object, then update
filetype <- options$input2
if(filetype == "scanpy"){
  seurat_object <- ReadH5AD(file = options$input1)
  seurat_object <- UpdateSeuratObject(seurat_object)

} else if (filetype == "seurat"){
  seurat_object <- readRDS(options$input1)
  seurat_object <- UpdateSeuratObject(seurat_object)
}

## Create data frame containing embeddings, metadata and expression values for feather file interaction

## Check which type of embedding the user wants to use
embeddings <- options$input3
if(embeddings == "umap"){
  cell_embeddings <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
} else if (embeddings == "tsne"){
  cell_embeddings <- as.data.frame(seurat_object@reductions$tsne@cell.embeddings)
}


## Merge metadata, embeddings and normalised expression
metadata <- seurat_object@meta.data
metadata$cell_classification <- seurat_object@active.ident

## Select the normalization method the user chose
normalization <- options$input4
if(normalization == "RNA"){
  norm_data <- t(as.data.frame(as.matrix(seurat_object@assays$RNA@data)))
}else if(normalization == "SCT"){
  norm_data <- t(as.data.frame(as.matrix(seurat_object@assays$SCT@data)))
}

cell_embeddings_with_expression <- merge(cell_embeddings,metadata,by=0)
cell_embeddings_with_expression$cell_id <- cell_embeddings_with_expression$Row.names
rownames(cell_embeddings_with_expression) <- cell_embeddings_with_expression$Row.names
cell_embeddings_with_expression <- cell_embeddings_with_expression[2:ncol(cell_embeddings_with_expression)]
cell_embeddings_with_expression <- merge(cell_embeddings_with_expression,norm_data,by=0)

## get gene names for shiny
gene_names <- rownames(seurat_object@data)
gene_names_df <- data.frame("genes" = gene_names)

## Prepare sparse matrix for marker calculation
cell_embeddings_with_expression_genes <- cell_embeddings_with_expression[,gene_names_df$genes]
cell_embeddings_with_expression_genes_transposed_sparse <- as(t(cell_embeddings_with_expression_genes), "sparseMatrix")

## prepare a small dataset only containing cell IDs and initial clustering for shiny
shiny_user_clustering <- cell_embeddings_with_expression_genes %>%
  select(cell_id,cell_classification)

####
# Write output files

## 1) Feather file containing clustering and metadata
write_feather(cell_embeddings_with_expression,
              path = options$output1)

## 2) data table containing gene names
fwrite(gene_names_df,
       file = options$output2)

## 3) File containing the clustering for user defined cluster saving
write_feather(shiny_user_clustering,
              path = options$output3)

## 4) sparseMatrix for genes as .rds object to use for presto marker calculation
saveRDS(cell_embeddings_with_expression_genes_transposed_sparse,
        file = options$output4)

cat("\n Successfully transformed data! \n")

