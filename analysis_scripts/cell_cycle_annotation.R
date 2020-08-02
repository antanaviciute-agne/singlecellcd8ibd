library(scran)
library(stringr)

#read in seurat object
cd8.seurat <- readRDS("cd8.seurat.RDS")
##read in cell cycle marker pairs
cc.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
##read in 10x cellranger features table to convert to ensembl identifiers
genes <- read.table("features.tsv.gz", header=FALSE)
##Subset to keep only mRNA tables and not feature barcoding/antibodies
genes <- genes[genes$V4 == "Expression", ]

##format gene names
gene.names <- genes$V1
names(gene.names) <- genes$V2
gene.names <- str_replace(gene.names, "\\.\\d+", "")

##run cell cycle predictions
cc <- cyclone( cd8.seurat@assays$RNA@data, pairs=cc.pairs, verbose=T,  gene.names=gene.names)
##store in seurat object
cd8.seurat$phases <- cc$phases
cd8.seurat$G1_score <- cc$normalized.scores$G1
cd8.seurat$$G2M_score <- cc$normalized.scores$G2M
cd8.seurat$S_score <- cc$normalized.scores$S
cd8.seurat$G1_score_raw <- cc$scores$G1
cd8.seurat$G2M_score_raw <- cc$scores$G2M
cd8.seurat$$S_score_raw <- cc$scores$S

saveRDS(cd8.seurat, file="cd8.seurat.RDS")