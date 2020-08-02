library(AUCell)
library(clusterProfiler)
library(ggplot2)

##seurat RDS object
cd8.seurat <- readRDS("cd8.seurat.RDS")
cells_rankings <- AUCell_buildRankings(cd8.seurat@data)

##load gene set, e.g. GSEA lists from BROAD
c5 <- read.gmt("\c5.all.v6.2.symbols.gmt") ## ALL GO

geneSets <- lapply(unique(c5$ont), function(x){print(x);c5$gene[c5$ont == x]})
names(geneSets) <- unique(c5$ont)

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

##set gene set of interest here for plotting
geneSet <- ""

aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
cd8.seurat$AUC <- aucs
png("Gene Set Activity.png", res=300, height=2000, width=4000)
ggplot(data.frame(cd8.seurat@meta.data, cd8.seurat@dr$umap@cell.embeddings), aes(UMAP1, UMAP2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 26) + facet_grid(.~Type)

dev.off()
