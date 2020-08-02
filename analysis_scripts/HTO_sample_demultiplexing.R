options(stringsAsFactors = FALSE)

library(Seurat)
library(ggplot2)
library(DropletUtils)

dirs <- c("/path/to/Pool1/outs/raw_feature_bc_matrix",
			"/path/to/Pool2/outs/raw_feature_bc_matrix",
			"/path/to/Pool3/outs/raw_feature_bc_matrix",
			"/path/to/Pool4/outs/raw_feature_bc_matrix" )

for ( i in 1:4){

	data.dir <- dirs[i]

	raw.mat <- Read10X(dat.dir)

	cell.calls <- DropletUtils::emptyDrops(raw.mat[[1]])

	Pool <- CreateSeuratObject(raw.mat[[1]], project = paste0("Pool", i), meta.data = as.data.frame(cell.calls))

	#Hashtag antibodies 1-5 in first matrix positions
	Pool[["HTO"]] <- CreateAssayObject(raw.mat[[2]][c(1:5), ])
	#CITE-seq protein expression panel in positions 6-onwards in matrix - split off into separate assay for later
	Pool[["ADT"]] <- CreateAssayObject(raw.mat[[2]][6:19, ])

	Pool <- subset(Pool, nCount_HTO > 0 & Pool@meta.data$FDR < 0.05 & Pool@meta.data$nCount_RNA > 500)

	Pool <- NormalizeData(object = Pool, assay = "HTO", normalization.method = "CLR")

	Pool <- HTODemux(Pool, assay = "HTO", positive.quantile = .99)

	png(paste0("Pool", i, "_heatmap.png")
	print(HTOHeatmap(Pool, assay = "HTO"))
	dev.off()

	Idents(Pool) <- "HTO_classification"

	Pool.sub <- subset(Pool, idents = "Negative", invert = TRUE)

	Pool.sub <- Pool.sub[, sample(size=2000, 1:length(Cells(Pool.sub)))]

	hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = Pool.sub, assay = "HTO"))))

	Pool.sub <- RunTSNE(Pool.sub, distance.matrix = hto.dist.mtx, perplexity = 100)

	png(paste0("Pool", i, "_tsne.png")
	print(DimPlot(Pool.sub, group.by = "hash.ID", reduction="tsne"))
	dev.off()

	png(paste0("Pool", i, "_all_hashtags_tsne.png")
	print(FeaturePlot(Pool.sub, rownames(GetAssayData(object = Pool.sub, assay = "HTO"))))
	dev.off()
	
	saveRDS(Pool, file=paste0("Pool", i, "_dehashed.RDS")

}

