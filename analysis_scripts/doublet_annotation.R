library(DoubletFinder)
library(ggplot2)
options(stringsAsFactors=FALSE)

#Quick param grid search for identifying within-sample doublets in hashed pools
Pool <- readRDS("Pool1.RDS")

dfs <- list()

expected <- c( 0.05, 0.1, 0.15, 0.20, .25, .3, .35, .4, .5, 0.8, .99)

pKs <- c( 0.0001, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2)

for ( i in expected){
  
	homotypic.prop <- modelHomotypic(Pool@active.ident)
	nExp_poi <- round(i*length(Cells(Pool)))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

	for (k in pKs){
	  
		Pool <- doubletFinder_v3(Pool, pN = 0.25, pK = k, nExp = nExp_poi, reuse.pANN = FALSE, PCs=1:10)
		Pool <- doubletFinder_v3(Pool, pN = 0.25, pK = k, nExp = nExp_poi.adj,PCs=1:10)

		Pool@meta.data[,"DF_hi.lo"] <- Pool@meta.data [, paste0("DF.classifications_0.25_", k, "_", nExp_poi.adj)]  
		Pool@meta.data$DF_hi.lo[which(Pool@meta.data$DF_hi.lo == "Doublet" & Pool@meta.data[, paste0("DF.classifications_0.25_", k, "_", nExp_poi)] == "Singlet")] <- "Doublet_lo"
		Pool@meta.data$DF_hi.lo[which(Pool@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"


		df <- as.data.frame(table( Pool$HTO_classification.global,Pool$DF_hi.lo) / colSums(table(Pool$DF_hi.lo, Pool$HTO_classification.global)))
		df$i <- i
		df$k <- k
		dfs[[as.character(paste0(i, "_", k))]] <- df
	}

}

dfs <- do.call(rbind, dfs)
ggplot(dfs[ dfs$Var1 != "Negative" & dfs$Var2 == "Doublet_hi", ], 
       aes(i, Freq, color=Var1, shape=factor(k))) + geom_point() + geom_line() + labs(x="Expected fraction of doublets",
                                                                     y="Percentage identified")