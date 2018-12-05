path = "/Volumes/broad_macosko/data/clusters/atlas_ica/"
cerebellum_path = paste0(path, "F_GRCm38.81.P60Cerebellum_ALT/dge/F_GRCm38.81.P60Cerebellum_ALT.filtered.raw.dge.RDS")

cerebellum_dge = readRDS(cerebellum_path)
dim(cerebellum_dge)

clusterMarkers = readRDS("/Volumes/broad_macosko/data/clusters/atlas_ica/F_GRCm38.81.P60Cerebellum_ALT/F_GRCm38.81.P60Cerebellum_ALT_clusterMarkers.RDS")
length(clusterMarkers)
head(clusterMarkers)

selectedgenes = read.csv("/Volumes/broad_macosko/data/clusters/atlas_ica/F_GRCm38.81.P60Cerebellum_ALT/F_GRCm38.81.P60Cerebellum_ALT.selected_genes.txt", header = F)

dim(selectedgenes)
head(selectedgenes)

selectedgenesvec = as.vector(t(selectedgenes))

head(selectedgenesvec)
write.csv(selectedgenesvec, "/Volumes/broad_macosko/aleks/SpatialExpressionPatterns/data/Cerebellum_ALT_metacells/highlyvargenes.csv")

metacells = readRDS("/Volumes/broad_macosko/data/clusters/atlas_ica/F_GRCm38.81.P60Cerebellum_ALT/metacells/F_GRCm38.81.P60Cerebellum_ALT.metacells.RDS")
head(metacells)
dim(metacells)

metacells_hvgs = metacells[selectedgenesvec, ]
dim(metacells_hvgs)
head(metacells_hvgs)


cerebellum_dge_hvgs = cerebellum_dge[selectedgenesvec, ]
dim(cerebellum_dge_hvgs)
head(cerebellum_dge_hvgs)
cerebellum_dge_hvgs_t = t(cerebellum_dge_hvgs)
dim(cerebellum_dge_hvgs_t)
head(cerebellum_dge_hvgs_t)

cerebellum_dge_hvgs_t_m = as.matrix(cerebellum_dge_hvgs_t)
dim(cerebellum_dge_hvgs_t_m)
cerebellum_dge_hvgs_t_m[1:5, 1:5]

data_path_hvgs = "/Volumes/broad_macosko/aleks/NMFreg/data/cerebellum_dge_hvgs.csv"
#write.csv(selectedgenesvec, "/Volumes/broad_macosko/aleks/SpatialExpressionPatterns/data/Cerebellum_ALT_metacells/highlyvargenes.csv")
write.table(cerebellum_dge_hvgs_t_m, file = data_path_hvgs, sep = ",", row.names = TRUE, col.names = TRUE, quote=FALSE)


path = "/Volumes/broad_macosko/aleks/NMFreg/data/"
cluster_path = paste0(path, "F_GRCm38.81.P60Cerebellum_ALT.cell_cluster_outcomes.RDS")
cell_cluster_outcome = readRDS(cluster_path)
dim(cell_cluster_outcome)
cell_cluster_outcome_clean = cell_cluster_outcome[rownames(cerebellum_dge_hvgs_t_m),]
dim(cell_cluster_outcome_clean)

data_path_cluster = "/Volumes/broad_macosko/aleks/NMFreg/data/cell_cluster_outcome.csv"
write.table(cell_cluster_outcome_clean, file = data_path_cluster, sep = ",", row.names = TRUE, col.names = TRUE, quote=FALSE)
