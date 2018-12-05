path = "/broad/macosko/data/clusters/atlas_ica/"
hippocampus_path = paste0(path, "F_GRCm38.81.P60Hippocampus/dge/F_GRCm38.81.P60Hippocampus.filtered.raw.dge.RDS")

hippocampus_dge = readRDS(hippocampus_path)
dim(hippocampus_dge)
clusterMarkers = readRDS("/broad/macosko/data/clusters/atlas_ica/F_GRCm38.81.P60Hippocampus/GRCm38.81.P60Hippocampus_clusterMarkers.RDS")
length(clusterMarkers)
head(clusterMarkers)

selectedgenes = read.csv("/broad/macosko/data/clusters/atlas_ica/F_GRCm38.81.P60Hippocampus/F_GRCm38.81.P60Hippocampus.selected_genes.txt", header = F)

dim(selectedgenes)
head(selectedgenes)

selectedgenesvec = as.vector(t(selectedgenes))

head(selectedgenesvec)

metacells = readRDS("/broad/macosko/data/clusters/atlas_ica/F_GRCm38.81.P60Hippocampus/metacells/F_GRCm38.81.P60Hippocampus.metacells.RDS")
head(metacells)
dim(metacells)

metacells_hvgs = metacells[selectedgenesvec, ]
dim(metacells_hvgs)
head(metacells_hvgs)

str(hippocampus_dge)
hippocampus_dge_hvgs = hippocampus_dge[selectedgenesvec,]
dim(hippocampus_dge_hvgs)
head(hippocampus_dge_hvgs)
hippocampus_dge_hvgs_t = t(hippocampus_dge_hvgs)
dim(hippocampus_dge_hvgs_t)
head(hippocampus_dge_hvgs_t)

hippocampus_dge_hvgs_t_m = as.matrix(hippocampus_dge_hvgs_t)
dim(hippocampus_dge_hvgs_t_m)
hippocampus_dge_hvgs_t_m[1:5, 1:5]

data_path_hvgs = "/broad/macosko/bstickels/data/slideseq/NMFreg/hippocampus/hippocampus_dge_hvgs.csv"

write.table(hippocampus_dge_hvgs_t_m, file = data_path_hvgs, sep = ",", row.names = TRUE, col.names = TRUE, quote=FALSE)


path = "/broad/macosko/data/clusters/atlas_ica/F_GRCm38.81.P60Hippocampus/assign/"
cluster_path = paste0(path, "F_GRCm38.81.P60Hippocampus.cell_cluster_outcomes.RDS")
cell_cluster_outcome = readRDS(cluster_path)
dim(cell_cluster_outcome)
cell_cluster_outcome_clean = cell_cluster_outcome[rownames(hippocampus_dge_hvgs_t_m),]
dim(cell_cluster_outcome_clean)

data_path_cluster = "/broad/macosko/bstickels/data/slideseq/NMFreg/hippocampus/cell_cluster_outcome.csv"
write.table(cell_cluster_outcome_clean, file = data_path_cluster, sep = ",", row.names = TRUE, col.names = TRUE, quote=FALSE)
