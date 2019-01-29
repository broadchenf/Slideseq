library(Seurat)
library(Matrix)
library(FNN)

source("\\\\iodine-cifs\\broad_macosko\\evan\\SeuratExtrafunctions.R")

ClusterToGrab=5;

path="\\\\iodine-cifs\\broad_macosko\\data\\clusters\\atlas_ica\\F_GRCm38.81.P60Hippocampus"


cere=AtlasToSeuratWindows(path)

clusterpath='\\\\iodine-cifs\\broad_macosko\\data\\clusters\\atlas_ica\\F_GRCm38.81.P60Hippocampus\\assign\\F_GRCm38.81.P60Hippocampus.cluster.assign.RDS'
dropseq_clusts = readRDS(clusterpath)

#Subset the DGE:
atlas.dge = as.matrix(cere@raw.data)

Clust.DGE=atlas.dge[,which(colnames(atlas.dge) %in% names(dropseq_clusts[which(dropseq_clusts==ClusterToGrab)]))]

atlas.dge=atlas.dge[which(rowSums(atlas.dge)>0),]
atlas.dge=atlas.dge[,which(colSums(atlas.dge)>ReadCutoff)]
