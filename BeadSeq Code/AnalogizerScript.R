
#Arguments are:
#1) Path to the slideseq data
# \\iodine-cifs\broad_macosko\\data\\clusters\\atlas_ica\\F_GRCm38.81.P60Cerebellum_ALT
#2) Path to the dropseq DGE, without escape characters
#   e.g. \\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT
#3) Path to the dropseq clusters
#   e.g. \\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT\assign\F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS
#4) Read cutoff 
#
#
#
#
#
#
#
#
#
#
#
#
#
#

args = commandArgs(trailingOnly=TRUE)

ReadCutoff=args[[4]]
DropseqDGEPath=args[[2]]
DropseqClusterPath=args[[3]]
downsamplecelltypes=0
fileseparator='\\'
SlideseqBasePath="C:\\Users\\sgr\\Dropbox (MIT)\\Project - SlideSeq\\Pucks\\Barcodes"
puckpath=file.path(SlideseqBasePath,args[[1]],fsep=fileseparator)

library(Analogizer)
library(Seurat)
library(Matrix)
library(FNN)

source("\\\\iodine-cifs\\broad_macosko\\evan\\SeuratExtrafunctions.R")


AnalogizerScript<-function(puckpath,AtlasDGEPath,AtlasClustersPath,ReadCutoff) {

cere = AtlasToSeuratWindows(DropseqDGEPath)
#setwd("/Users/rstickel/Documents/mac_lab/slideseq/data/seq_data/cell_alignment/180724_cerebellum") #change this to whatever directory you want the output to be
#this will read in the atlas cluster assingments for each bead
dropseq_clusts = readRDS(DropseqClusterPath)

if (downsamplecelltypes==1) {
  #this function will downsample the number of atlas beads for each cell type to a maximum
  #of the number the user defines in num_per_cluster
  stratified_sample = function(clusters,num_per_cluster,use_these_clusters=levels(clusters))
  {
    subset = c()
    for (i in use_these_clusters)
    {
      num_to_sample = min(sum(clusters==i),num_per_cluster)
      subset = c(subset,sample(names(clusters)[clusters==i],num_to_sample))
    }
    return(subset)
  }

  cere_subset = stratified_sample(dropseq_clusts[colnames(cere@raw.data)],2500)
  atlas.dge = as.matrix(cere@raw.data[,cere_subset])
  atlas.dge=atlas.dge[which(rowSums(atlas.dge)>0),]
  atlas.dge=atlas.dge[,which(colSums(atlas.dge)>ReadCutoff)]
} else {
#non subsampled atlas DGE
atlas.dge = as.matrix(cere@raw.data)
atlas.dge=atlas.dge[which(rowSums(atlas.dge)>0),]
atlas.dge=atlas.dge[,which(colSums(atlas.dge)>ReadCutoff)]
}



bijective<-read.table(file.path(puckpath,"BeadLocationsForR.csv",fsep=fileseparator), header= T,sep = ',', row.names = 1)
cere.slideseq.dge<-read.table(file.path(puckpath,"MappedDGEForR.csv",fsep=fileseparator),header = T,sep=',',row.names=1)
cere.slideseq.dge=cere.slideseq.dge[which(rowSums(cere.slideseq.dge)>0),]
cere.slideseq.dge=cere.slideseq.dge[,which(colSums(cere.slideseq.dge)>ReadCutoff)]
#NOTE: Bijective has to be filtered in the same way as cere.slideseq.dge

analogy = Analogizer(list(atlas = as.matrix(atlas.dge), ss= as.matrix(cere.slideseq.dge)))
analogy = normalize(analogy) #Normalize so each cell has unit expression
analogy = selectGenes(analogy,varthresh = 0.2)
analogy = scaleNotCenter(analogy) #Scale so variance of each gene is equal to 1. Because this is NMF, so everything has to be non-negatize.


#iNMF step via alternating least squares. Two parameters: k and lambda. k is the number of factors and lambda is the lagrange multiplier parameter in the ridge regression. Note that ALS is alternating least squares: you optimize V or K once, and then the other.
analogy = optimizeALS(analogy,k=30)  ##nrep=1 is the default and gives a good idea of the results. Recommend more than one initialization (nrep=10) for final analysesnrep=1 is the default and gi 
#iNMF puts out a list of H vectors, which have dimensions of cells by number of factors, but they are not completely comparable. If you run a tSNE on unaligned Hthey just separate by dataset.
#The goal here is to take cells that have similar factor expressions, and normalize cells with similar factor expressions so they will cluster together in the tSNE.
#To do that, given a cell X, you take 20 nearest neighbors y_i of X, and find the maximum factor weight for y_i. You then assign a vector to X which is the frequency of the maximum factors among the y_i. That vector is now used in graph clustering.
#The raw factor weights is @H.
analogy = quantile_align_SNF(analogy,ref_dataset='atlas') #SNF clustering and quantile alignment


analogy = run_tSNE(analogy)
#pdf("aligned_tsne_2.pdf")
#plotByDatasetAndCluster(analogy) #Can also pass in different set of cluster labels to plot
#dev.off()
#pdf("word_clouds.pdf")
#plot_word_clouds(analogy)
#dev.off()

#analogy = optimizeNewK(analogy,k=30) #Can also decrease K
#analogy = quantile_align_SNF(analogy) #SNF clustering and quantile alignment
#analogy = run_tSNE(analogy)
#plotByDatasetAndCluster(analogy) #Can also pass in different set of cluster labels to plot
#pdf("word_clouds_2.pdf")
#plot_word_clouds(analogy)
#dev.off()


dropseq_clusts = readRDS(DropseqClusterPath)
pdf(file.path(puckpath,"Analogizer_atlas_assignments_overlay_tsne.pdf"))
plotByDatasetAndCluster(analogy,clusters=dropseq_clusts)
dev.off()


analogy=AssignByAtlas(analogy,dropseq_clusts)
pdf(file.path(puckpath,"Analogizer_atlas_assignments_overlay_tsne_beadsAssigned.pdf"))
plotByDatasetAndCluster(analogy)
dev.off()


#Somehow the "..." in the function definition seems important?
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

#We want to find the significantly variable genes or significantly expressed genes for cluster k
NumClusts=length(levels(dropseq_clusts))
ExpressionByClusterDropseq=Matrix(0,dim(atlas.dge)[[1]],NumClusts)
rownames(ExpressionByClusterDropseq)<-dimnames(atlas.dge)[[1]]
VarianceByClusterDropseq=Matrix(0,dim(atlas.dge)[[1]],NumClusts)
rownames(VarianceByClusterDropseq)<-dimnames(atlas.dge)[[1]]


for (z in 1:length(levels(dropseq_clusts))) {
  ClustZ.DGE=atlas.dge[,which(colnames(atlas.dge) %in% names(dropseq_clusts[which(dropseq_clusts==z)]))]
  ExpressionByClusterDropseq[,z]=rowMeans(ClustZ.DGE)
  VarianceByClusterDropseq[,z]=RowVar(ClustZ.DGE)
}

write.csv(as.matrix(ExpressionByClusterDropseq),file=file.path(puckpath,"AnalogizerExpressionByDropseqCluster.csv",fsep=fileseparator))
write.csv(as.matrix(VarianceByClusterDropseq),file=file.path(puckpath,"AnalogizerVarianceByDropseqCluster.csv",fsep=fileseparator))
write.csv(analogy@clusters[which(names(analogy@clusters) %in% colnames(cere.slideseq.dge))],file=file.path(puckpath,"AnalogizerClusterAssignments.csv",fsep=fileseparator))

#END HERE
}

AnalogizerScript(puckpath,AtlasDGEPath,AtlasClustersPath,ReadCutoff)