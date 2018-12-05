function tSNEandCluster(DGE,GenesToUse)

GeneCountSums=sum(DGE,2);
GeneCountSumsSorted=sort(GeneCountSums,'descend');

GenesToKeep=GeneCountSums>50;%Have to keep only the highly variable genes
BeadsToKeep=sum(DGE,1)>50;
DGEReduced=DGE(:,BeadsToKeep);

tSNEofDGE=tsne(DGEReduced(GenesToKeep,:)','NumDimensions',2,'NumPCAComponents',10);