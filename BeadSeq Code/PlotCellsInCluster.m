function BeadImage=PlotCellsInCluster(PuckName,Cluster,FigNum)

[GeneNames,UniqueMappedBeads,UniqueMappedDGE,UniqueMappedBarcodes]=LoadBijectiveMapping(PuckName);

ClusterPath=fullfile(GetPuckDirectory(PuckName),FindMostRecentMapping(GetPuckDirectory(PuckName)),'AnalogizerClusterAssignments.csv');
opts=detectImportOptions(ClusterPath);
opts=setvartype(opts,'x','double');
ClusterAssignments=readtable(ClusterPath,opts,'ReadVariableNames',true);
ClusterAssignments.Properties.VariableNames={'Barcode','Cluster'};

%NOTE THAT WE ARE FILTERING HERE: we have thrown out slideseq beads with
%<100 reads in the R program
[C,ia,ib]=intersect(ClusterAssignments.Barcode,UniqueMappedBarcodes);
assert(length(C)==length(ClusterAssignments.Barcode))
ClusterAssignments=ClusterAssignments(ia,:);
DGE=UniqueMappedDGE(:,ib);
Beads=UniqueMappedBeads(ib);
%We now need to make sure the ClusterAssignments is in the same order as
%the DGE
BeadsInCluster=ClusterAssignments.Cluster==Cluster;
DGE=DGE(:,BeadsInCluster);
Beads=Beads(BeadsInCluster);
figure(FigNum)
clf
DGEsums=sum(DGE,1);
for k=1:length(UniqueMappedBeads)
    if k/5==floor(k/5)
        rectangle('Position',[UniqueMappedBeads(k).Locations(1),UniqueMappedBeads(k).Locations(2),35,35],...
          'Curvature',[1,1], 'FaceColor','g','EdgeColor','None')
    end
end

for k=1:length(DGEsums)
    rectangle('Position',[Beads(k).Locations(1),Beads(k).Locations(2),25*min(500,DGEsums(k))/mean(DGEsums),25*min(DGEsums(k),500)/mean(DGEsums)],...
      'Curvature',[1,1], 'FaceColor','b')
end
axis([0,6000,0,6000]);

BeadImage=1;