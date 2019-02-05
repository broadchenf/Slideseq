function BeadIndices=GetBeadsInRegion(DGE,Beads,GeneNames)
    ImageSize=[6030,6030];
    Malat1Image=PlotGeneFromName({'Malat1'},GeneNames,DGE,Beads,'PlotStyle','Default');

    disp('Draw mask')
    imagesc(Malat1Image)
    h = imfreehand; %draw something 
    mymask = h.createMask();        
    Locs=[Beads.Locations];
    BeadIndices=mymask(sub2ind(ImageSize,round(Locs(2,:)),round(Locs(1,:))))==1;