function [ClusterUniqueMappedDGE,ClusterUniqueMappedBeads,ClusterUniqueMappedIlluminaBarcodes,GeneNames, FactorWeight]=DGEByCutoff(PuckName,Cluster,Cutoff,varargin)

    try
        PuckDirectory=GetPuckDirectory(PuckName);
    catch
        MySplit=split(PuckName,'\');
        PuckName=MySplit{length(MySplit)};
        PuckDirectory=GetPuckDirectory(PuckName);
    end

    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadMappingFile"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadMappingFile=varargin{index+1};
    else
        BeadMappingFile=FindMostRecentMapping(PuckDirectory);
    end
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="CropSuffix"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        CropSuffix=varargin{index+1};
    else
        CropSuffix='';
    end
    
    [FactorLoadings,DGE,Beads,Barcodes,GeneNames]=GetFactorLoadings(PuckName,'CropSuffix',CropSuffix,'BeadMappingFile',BeadMappingFile);
    TotalFactors=sqrt(sum(FactorLoadings.^2,2));
    CellsMat=(FactorLoadings./TotalFactors)>Cutoff;
    GoodBeads=find(sum(CellsMat(:,Cluster),2)>0);
    ClusterUniqueMappedDGE=DGE(:,GoodBeads);
    ClusterUniqueMappedBeads=Beads(GoodBeads);
    ClusterUniqueMappedIlluminaBarcodes=Barcodes(GoodBeads);
    FactorWeight=FactorLoadings(:,Cluster);