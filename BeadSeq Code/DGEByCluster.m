function [ClusterUniqueMappedDGE,ClusterUniqueMappedBeads,ClusterUniqueMappedIlluminaBarcodes,GeneNames]=DGEByCluster(PuckDirectory,ClusterToAnalyze,varargin)

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

    if length(ClusterToAnalyze)==1
        ClusterDataOutputFile=fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_UniqueMappedBeads.mat']);
    end
    load(fullfile(PuckDirectory,BeadMappingFile,['BijectiveMapping',CropSuffix,'.mat']))
    ClusterPath=fullfile(PuckDirectory,BeadMappingFile,'AnalogizerClusterAssignments.csv');
    opts=detectImportOptions(ClusterPath);
    opts=setvartype(opts,'x','double');
    ClusterAssignments=readtable(ClusterPath,opts,'ReadVariableNames',true);
    ClusterAssignments.Properties.VariableNames={'Barcode','Cluster'};

    if ~exist('GeneNames') %This is basically only an issue for the 180413 pucks
        tableinfo = detectImportOptions('C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\Puck_180430_6\DGE.csv');
        tableinfo.SelectedVariableNames={'GENE'};
        T=readtable('C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\Puck_180430_6\DGE.csv',tableinfo);
        GeneNames=string(T.GENE);
    end

    %NOTE THAT WE ARE FILTERING HERE: we have thrown out slideseq beads with
    %<100 reads in the R program
    [C,ia,ib]=intersect(ClusterAssignments.Barcode,UniqueMappedIlluminaBarcodes);
    ClusterAssignments=ClusterAssignments(ia,:);
    UniqueMappedDGE=UniqueMappedDGE(:,ib);
    UniqueMappedBeads=UniqueMappedBeads(ib);
    UniqueMappedIlluminaBarcodes=UniqueMappedIlluminaBarcodes(ib);
    %We now need to make sure the ClusterAssignments is in the same order as
    %the DGE
    BeadsInCluster=ismember(ClusterAssignments.Cluster,ClusterToAnalyze);
    ClusterUniqueMappedDGE=UniqueMappedDGE(:,BeadsInCluster);
    ClusterUniqueMappedBeads=UniqueMappedBeads(BeadsInCluster);
    ClusterUniqueMappedIlluminaBarcodes=UniqueMappedIlluminaBarcodes(BeadsInCluster);
    if length(ClusterToAnalyze)==1
        save(ClusterDataOutputFile,'ClusterUniqueMappedDGE','ClusterUniqueMappedBeads','ClusterUniqueMappedIlluminaBarcodes','GeneNames')
    end