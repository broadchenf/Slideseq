function [FactorWeight,ClusterUniqueMappedDGE,ClusterUniqueMappedBeads,ClusterUniqueMappedIlluminaBarcodes,GeneNames]=GetFactorLoadings(PuckName,varargin)
    PuckDirectory=GetPuckDirectory(PuckName);

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
        
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="SaveOutput"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        SaveOutput=varargin{index+1};
    else
        SaveOutput=0;
    end
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Cluster"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Cluster=varargin{index+1};
    else
        Cluster=0; %Note, if cluster is specified, we subset the final DGE to the cluster(s) specified.
    end
    ClusterDataOutputFile=fullfile(PuckDirectory,BeadMappingFile,['FactorWeights',CropSuffix,'.mat']);
    load(fullfile(PuckDirectory,BeadMappingFile,['BijectiveMapping',CropSuffix,'.mat']))
    ClusterPath=fullfile(PuckDirectory,BeadMappingFile,'AnalogizerClusterAssignmentsOriginal.csv');
    opts=detectImportOptions(ClusterPath);
    opts=setvartype(opts,'max_factor','double');
    opts=setvartype(opts,'atlas_cluster','double');
    ClusterAssignments=readtable(ClusterPath,opts,'ReadVariableNames',true);
    ClusterAssignments.Properties.VariableNames={'var1', 'Barcode', 'max_factor', 'Cluster'};

   
    %NOTE THAT WE ARE FILTERING HERE: NMFReg will have thrown out some slideseq beads.
    [C,ia,ib]=intersect(ClusterAssignments.Barcode,UniqueMappedIlluminaBarcodes);
    ClusterAssignments=ClusterAssignments(ia,:);
    ClusterUniqueMappedDGE=UniqueMappedDGE(:,ib);
    ClusterUniqueMappedBeads=UniqueMappedBeads(ib);
    ClusterUniqueMappedIlluminaBarcodes=UniqueMappedIlluminaBarcodes(ib);
    %We now need to make sure the ClusterAssignments is in the same order as
    %the DGE
    if Cluster>0
        BeadsInCluster=ismember(ClusterAssignments.Cluster,Cluster);
        ClusterUniqueMappedDGE=ClusterUniqueMappedDGE(:,BeadsInCluster);
        ClusterUniqueMappedBeads=ClusterUniqueMappedBeads(BeadsInCluster);
        ClusterUniqueMappedIlluminaBarcodes=ClusterUniqueMappedIlluminaBarcodes(BeadsInCluster);
    end
         
    for j = 1:max(ClusterAssignments.max_factor)+1

       Assignment(j) = ClusterAssignments.Cluster( find(ClusterAssignments.max_factor == j-1,1));

    end
    
    
    
    Hsmatrix = readtable(fullfile(PuckDirectory,BeadMappingFile,'Analogizer_NMFreg_output\Hs30_0_0_17.csv'));
   
    [C,ia,ib] = intersect(Hsmatrix.Var1,ClusterUniqueMappedIlluminaBarcodes);
    Hsmatrix = Hsmatrix(ia,:);
    FactorWeight=zeros(size(Hsmatrix,1),max(Assignment));
    for k=1:max(Assignment)
        FactorWeight(:,k)  = sqrt(sum(table2array(Hsmatrix(:,find(Assignment == k)+1)).^2,2));
    end
        
    if SaveOutput
        save(ClusterDataOutputFile,'ClusterUniqueMappedDGE','ClusterUniqueMappedBeads','ClusterUniqueMappedIlluminaBarcodes','GeneNames')
    end