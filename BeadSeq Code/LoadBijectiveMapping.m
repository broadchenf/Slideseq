function [GeneNames,UniqueMappedBeads,UniqueMappedDGE,UniqueMappedIlluminaBarcodes]=LoadBijectiveMapping(PuckName,varargin)
    CropSuffix='';
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="CropSuffix"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        CropSuffix=varargin{index+1};
    end
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadMappingFile"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadMappingFile=varargin{index+1};
    else
        BeadMappingFile=FindMostRecentMapping(GetPuckDirectory(PuckName));
    end
    load(fullfile(GetPuckDirectory(PuckName),BeadMappingFile,['BijectiveMapping',CropSuffix,'.mat']))