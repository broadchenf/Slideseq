function [GeneNames,UniqueMappedBeads,UniqueMappedDGE,UniqueMappedIlluminaBarcodes]=LoadBijectiveMapping(PuckName,varargin)
    CropSuffix='';
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="CropSuffix"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        CropSuffix=varargin{index+1};
    end

    load(fullfile(GetPuckDirectory(PuckName),FindMostRecentMapping(GetPuckDirectory(PuckName)),['BijectiveMapping',CropSuffix,'.mat']))