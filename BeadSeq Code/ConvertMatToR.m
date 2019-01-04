function ConvertMatToR(PuckDirectory,varargin)
    %PuckDirectory should point to the root directory for the puck, e.g. C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180506_1
    %containing BijectiveMapping.mat
    %We take one variable argument, which is the name of the beadmapping
    %file. If left blank, we will use the most recent beadmapping file

    BeadMappingFile=FindMostRecentMapping(PuckDirectory);
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadMappingFile"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadMappingFile=varargin{index+1};
    end    
    
    CropSuffix='';
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="CropSuffix"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        CropSuffix=varargin{index+1};
    end    
    load(fullfile(PuckDirectory,BeadMappingFile,['BijectiveMapping',CropSuffix,'.mat']))
    
    coords=[UniqueMappedBeads.Locations]';
    barcode_locations = [UniqueMappedIlluminaBarcodes',mat2cell(coords(:,1),ones(size(coords,1),1)),mat2cell(coords(:,2),ones(size(coords,1),1))];
    barcode_locations = cell2table(barcode_locations, 'VariableNames',{'barcodes','xcoord','ycoord'});
    DGETable=array2table(UniqueMappedDGE,'VariableNames',UniqueMappedIlluminaBarcodes,'RowNames',cellstr(GeneNames'));
    writetable(barcode_locations, fullfile(PuckDirectory,BeadMappingFile,['BeadLocationsForR',CropSuffix,'.csv']));
    writetable(DGETable, fullfile(PuckDirectory,BeadMappingFile,['MappedDGEForR',CropSuffix,'.csv']),'WriteRowNames',true);
%Also write out UniqueMappedDGE