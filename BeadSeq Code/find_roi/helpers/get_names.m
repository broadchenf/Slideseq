function [mapName,queryName] = get_names()

% Find the names of the last two modified image files in the directory.

D = dir;
[~,nameidx] = sort([D.datenum],'descend');
nameList = {D(nameidx).name};
imgName = {};
for i = 1:numel(nameList)
    name = lower(nameList{i});
    if numel(name) >= 4 && strcmp(name(end-3:end),'.tif') && ...
            ~strcmp(name,'result.tif')
        imgName{end+1} = name;
    end
end
if numel(imgName) < 2
    error('Images not found in folder.');
end

info1 = imfinfo(imgName{1});
info2 = imfinfo(imgName{2});
if info1.FileSize > info2.FileSize
    mapName = imgName{1};
    queryName = imgName{2};
else
    mapName = imgName{2};
    queryName = imgName{1};
end

end