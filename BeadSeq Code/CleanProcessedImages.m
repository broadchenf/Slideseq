%This script crawls through the "processed" image directory and deletes any
%files that aren't "stitched"

ProcessedFolder='D:\Slideseq\Processed';
PucksToSkip={'180618'};
ProcessedFolderDirectory=dir(ProcessedFolder);

for k=3:length(ProcessedFolderDirectory) %skip . and ..
    thisdir=[ProcessedFolder,'\',ProcessedFolderDirectory(k).name];
    if isdir(thisdir)&& ~contains(thisdir,PucksToSkip)
        filelist=dir([thisdir,'\*.tif']);
        for j=1:length(filelist)
            if ~contains(filelist(j).name,'transform')&&contains(filelist(j).name,'Stitched')
                delete([thisdir,'\',filelist(j).name])
            end
        end
    end
end
    