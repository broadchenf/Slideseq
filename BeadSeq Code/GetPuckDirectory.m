function directory=GetPuckDirectory(PuckName)
try
    directory=fullfile('D:\Sam\Dropbox (Personal)\Projects\Project - SlideSeq\Pucks\Barcodes',PuckName);
    assert(exist(directory)>0);
catch
    directory=fullfile('C:\Users\sgr\Dropbox\Projects\Project - SlideSeq\Pucks\Barcodes',PuckName);
    assert(exist(directory)>0);
end