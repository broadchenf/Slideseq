function directory=GetPuckDirectory(PuckName)
directory=fullfile('\\iodine-cifs\broad_macosko\data\Slideseq\Barcodes',PuckName);
assert(exist(directory)>0);