function ImageOut=FindTranslationXCorr(ImageIn,varargin)
%We are given a stack of images
%We align them using xcorr.
%We then output the aligned stack.
%Input and output should both be uint16.
ImageOut=zeros(size(ImageIn));
ImageOut(:,:,1)=ImageIn(:,:,1);

b=[1,size(ImageIn,1),1,size(ImageIn,2)];
index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="XCorrBounds"), varargin, 'UniformOutput', 1));
if ~isempty(index)
    b=varargin{index+1};
end

PixelCutoff=300;
index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="PixelCutoff"), varargin, 'UniformOutput', 1));
if ~isempty(index)
    PixelCutoff=varargin{index+1};
end


for j=2:size(ImageIn,3)
    Cor=xcorr2(uint8(ImageIn(b(1):b(2),b(3):b(4),1)>PixelCutoff),uint8(ImageIn(b(1):b(2),b(3):b(4),j)>PixelCutoff));
    tmp=find(Cor(:)==max(Cor(:)));
    [y,x]=ind2sub(size(Cor),tmp(1));
    offset=[-((b(4)-b(3)+1-x)),-(b(2)-b(1)+1-y)];
%    figure(8)
%    imshowpair(ImageIn(b(1):b(2),b(3):b(4),1)>PixelCutoff,ImageIn(b(1):b(2),b(3):b(4),j)>PixelCutoff)
%    figure(9)
%    imshowpair(ImageIn(b(1):b(2),b(3):b(4),1)>PixelCutoff,imtranslate(ImageIn(b(1):b(2),b(3):b(4),j)>PixelCutoff,offset))
    ImageOut(:,:,j)=uint16(imtranslate(ImageIn(:,:,j),offset));    
end