clear all
cd 'C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\Data\170405 Phasing'


%path='First base - normal 10x.tif'
%path='10x beads control.tiff';
%path='10x beads CIP only.tif';
path='10x beads dark oligo only.tif';
%path='10x beads dark oligo and CIP.tif';
%camerapixelx=2048;
%camerapixely=2048;
comparisoncolors=[1,3];

info2=imfinfo(path);
stackheight=numel(info2);
imwidth=info2.Width;
imwidth=imwidth(1);
imheight=info2.Height;
imheight=imheight(1);

images=zeros(imheight,imwidth,stackheight);


for i=1:stackheight
    images(:,:,i)=imread(path,'index',i);
end
x=squeeze(reshape(images(:,:,comparisoncolors(1)),[1,imheight*imwidth]));
y=squeeze(reshape(images(:,:,comparisoncolors(2)),[1,imheight*imwidth]));

interestingx=x(x>(mean(x)));
interestingy=y(x>(mean(x)));
randomnumbers=floor(length(interestingx)*rand(1,50000));
randomx=interestingx(randomnumbers);
randomy=interestingy(randomnumbers);


scatter(randomx,randomy);