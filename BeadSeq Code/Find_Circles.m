cd 'C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code'
%green1 = imread('TestImage10x.tif','index',1);

green1 = imread('TestImage10xpng.png');

[centers,radii,metric]=imfindcircles(green1,[5, 16],'ObjectPolarity','bright','Sensitivity',0.93);

green2 = imread('TestImage10xpng.png');

[centers2,radii2,metric2]=imfindcircles(green1,[5, 16],'ObjectPolarity','bright','Sensitivity',0.93);

imshow(green1)
hold on;
viscircles(centers,radii,'color','b')
hold off;
%Read: https://www.mathworks.com/help/images/ref/imfindcircles.html
%Read: https://www.mathworks.com/help/images/examples/detect-and-measure-circular-objects-in-an-image.html