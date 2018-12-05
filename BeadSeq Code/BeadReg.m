%As input, we want a list of files that are of the form X1, X2, X3, etc.,
%where X is the BaseNam.

%CURRENTLY: it's not working here, but it's also not working if I save the
%images and then use the original find_roi on them. So something is wrong
%about the way I am treating the images.
clear all
close all;
clearvars -except f1s d1s database;
cd find_roi
run('helpers\vlfeat-0.9.20\toolbox\vl_setup.m');
addpath helpers;
cd ..

channelnum=4; %number of channels

BaseName='rSAP - Puck 1 - Base ';

l=0;
while true
    if exist([BaseName,int2str(l+1),'.tif'],'file')
        l=l+1;
    else
        break
    end
end

%Take the max project of the first image, across channels
path = [BaseName,'1','.tif'];

info=imfinfo(path);
stackheight=numel(info);
imwidth=info.Width;
imwidth=imwidth(1);
imheight=info.Height;
imheight=imheight(1);
image=zeros(imheight,imwidth,channelnum);
for j=1:channelnum
    image(:,:,j)=imread(path,'index',j);
            %WE CAN NORMALIZE THE IMAGE HERE IF NECESSARY, to facilitate
            %feature finding
end
firstimagemax=max(image,[],3);

registeredimages=zeros(imheight,imwidth,channelnum,l);

for i=1:l
    %maxproject across channels
        info=imfinfo(path);
        stackheight=numel(info);
        imwidth=info.Width;
        imwidth=imwidth(1);
        imheight=info.Height;
        imheight=imheight(1);
        image=zeros(imheight,imwidth,channelnum);
        for j=1:channelnum
            image(:,:,j)=imread(path,'index',j);
            %WE CAN NORMALIZE THE IMAGE HERE IF NECESSARY, to facilitate
            %feature finding
        end
        imagemax=max(image,[],3);
        [tform, mask]=find_roi_new(imagemax,firstimagemax,0);
        for k=1:channelnum
            registeredimages(:,:,k,i)=imtransform(image(:,:,k),tform);
        end
    
end


%[cropped_roi] = imtransform(cropped_roi,tform);
%mask = imtransform(mask,tform);

%we output a 4D tif, where the first two dimensions are spatial, the third
%dimension is the color, and the 4th dimension is the imaging cycle.