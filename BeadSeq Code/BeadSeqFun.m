function [Bead BeadImage]=BeadSeqFun(BaseName, suffix, OutputFolder,BeadZeroThreshold,BarcodeSequence,NumPar)
%in addition to its outputs, the function will write images of the
%basecalls and of the microscope images to the OutputFolder directory


%clear all
%close all

%Currently, backgroundthreshold is determined *manually* from the images.
%This should be the average intensity of the dark pixels *between beads*.
%For the confocal 150 is usually good. Thanks, Marvin.
    starttime=tic();
backgroundthreshold=[150, 150, 150, 150]; %However, this can break if there is too much background in any one image

calibrate=0;
calibratechannels=[13,14]; %typically these should be successive, and calibratechannels(1) should be odd, so you're looking at a first and second ligation
channelnum=4; %number of channels
numchannels=4;
BeadSizeThreshold=50; %beads are clusters of at least 40 identical barcodes with no zeros.
ZScoreThreshold=1;
Cy3TxRMixing=0.5; %we subtract this factor *Cy3 Z score from TxR Z score. NOTE: This is empirical! If you change parameters you have to change this. NOTE SHOULD CHANGE THIS BACK I THINK IT WAS 0.5!!
PreviousRoundMixing=0.4; %we subtract this factor *the previous round's Z scores from the current round's Z scores on Even ligations only
loverride=0; %This is the maximum value of l to use when processing the barcodes. If loverride=0, we just use l.
%LOVE RIDE!

%At the moment we use all 20 ligations for the analysis, but this might not
%be necessary.



l=0;
while true
    if exist([BaseName,pad(num2str(l+1),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif'],'file')
        l=l+1;
    else
        break
    end
end

info=imfinfo([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif']);

%FOR DEBUGGING: Reduce the size of the area to be analyzed here
ROI=[[1,info.Height];[1,info.Width]];

ROIHeight=ROI(1,2)-ROI(1,1)+1;
ROIWidth=ROI(2,2)-ROI(2,1)+1;

%puckimagefull=zeros(info.Height,info.Width,numchannels);
puckimage=zeros(ROIHeight,ROIWidth,numchannels);
pixelvals=zeros(ROIHeight*ROIWidth,numchannels);
puckzscores=zeros(ROIHeight,ROIWidth,numchannels);
PreviousRoundZScores=zeros(ROIHeight,ROIWidth,numchannels);
%maxpixelvals=zeros(ROIHeight*ROIWidth,l);
CertaintyMap=zeros(ROIHeight,ROIWidth,l);
MaxOfPuckImage=zeros(ROIHeight,ROIWidth,l);

pixelvalplot=zeros(20,4);
pixelzscoreplot=zeros(20,4);

if loverride>0
    l=loverride;
end

%% Loading in the data and calling the bases
for m=1:20
    m
    for k=1:numchannels
%        puckimagefull(:,:,k)=imread([BaseName,int2str(m),' channel ',int2str(k),' transform.tif']);
%        puckimage(:,:,k)=puckimagefull(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),k);
        puckimage(:,:,k)=imread([BaseName,pad(num2str(m),2,'left','0'),' channel ',int2str(k),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)});
        %We find the standard deviation and mean of the bottom 75% of pixel
        %values
        pixelvals(:,k)=sort(reshape(puckimage(:,:,k),ROIHeight*ROIWidth,1));
        pixelvalsabovethreshold=pixelvals(pixelvals(:,k)>backgroundthreshold(k),k);
        pixelvalplot(m,k)=puckimage(1935,1356,k);
        %Dark beads are still much brighter than the background. So we can
        %use a threshold here to isolate the puck pixels effectively.
        pixelvalsmidMean(k)=mean(pixelvalsabovethreshold(1:floor(0.66*length(pixelvalsabovethreshold)))); %to de-emphasize dark beads, we take the mean and stdev of the middle 50% of pixels
        pixelvalsmidStd(k)=std(pixelvalsabovethreshold(1:floor(0.66*length(pixelvalsabovethreshold))));    
        puckzscores(:,:,k)=(double(puckimage(:,:,k))-pixelvalsmidMean(k))/pixelvalsmidStd(k);
        
        pixelzscoreplot(m,k)=puckzscores(1935,1356,k);
    end
    if m==calibratechannels(1) && calibrate
        histogram(pixelvals(:,1))
        %Use this to determine what the mixing is between cy3 and TxR
        figure(2)
        scatter(reshape(puckzscores(:,:,2),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,3),ROIHeight*ROIWidth,1))
        figure(3)
        scatter(reshape(puckimage(:,:,4),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,4),ROIHeight*ROIWidth,1))
        figure(4)
        scatter(reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,4),ROIHeight*ROIWidth,1))
    end
    
    if m==calibratechannels(2) && calibrate %this is one way to look at phasing: compare the mixing between 488 and 647 in round 1 and round 2. Then compare with the background subtracted
        figure(5)
        scatter(reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,4),ROIHeight*ROIWidth,1))        
    end

    puckzscores(:,:,3)=puckzscores(:,:,3)-Cy3TxRMixing*(puckzscores(:,:,2)>2).*puckzscores(:,:,2); %We subtract Cy3 from TxR, only if the Cy3 signal is significant.

    if m==calibratechannels(1) && calibrate %Use this to look at the mixing between Cy3 and TxR after correction
        figure(7)
        scatter(reshape(puckzscores(:,:,2),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,3),ROIHeight*ROIWidth,1))
    end
    
    if m==calibratechannels(2) && calibrate
        figure(8) %This, along with figure 9, is intended to look at whether subtracting off some fraction of the previous round's Z score would cause beads that are bright in both rounds not to be called correctly
        scatter(reshape(PreviousRoundZScores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1))
    end    
    
    if mod(m,2)==1
        PreviousRoundZScores=puckzscores;
    else
        SubtractLastRound=(max(puckzscores,[],3)./sum(puckzscores,3)<0.8) & max(puckzscores-PreviousRoundZScores,[],3)>2; %We only subtract out the last round's Z scores if there are other Z scores that are reasonably large compared to the brightest; and if at least once channel increased significantly in brightness
        for k=1:4 %this is ugly, but I couldn't work out a better way to do it
            puckzscores(:,:,k)=puckzscores(:,:,k)-PreviousRoundMixing*(PreviousRoundZScores(:,:,k)>2 & SubtractLastRound).*PreviousRoundZScores(:,:,k);
        end        
    end
    
    if m==calibratechannels(2) && calibrate %this is one way to look at phasing: compare the mixing between 488 and 647 in round 1 and round 2. Then compare with the background subtracted
        figure(6)
        scatter(reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,4),ROIHeight*ROIWidth,1))        
        figure(9)
        scatter(reshape(PreviousRoundZScores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1))
        
    end
    

        
    [M,I]=max(puckzscores,[],3); %I is the matrix of indices.
    
    sortedZscores=sort(puckzscores,3);
    BottomZScorestDev=std(sortedZscores(:,:,1:3),[],3);
    CertaintyMap(:,:,m)=(sortedZscores(:,:,4)-sortedZscores(:,:,3))./BottomZScorestDev;
    MaxOfPuckImage(:,:,m)=max(puckimage,[],3);
%    CertaintyMap(:,:,m)= ;
%    maxpixelvals(:,m)=max(pixelvals,[],2);
    
    I(M<ZScoreThreshold)=0; %This is the z score threshold -- z scores less than ZScoreThreshold are treated as 0.
    Indices(:,:,m)=uint8(I);
    
%    if m==2 %to look at phasing, we want to see the correlation between the previous round's call and this round's call. So we 
%        figure(8)
%        heatmap(reshape(Indices(:,:,1:2),ROIHeight*ROIWidth,2))
%    end

    %NOTE: This does not account for 1) knowledge of the previous base or
    %2) certainty, which can be incorporated based on the distance between
    %the two highest Z scores, or the Z score of the Z score, i.e. how much
    %higher the max Z score is than the other Z scores.
    
end

%% Show a plot of the pixel z scores for the pixel that is chosen in the previous section
if 0
    b=bar(pixelzscoreplot);
    b(1).FaceColor='b';
    b(2).FaceColor='g';
    b(3).FaceColor='y';
    b(4).FaceColor='r';
    
    b=bar(pixelvalplot);
    b(1).FaceColor='b';
    b(2).FaceColor='g';
    b(3).FaceColor='y';
    b(4).FaceColor='r';
    
    
end


%% Output images of the base calls:
if 1
%    OutputFolder='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\InputFolder-Puck85-170818\Position A1 Base Calls - Params 7\';
    %NOTE: if you are in matlab and you try to write a binary image as a
    %tiff, imagej won't be able to open it for some weird reason.
    imwrite(256*uint16(Indices(:,:,1)==1),[OutputFolder,'Channel1Calls.tiff']);
    imwrite(256*uint16(Indices(:,:,1)==2),[OutputFolder,'Channel2Calls.tiff']);
    imwrite(256*uint16(Indices(:,:,1)==3),[OutputFolder,'Channel3Calls.tiff']);
    imwrite(256*uint16(Indices(:,:,1)==4),[OutputFolder,'Channel4Calls.tiff']);
    for jm = 2:size(Indices,3)
        imwrite(256*uint16(Indices(:,:,jm)==1),[OutputFolder,'Channel1Calls.tiff'],'WriteMode','append');
        imwrite(256*uint16(Indices(:,:,jm)==2),[OutputFolder,'Channel2Calls.tiff'],'WriteMode','append');
        imwrite(256*uint16(Indices(:,:,jm)==3),[OutputFolder,'Channel3Calls.tiff'],'WriteMode','append');
        imwrite(256*uint16(Indices(:,:,jm)==4),[OutputFolder,'Channel4Calls.tiff'],'WriteMode','append');
    end

    if 0 %this is currently not executed because the maximum tiff file size is exceeded
    imwrite(imread([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel1Image.tiff']);
    imwrite(imread([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(2),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel2Image.tiff']);
    imwrite(imread([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(3),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel3Image.tiff']);
    imwrite(imread([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(4),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel4Image.tiff']);
    pause(2);
    for jm = 2:size(Indices,3)
        imwrite(imread([BaseName,pad(num2str(jm),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel1Image.tiff'],'WriteMode','append');
        imwrite(imread([BaseName,pad(num2str(jm),2,'left','0'),' channel ',int2str(2),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel2Image.tiff'],'WriteMode','append');
        imwrite(imread([BaseName,pad(num2str(jm),2,'left','0'),' channel ',int2str(3),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel3Image.tiff'],'WriteMode','append');
        imwrite(imread([BaseName,pad(num2str(jm),2,'left','0'),' channel ',int2str(4),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel4Image.tiff'],'WriteMode','append');
        pause(2);
    end
    end

end


%% Some analysis of the certainty:
if 0
RandomIndices=ceil(10752*10752*rand(2000));
ReshapedCertaintyMap=reshape(CertaintyMap,ROIHeight*ROIWidth,l);
scatter(ReshapedCertaintyMap(RandomIndices,1),ReshapedCertaintyMap(RandomIndices,3));
end

%scatter(maxpixelvals(RandomIndices,1),maxpixelvals(RandomIndices,3));

if 1 %this is the part where we call barcodes. We need to tell it which ligations to use for barcode calling.

%To identify bead locations, we have to convert Indices into integer barcodes.
%We use base 6. I.e. the barcode '142342' gets converted to 2*6^0 + 4*6^1 +
%3*6^2 + etc. Note that 0 indicates that no base was read.

    
%% Calling barcodes from the bases called above
FlattenedBarcodes=uint64(zeros(ROIHeight,ROIWidth));
NumSkippedBases=0;
for mm=1:l
    if BarcodeSequence(mm)==0
        NumSkippedBases=NumSkippedBases+1;
        continue
    end
    m=BarcodeSequence(mm);
    FlattenedBarcodes=FlattenedBarcodes+uint64(Indices(:,:,mm))*5^(m-1);
end
PresentBarcodes=unique(FlattenedBarcodes);

BarcodeOccCounts=histogram(FlattenedBarcodes,PresentBarcodes);
manypixelbarcodeswithzeros=PresentBarcodes(BarcodeOccCounts.Values>BeadSizeThreshold);

figure(8)
histogram(BarcodeOccCounts.Values,0:2:500)
axis([0,500,0,1])
axis 'auto y'
set(gca,'yscale','log')
title('Pixels per Barcode')
export_fig([OutputFolder,'Report.pdf'],'-append');    %export_fig([OutputFolder,'Report.pdf'],'-append');

%This is for the analysis of zeros:
BaseBalanceBarcodes=manypixelbarcodeswithzeros(cellfun(@numel,strfind(string(dec2base(manypixelbarcodeswithzeros,5,(l-NumSkippedBases))),'0'))<=7);
%The base 5 representations of the basecalls are:
BaseBalanceBase5Barcodes=cellfun(@(x) reverse(string(x)),{dec2base(BaseBalanceBarcodes,5,(l-NumSkippedBases))},'UniformOutput',false);
BaseBalanceBase5Barcodes=BaseBalanceBase5Barcodes{1};

BaseBalanceMatrix=zeros(5,l-NumSkippedBases);
for jp=1:(l-NumSkippedBases)
    testcmp0(jp)='0';
    testcmp1(jp)='1';
    testcmp2(jp)='2';
    testcmp3(jp)='3';
    testcmp4(jp)='4';
end
BaseBalanceMatrix(1,:)=sum(char(BaseBalanceBase5Barcodes)==testcmp0,1);
BaseBalanceMatrix(2,:)=sum(char(BaseBalanceBase5Barcodes)==testcmp1,1);
BaseBalanceMatrix(3,:)=sum(char(BaseBalanceBase5Barcodes)==testcmp2,1);
BaseBalanceMatrix(4,:)=sum(char(BaseBalanceBase5Barcodes)==testcmp3,1);
BaseBalanceMatrix(5,:)=sum(char(BaseBalanceBase5Barcodes)==testcmp4,1);
figure(77)
bar(BaseBalanceMatrix')
title('For Barcodes with 7 nonzero entries, the base balance per ligation');
export_fig([OutputFolder,'Report.pdf'],'-append');
%export_fig([OutputFolder,'Report.pdf'],'-append');

NumZerosPlot=zeros(1,l-NumSkippedBases+1);
for kl=1:(l-NumSkippedBases+1)
    NumZerosPlot(kl)=size(manypixelbarcodeswithzeros(cellfun(@numel,strfind(string(dec2base(manypixelbarcodeswithzeros,5,(l-NumSkippedBases))),'0'))==kl-1),1);
end
    figure(76);
    bar(0:(l-NumSkippedBases),NumZerosPlot)
    title('Number of barcodes with a given number of 0s in them');
    export_fig([OutputFolder,'Report.pdf'],'-append');    %export_fig([OutputFolder,'Report.pdf'],'-append');

%OLD CODE from the first analysis of Puck 8-5:
%nonzerobarcodes=PresentBarcodes(cellfun(@numel,strfind(string(dec2base(PresentBarcodes,5,(l-NumSkippedBases))),'0'))<=BeadZeroThreshold);
%beadhist=histogram(FlattenedBarcodes,nonzerobarcodes);% THIS DOESN'T WORK***
%manypixelbarcodestmp=nonzerobarcodes(beadhist.Values>BeadSizeThreshold);


%We make sure the length of manypixelbarcodes is divisible by NumPar to
%facilitate parallelization
manypixelbarcodestmp=manypixelbarcodeswithzeros(cellfun(@numel,strfind(string(dec2base(manypixelbarcodeswithzeros,5,(l-NumSkippedBases))),'0'))<=BeadZeroThreshold);
manypixelbarcodes=zeros(1,ceil(length(manypixelbarcodestmp)/NumPar)*NumPar);
manypixelbarcodes(1:length(manypixelbarcodestmp))=manypixelbarcodestmp; %We have to be careful here to get the indexing right.
%This is almost working but manypixelbarcodes still contains an element
%with 29 pixels rather than 30
disp(['There are ',num2str(length(manypixelbarcodestmp)),' barcodes passing filter.'])


%% Identifying which beads have significant clusters

%To parallelize, we break BeadImage up into a 3D array, with ParNum slices
%in the 3rd dimension. We reshape manypixelbarcodes so that it is a 2D
%array with ParNum slices in the 2nd dimension. And we likewise make
%BeadBarcodes and BeadLocations be 2d arrays. Then each worker gets its own
%row in BeadBarcodes, BeadLocations, and manypixelbarcodes, and runs the
%algorithm. At the end, we sum beadimage over the 3rd dimension and
%reshape BeadBarcodes and BeadLocations
BeadImage=false(ROIHeight,ROIWidth);
totalbarcodes=0;
delete(gcp('nocreate'));
pool=parpool(NumPar);

manypixelbarcodesforpar=reshape(manypixelbarcodes,ceil(length(manypixelbarcodes)/NumPar),NumPar);
BeadBarcodeCell={};
BeadLocationCell={};
BeadPixCell={};

parfor parnum=1:NumPar
    pp=0;
    LocalManyPixelBarcodes=manypixelbarcodesforpar(:,parnum);
    LocalBeadImage=false(ROIHeight,ROIWidth);
    LocalBeadBarcodes=zeros(1,nnz(manypixelbarcodes(:,parnum)));
    LocalBeadLocations=zeros(2,nnz(manypixelbarcodes(:,parnum)));
    LocalBeadPix=cell(1,nnz(manypixelbarcodes(:,parnum)));
    for qq=1:nnz(manypixelbarcodesforpar(:,parnum))
        if qq/100==floor(qq/100)
            disp(['Worker ',num2str(parnum),' is on barcode ',num2str(qq)])
        end
        connected=bwconncomp(FlattenedBarcodes==LocalManyPixelBarcodes(qq));
        centroids=regionprops(connected,'Centroid');
    if max(cellfun(@numel,connected.PixelIdxList))>BeadSizeThreshold
        for t=1:length(connected.PixelIdxList)
            if numel(connected.PixelIdxList{t})>BeadSizeThreshold
                LocalBeadImage(connected.PixelIdxList{t})=true; 
                pp=pp+1;
                LocalBeadBarcodes(pp)=LocalManyPixelBarcodes(qq);
                LocalBeadLocations(:,pp)=centroids(t).Centroid;
                LocalBeadPix{pp}=connected.PixelIdxList{t};
            end
        end
    end
    end
    imwrite(LocalBeadImage,[OutputFolder,'Worker_',num2str(parnum),'_LocalBeadImage.tif']);
    BeadBarcodeCell{parnum}=LocalBeadBarcodes;
    BeadLocationCell{parnum}=LocalBeadLocations;
    BeadPixCell{parnum}=LocalBeadPix;
end

BeadBarcodeLength=0;
for k=1:NumPar
    BeadBarcodeLength=BeadBarcodeLength+length(BeadBarcodeCell{k});
end
BeadBarcodes=zeros(1,BeadBarcodeLength);
BeadLocations=zeros(2,BeadBarcodeLength);
BeadPixCelljoined=cell(1,BeadBarcodeLength);

delete(pool);

BeadBarcodeIndex=1;
for k=1:NumPar
    BeadBarcodes(BeadBarcodeIndex:(BeadBarcodeIndex+length(BeadBarcodeCell{k})-1))=BeadBarcodeCell{k};
    BeadLocations(:,BeadBarcodeIndex:(BeadBarcodeIndex+length(BeadBarcodeCell{k})-1))=BeadLocationCell{k};
    tmpcell=BeadPixCell{k};
    for ll=1:length(BeadBarcodeCell{k})
        BeadPixCelljoined{BeadBarcodeIndex+ll-1}=tmpcell{ll}; %if this is too long, we could also just make the beadpix cell array within the parfor above
    end
    BeadBarcodeIndex=BeadBarcodeIndex+length(BeadBarcodeCell{k});
    BeadImage=BeadImage + imread([OutputFolder,'Worker_',num2str(k),'_LocalBeadImage.tif']);
end
Bead=struct('Barcodes',num2cell(BeadBarcodes),'Locations',num2cell(BeadLocations,1),'Pixels',BeadPixCelljoined);


    imwrite(BeadImage,[OutputFolder,'BeadImage.tif']);


%% This is the non-parallel version
if 0
BeadImage=false(ROIHeight,ROIWidth);
BeadBarcodes=zeros(1,length(manypixelbarcodes));
BeadLocations=zeros(2,length(manypixelbarcodes));
pp=0;
totalbarcodes=0;

for qq=1:length(manypixelbarcodes)
    if qq/100==floor(qq/100)
        qq
    end
    connected=bwconncomp(FlattenedBarcodes==manypixelbarcodes(qq));
    centroids=regionprops(connected,'Centroid');
    if max(cellfun(@numel,connected.PixelIdxList))>BeadSizeThreshold
        for t=1:length(connected.PixelIdxList)
            if numel(connected.PixelIdxList{t})>BeadSizeThreshold
                BeadImage(connected.PixelIdxList{t})=true;
                pp=pp+1;
                BeadBarcodes(pp)=manypixelbarcodes(qq);
                BeadLocations(:,pp)=centroids(t).Centroid;
            end
        end
    end
end
%for q=1:length(BeadBarcodes)
%    if q/100==floor(q/100)
%        q
%    end
%    BeadImage=BeadImage | (FlattenedBarcodes==BeadBarcodes(q)); %this is not exactly correct.
%end
end
figure(78)
imshow(BeadImage)

save([OutputFolder,'AnalysisOutputs-selected'],'BeadImage','FlattenedBarcodes','Indices','Bead','BaseBalanceMatrix','NumZerosPlot','-v7.3');


%export_fig([OutputFolder,'\Report.pdf'],'-append');
fileid=fopen([OutputFolder,'Metrics.txt'],'w');
fprintf(fileid,['The total runtime for basecalling was ',num2str(toc(starttime)/60),' minutes.\n',...
    'There are ',num2str(length(manypixelbarcodestmp)),' barcodes passing filter.\n',...
    'BeadSizeThreshold=',num2str(BeadSizeThreshold),'.\n',...
    'ZScoreThreshold=',num2str(ZScoreThreshold),'.\n',...
    'Cy3TxRMixing=',num2str(Cy3TxRMixing),'.\n',...
    'PreviousRoundMixing=',num2str(PreviousRoundMixing),'.\n',...
    'BeadZeroThreshold=',num2str(BeadZeroThreshold)...
    ]);
fclose(fileid);

end