%Registration
cd 'C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code'

clear all
close all

BaseName='find_roi\InputFolder-Puck85-170818\Puck 85 Ligation ';
suffix=' Position A1'; %this is '' by default

%Currently, backgroundthreshold is determined *manually* from the images.
%This should be the average intensity of the dark pixels *between beads*.
%For the confocal 150 is usually good. Thanks, Marvin.
backgroundthreshold=[150, 150, 150, 150]; %However, this can break if there is too much background in any one image

calibrate=0;
channelnum=4; %number of channels
numchannels=4;
BeadZeroThreshold=0; %Number of 0s allowed in a barcode
BeadSizeThreshold=20; %beads are clusters of at least 20 identical barcodes with no zeros.
ZScoreThreshold=1;
Cy3TxRMixing=0.2; %we subtract this factor *Cy3 Z score from TxR Z score. NOTE: This is empirical! If you change parameters you have to change this. NOTE SHOULD CHANGE THIS BACK I THINK IT WAS 0.5!!
PreviousRoundMixing=0.4; %we subtract this factor *the previous round's Z scores from the current round's Z scores on Even ligations only
loverride=0; %This is the maximum value of l to use when processing the barcodes. If loverride=0, we just use l.
%LOVE RIDE!

%At the moment we use all 20 ligations for the analysis, but this might not
%be necessary.



l=0;
while true
    if exist([BaseName,int2str(l+1),' channel ',int2str(1),suffix,' transform.tif'],'file')
        l=l+1;
    else
        break
    end
end

info=imfinfo([BaseName,int2str(1),' channel ',int2str(1),suffix,' transform.tif']);

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
for m=13:14
    m
    for k=1:numchannels
%        puckimagefull(:,:,k)=imread([BaseName,int2str(m),' channel ',int2str(k),' transform.tif']);
%        puckimage(:,:,k)=puckimagefull(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),k);
        puckimage(:,:,k)=imread([BaseName,int2str(m),' channel ',int2str(k),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)});
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
    if m==1 && calibrate
        histogram(pixelvals(:,1))
        %Use this to determine what the mixing is between cy3 and TxR
        figure(2)
        scatter(reshape(puckzscores(:,:,2),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,3),ROIHeight*ROIWidth,1))
        figure(3)
        scatter(reshape(puckimage(:,:,4),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,4),ROIHeight*ROIWidth,1))
        figure(4)
        scatter(reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,4),ROIHeight*ROIWidth,1))
    end
    
    if m==2 && calibrate %this is one way to look at phasing: compare the mixing between 488 and 647 in round 1 and round 2. Then compare with the background subtracted
        figure(5)
        scatter(reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,4),ROIHeight*ROIWidth,1))        
    end

    puckzscores(:,:,3)=puckzscores(:,:,3)-Cy3TxRMixing*(puckzscores(:,:,2)>2).*puckzscores(:,:,2); %We subtract Cy3 from TxR, only if the Cy3 signal is significant.

    if m==1 && calibrate %Use this to look at the mixing between Cy3 and TxR after correction
        figure(7)
        scatter(reshape(puckzscores(:,:,2),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,3),ROIHeight*ROIWidth,1))
    end
    
    if m==2 && calibrate
        figure(8) %This, along with figure 9, is intended to look at whether subtracting off some fraction of the previous round's Z score would cause beads that are bright in both rounds not to be called correctly
        scatter(reshape(PreviousRoundZScores(:,:,1),ROIHeight*ROIWidth,1),reshape(puckzscores(:,:,1),ROIHeight*ROIWidth,1))
    end    
    
    if mod(m,2)==1
        PreviousRoundZScores=puckzscores;
    else
        puckzscores=puckzscores-PreviousRoundMixing*(PreviousRoundZScores>2).*PreviousRoundZScores;
        
    end
    
    if m==2 %this is one way to look at phasing: compare the mixing between 488 and 647 in round 1 and round 2. Then compare with the background subtracted
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
if 0
    OutputFolder='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\InputFolder-Puck85-170818\Position A1 Base Calls - Params 3\';
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


    imwrite(imread([BaseName,int2str(1),' channel ',int2str(1),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel1Image.tiff']);
    imwrite(imread([BaseName,int2str(1),' channel ',int2str(2),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel2Image.tiff']);
    imwrite(imread([BaseName,int2str(1),' channel ',int2str(3),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel3Image.tiff']);
    imwrite(imread([BaseName,int2str(1),' channel ',int2str(4),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel4Image.tiff']);
    pause(2);
    for jm = 2:size(Indices,3)
        imwrite(imread([BaseName,int2str(jm),' channel ',int2str(1),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel1Image.tiff'],'WriteMode','append');
        imwrite(imread([BaseName,int2str(jm),' channel ',int2str(2),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel2Image.tiff'],'WriteMode','append');
        imwrite(imread([BaseName,int2str(jm),' channel ',int2str(3),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel3Image.tiff'],'WriteMode','append');
        imwrite(imread([BaseName,int2str(jm),' channel ',int2str(4),suffix,' transform.tif'],'PixelRegion',{ROI(1,1:2),ROI(2,1:2)}),[OutputFolder,'Channel4Image.tiff'],'WriteMode','append');
        pause(2);
    end

end


%% Some analysis of the certainty:
if 0
RandomIndices=ceil(10752*10752*rand(2000));
ReshapedCertaintyMap=reshape(CertaintyMap,ROIHeight*ROIWidth,l);
scatter(ReshapedCertaintyMap(RandomIndices,1),ReshapedCertaintyMap(RandomIndices,3));
end

%scatter(maxpixelvals(RandomIndices,1),maxpixelvals(RandomIndices,3));

if 0 %this is the part where we call barcodes. We need to tell it which ligations to use for barcode calling.

%To identify bead locations, we have to convert Indices into integer barcodes.
%We use base 6. I.e. the barcode '142342' gets converted to 2*6^0 + 4*6^1 +
%3*6^2 + etc. Note that 0 indicates that no base was read.

    
%Ultimately we are reading 12 bases and 6^12 ~ 2^31, so FlattenedBarcodes will be uint32.
%% Calling barcodes from the bases called above
FlattenedBarcodes=uint32(zeros(ROIHeight,ROIWidth));
for m=1:l
    FlattenedBarcodes=FlattenedBarcodes+uint32(Indices(:,:,m))*5^(m-1);
end
PresentBarcodes=unique(FlattenedBarcodes);
%This is if you only want ones with no '0's:
%nonzerobarcodes=PresentBarcodes(find(~contains(string(dec2base(PresentBarcodes,5,l)),'0')));
%if you're using a threshold, use:
nonzerobarcodes=PresentBarcodes(cellfun(@numel,strfind(string(dec2base(PresentBarcodes,5,l)),'0'))<=BeadZeroThreshold);

beadhist=histogram(FlattenedBarcodes,nonzerobarcodes);

manypixelbarcodes=nonzerobarcodes(beadhist.Values>BeadSizeThreshold); %We have to be careful here to get the indexing right.
%This is almost working but manypixelbarcodes still contains an element



%% Identifying which beads have significant clusters
%with 29 pixels rather than 30.
BeadImage=false(ROIHeight,ROIWidth);
BeadBarcodes=[];
pp=0;
totalbarcodes=0;
for qq=1:length(manypixelbarcodes)
    if qq/100==floor(qq/100)
        qq
    end
    connected=bwconncomp(FlattenedBarcodes==manypixelbarcodes(qq));
    if max(cellfun(@numel,connected.PixelIdxList))>BeadSizeThreshold
        pp=pp+1;
        BeadBarcodes(pp)=manypixelbarcodes(qq);
        for t=1:length(connected.PixelIdxList)
            if numel(connected.PixelIdxList{t})>BeadSizeThreshold
                BeadImage(connected.PixelIdxList{t})=true;
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
imshow(BeadImage)

save([BaseName,'AnalysisOutputs-selected-Ligation9'],'BeadImage','FlattenedBarcodes','Indices','BeadBarcodes','CertaintyMap','-v7.3');
%% This was the rastering version of the puck identification module
%We are now going to try to identify the barcodes that have beads on them.
%We go through on a pixel by pixel basis:

if 0
SignificantBarcodes=[];
pp=0;
totalbarcodes=0;
BeadSizeThreshold=30;
    
barcodesseen=0;
barcodeschecked=[];
for pxa=1:ROIHeight
    pxa
    for pxb=1:ROIWidth
        if ~contains(dec2base(FlattenedBarcodes(pxa,pxb),5,l),'0') && isempty(find(barcodeschecked==FlattenedBarcodes(pxa,pxb),1)) %Note that conditions are tested left to right
            barcodesseen=barcodesseen+1;
            barcodeschecked(barcodesseen)=FlattenedBarcodes(pxa,pxb);
            samebarcodepix=FlattenedBarcodes==FlattenedBarcodes(pxa,pxb);
             if length(find(samebarcodepix))>BeadSizeThreshold
                connected=bwconncomp(samebarcodepix);
                if max(cellfun(@numel,connected.PixelIdxList))>BeadSizeThreshold
                    pp=pp+1;
                    SignificantBarcodes(pp)=FlattenedBarcodes(pxa,pxb);
                end
            end
        end
    end
end

end
        %We should also keep track of the biggest clusters here, so we can
        %save them and look at them later.
end




%% This is code that is good if you have more barcodes than beads
if 0
%Once we have more barcodes than beads, we will want to go through on a
%pixel by pixel basis. For each pixel, if it has a 0 in its barcode, or if you have seen it before, skip
%it. Otherwise, look for all other pixels with the same barcodes, and if
%you get bead-sized clusters save them and record the barcode.
for j=1:5^l %this is now the complex part
    if j/100==floor(j/100)
        j
        pp
        totalbarcodes-pp
    end
    if ~isempty(strfind(dec2base(j,5,l),'0'))
        continue
    end
    totalbarcodes=totalbarcodes+1;
    %We first break j into its base 5 representation. If any bits are 0 we
    %skip it.
    connected=bwconncomp(FlattenedBarcodes==j);

    if max(cellfun(@numel,connected.PixelIdxList))>35;
        pp=pp+1;
        SignificantBarcodes(pp)=j;
        %We should also keep track of the biggest clusters here, so we can
        %save them and look at them later.
    end

end
end
    %SOMEHOW: note, every barcode has very large numbers of isolated
    %pixels

%    for d=1:5^l
%        count1=count1+length(find(FlattenedBarcodes==d));
%    end
    
%    figure(1)
%    clf
%    imshow(I==1 & M>3)
%    figure(2)
%    clf
%    imshow(I==2 & M>3)
    
    %For each pixel we find the z score

    
%we first bin the image by a factor of x
%    for j=1:numchannels
%        BinnedImage(:,:,j)=imresize(puckimage(:,:,j),1/binningnumber,'method','bicubic');
%        %BinnedImage(:,:,k)=puckimage(1:2:end,1:2:end,k)+puckimage(1:2:end,1:2:end+1)+puckimage(1:2:end+1,1:2:end)+puckimage(1:2:end+1,1:2:end+1);
%    end

    
%For each time step, we first fit a mixture of Gaussians model to the image
%to identify bright pixels and dim pixels
%We then use the clustering algorithm to identify clusters of bright pixels
%in each channel in each time step.
%We then take the intersection of the clusters.



%An alternative method would be the following: for
%each pixel, we assemble the intensities in each channel at each round. For
%each channel, we take the minimum e.g. 4 intensities and get a standard deviation
%and mean from that. (So this won't work on barcodes which consist of 9 of the same base,
%of which there are 4 (choose the base that gets repeated) * 4^3 * 12 choose 3 = 56329, but there are 16777216
%barcodes total, so it will fail on 0.3% of barcodes, which is a lot... but
%it's okay for now. We might have to increase the number of barcodes. Note
%that if only use the bottom 3 intensities the number of barcodes we fail
%on goes to 4224 = 0.025%

%This will actually not be so terrible, since we can do it all with array
%mathematics.

%We then identify over each round which channel is the most standard deviations
%away from the mean for its channel, and throw and error if two channels are both
%very bright or if no channel is bright. We now have a barcode for each
%pixel, which we can represent in base 4, or base 10 whatever, and we can go through and use the clustering algorithm on each
%barcode. I'm not exactly sure how to do the latter step efficiently, but
%there is probably a way. It should be a special case of the bwconncomp
%function -- maybe there is a version of that function that includes
%colors.

%Ideally we would have a different 2D array for each channel, with the
%first dimension being pixel number and the 2nd dimension being intensity
%at each round. We then also keep a struct somewhere with a mapping between
%the linear index of the pixel and its index in 2d space.

%PixelVals488=zeros(NumPixels,NumRows);
%PixelValsCy3=zeros(NumPixels,NumRows);
%PixelValsTexR=zeros(NumPixels,NumRows);
%PixelVals647=zeros(NumPixels,NumRows);
%ALLOCATE MEMORY BEFOREHAND AND CLEAR OUT ARRAYS AS WE GO, OR WE WILL
%PROBABLY RUN OUT

%Assemble the PixelVals arrays, so that PixelValsX(i,j) is the intensity in
%channel X at pixel i at round j.
%for i=j:NumRounds
%    PixelVals488(:,j)=;
%    PixelValsCy3(:,j)=;
%    PixelValsTexR(:,j)=;
%    PixelVals647(:,j)=;
%end

%We sort the arrays. We should allocate this memory beforehand...
%PixelValsSorted488=sort(PixelVals488,2); %along the
%PixelValsLowMeans488=mean(PixelVals488(:,1:4),2); %take the mean along the rows
%PixelValsLowStd488=std(PixelVals488(:,1:4),2); %take the mean along the rows

    
%PixelValsSortedCy3=sort(PixelVals488,2); %along the
%PixelValsLowMeansCy3=mean(PixelVals488(:,1:4),2); %take the mean along the rows
%PixelValsLowStdCy3=std(PixelVals488(:,1:4),2); %take the mean along the rows

%ZScores=zeros(NumPixels,4);
%ZScoresSorted=zeros(NumPixels,4);
%ZScoresLowStd=zeros(NumPixels,4);
%for jj=1:NumRounds
%    ZScores(:,1)=(PixelVals488(:,jj)-PixelValsLowMeans488)/PixelValsLowStd488;%check that this is elementwise division
    %2 is Cy3, 3 is TexR, 4 is 647
    
%    PixelValBaseCalling(:,jj)=%we want the index of the max element along each row
    
%    ZScoresSorted=sort(ZScores,2);
%    ZScoresLowStd=std(ZScoresSorted(:,1:3),2);
%    Certainty(:,jj)=(ZScoresSorted(:,4)-ZScoresSorted(:,3))/ZScoresLowStd; %We want to define a metric of certainty

%%This is the null model for what the beads should look like if they are
%randomly distributed
%test73=zeros(256,256);
%for lm=1:256
%    for no=1:256
%        test73(lm,no)=floor(4*rand());
%    end
%end
%imagesc(test73==1)