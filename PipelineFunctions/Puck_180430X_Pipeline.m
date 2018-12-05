%%%%This pipeline function assumes that we are using the large image
%%%%feature and the XY feature in Nikon images, and that files ending in
%%%%xy1 are for puck 1. We do not do stitching.

%This function is for beads with 7 J bases in two groups.

%%%%SETUP:

%0) Make sure ImageJ is closed.

%1) When you are done with the code, you should make the raw data folder in
%Pucks and the InputFolder in find_roi online only via smart sync to save
%hard drive space. You should also delete the pucktmp directory.

%2) A number of paths are hardcoded in the code currently, e.g.
%"C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\helpers\vlfeat-0.9.20\toolbox\"
%in find_roi_stack_fun

%3) Change ImageSize to reflect something slightly smaller than the size of
%the final stitched images.

%4) If you have a digital gene expression matrix, you need to take the entire CSV and put it in the Illumina folder (C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\Puck_180106_3) and call it DGE.csv. I then copy out the top row with notepad++, delete the first entry (which is there as a placeholder for the row labels) and put it into a different file called IlluminaBarcodes.txt

%5) Each ligation sequence and primer set will require a bunch of modifications
%to MapLocationsFun. To accommodate for this, I make a different version of MapLocationsFun for
%each ligation sequence. You have to change which version is called below.

%6) NOTE: SlideseqSave.ijm is stored in C:\Fiji.app\macros, and the text
%is:
%params=getArgument()
%print("Opening file:")
%print(params)
%open(params);
%saveAs("Tiff","C:\\PuckOutputTmp.tif")
%print("Done.")
%exit();
%% Initialize
%Set up the paths to all of the various folders with functions in them that
%we will call:
clear all
close all

addpath('C:\Fiji.app\scripts','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\PipelineFunctions');
javaaddpath('C:\Program Files\MATLAB\R2017a\java\mij.jar');
javaaddpath(which('MatlabGarbageCollector.jar'))
%We assume that the nd2s have been exported to tiffs in the format:
%DescriptiveNametLYX, where DescriptiveName refers to one run of the microscope,  Y is a letter and X is a number, and L is the
%name of the ligation within that DescriptiveName file series.
%We convert the files to a final output format:
%Puck85 Ligation X Position AB
%BeadType="180402";
BeadType="ReversePhase";

SaveData=1;
IlluminaReadThreshold=10;
MaximumBarcodesToAnalyze=80000;
NumLigations=20; %We assume that the missing ligation is Truseq-4 if you are doing 13.
NumBases=14; %This is the number of bases sequenced
BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14];
%BarcodeSequence=[1,2,3,4,0,5,0,6,7,8,9,10,0,11,0,12,0,13]; %this determines how the numerical barcode is built from the ligations. Basically, the constant bases should be 0, and all other bases should be consecutive
FolderWithRawTiffs='D:\Slideseq\Raw\180522 - Pucks 180430-X\';
FolderWithProcessedTiffs='D:\Slideseq\Processed\';
tmpfolder='C:\Users\sgr\pucktmp\';
IndexFiles={'primers truseq through up','primers up-1 through up-4'}; %give the prefixes of each of the files. The next character after should be 't',
LigationToIndexFileMapping=[1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];%for ligations 1:20, which file number are the ligations found in?
%for ligations 1:20, which value of t are they, within their file? This
%could also be deduced from the LigationToIndexFileMapping array
tnumMapping=[1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8];
%Note that the order in PuckNames should match the order in the .nd2 file.
PuckNames={'Puck_180430_1','Puck_180430_2','Puck_180430_3','Puck_180430_4','Puck_180430_5','Puck_180430_6','Puck_180430_7','Puck_180430_8'}; %give the names of the pucks
OutputFolderRoot='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\';%\Puck_180106_1\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_2\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_3\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_4\'};
IlluminaFolderRoot='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\';
%SequencingFileName={'fei_S3_L001_R1_001.fastq'};
ImageSize=6030; %The image registration will output images that are not all the exact same size, because of stitching. So find_roi_stack_fun crops the images a bit. We can choose this value to be something like 0.95 * the size of the images. So e.g. for 3x3 it should be 0.95*2048*(2*(2/3)+1) = 4500. For 7x7 it is 0.95*2048*(1+6*2/3)=9728. I used 10400 previously though.
TopBarcodeNum=80000; %when it comes to analyzing barcodes from the illumina data, we only consider the TopBarcodeNum barcodes with the most reads. This should roughly be the number of beads on the puck.
BeadZeroThreshold=1;

%The illumina barcodes are assumed to be 13 bases long, consisting of the
%first 6 J bases and then the last 7 J bases. If the barcodes you are using
%are different, you have to change it in the bs2cs calls.

%And ligation sequnece:
%Tru_L1, Tru_L2, Tru-1_L1, Tru-1_L2, Tru-2_L1, Tru-2_L2, Tru-3_L1, Tru-3_L2, Tru-4_L1, Tru-4_L2
%The numbers in these variables are the ***second base being interrogated***
%Equivalently, it is the number associated with that ligation in figure 4a of the Fisseq nature protocols paper
PrimerNLigationSequence = [2, 7, 1, 6, 5, 4, 3]; %good for 14 ligations
%PrimerNLigationSequence = [2, 7, 1, 6, 5, 4]; %good for 13 ligations

%UP_L1, UP_L2, UP-1_L1, UP-1_L2, UP-2_L1, UP-2_L2, UP-3_L1, UP-3_L2, UP-4_L1, UP-4_L2

PrimerUPLigationSequence=[2, 7, 1, 6, 5, 4,3];

%Note that bs2cs will output the colors corresponding to each ligation in this ligation sequence, in the order specified here,
%So if you were only to do 6 ligations off primer N instead of 7, you would only need to
%remove the final element in the PrimerNLigationSequence

InverseLigationSequence=[3,1,7,6,5,4,2,10,8,14,13,12,11,9]; %Good for both 13 and 14 ligations.
%Before exporting, "N"s are added to the end of the ligation
%sequence, and must then be inserted into the correct location to stand in
%for the missing ligations. This code puts the N at position 7
%WhichLigationsAreMissing=[1,2,3,4,5,6,14,7,8,9,10,11,12,13];
WhichLigationsAreMissing=[1,2,3,4,5,6,7,8,9,10,11,12,13,14];

OutputFolders={};
for puck=1:length(PuckNames)
    ProcessedImageFolders{puck}=[FolderWithProcessedTiffs,PuckNames{puck},'\'];
    mkdir([FolderWithProcessedTiffs,PuckNames{puck}]);
end
for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
%We now send the command to do the registration with find_roi
    
    OutputFolders{puck}=[OutputFolderRoot,PuckNames{puck},'\'];
    mkdir(OutputFolders{puck});
end


display('Renaming Files');

%% Move files
for puck=1:length(PuckNames)
    mkdir([tmpfolder,PuckNames{puck}]);
    for ligation=1:NumLigations
        tnum=tnumMapping(ligation);
        if exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),'xy',num2str(puck),'.tif'])
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),'xy',num2str(puck),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),'xy',num2str(puck),'.tif'])
            filename = [FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),'xy',num2str(puck),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',num2str(tnum),'.tif'])
            filename = [FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',num2str(tnum),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',pad(num2str(tnum),2,'left','0'),'.tif'])
            filename = [FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',pad(num2str(tnum),2,'left','0'),'.tif'];
        else
            assert(1==0)
        end
            
            outputfilename=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'];
%                outputfilename=[OutputFolders{puck},PuckNames{puck},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Position_X_',pad(num2str(xpositionindex),2,'left','0'),'_Y_',pad(num2str(yposition),2,'left','0'),'.tif'];
        %command=['inputfile="',replace(filename,'\','\\'),'" outputfile="',replace(outputfilename,'\','\\'),'"'];
            commandfile=fopen('C:\FijiCommand.cmd','w');
            fwrite(commandfile,strcat('C:\Fiji.app\ImageJ-win64.exe --headless --console -macro SlideseqSave.ijm "',replace(filename,'\','\\'),'"'));
            fclose(commandfile);
            !C:/FijiCommand
            movefile('C:\PuckOutputTmp.tif',outputfilename)
    end
end


%% Registration
display('Image Registration')

for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
%We now send the command to do the registration with find_roi
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    Suffix='_Stitched';
    find_roi_stack_fun(BaseName,Suffix,ImageSize);
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end
%% Bead calling and sequencing
%we now run Bead_Seq itself. Again, where parallelization is trivial in
%Bead_Seq, it is implemented naively using parfor.

%If you have run the basecalling previously, and are rerunning it, you
%should rename that folder. For simplicity, we don't date the basecalling
%folder.

display('Base Calling')

%c=clock;
%starttimereadable=[num2str(c(2)),'-',num2str(c(3)),'_',pad(num2str(c(4)),2,'left','0'),pad(num2str(c(5)),2,'left','0')];

for puck=1:length(PuckNames)
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
	[Bead BeadImage]=BeadSeqFun2(BaseName,suffix,OutputFolders{puck},BeadZeroThreshold,BarcodeSequence,20,NumLigations,PuckNames{puck});
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end

%% Match Illumina barcodes - to rerun the bead mapping only, start here

display('Matching Illumina Barcodes')
c=clock;
mappingstarttimereadable=[num2str(c(2)),'-',num2str(c(3)),'_',pad(num2str(c(4)),2,'left','0'),pad(num2str(c(5)),2,'left','0')];

for puck=1:length(PuckNames)
    mkdir([OutputFolderRoot,PuckNames{puck},'\BeadMapping_',mappingstarttimereadable]);
    MappingOutputFolders{puck}=[OutputFolderRoot,PuckNames{puck},'\BeadMapping_',mappingstarttimereadable,'\'];

    %Load in the CSV:
    if ~(exist([IlluminaFolderRoot,PuckNames{puck},'\DGE.csv']))
        continue
    end
    DGETable=readtable([IlluminaFolderRoot,PuckNames{puck},'\DGE.csv'],'ReadVariableNames',true);

    %Special code if analyzing DGEs with the total transcript counts, rather than the breakdown by gene:
%    GeneNames=DGETable.Properties.VariableNames;
%    GeneNames=GeneNames(2:end);
%    if puck==1
%        IlluminaBarcodes=string(DGETable.CTATCATCCTCGN);
%    end
%    if puck==3
%        IlluminaBarcodes=string(DGETable.CCGTACCTACCTG);
%    end
%    if puck==2
%        IlluminaBarcodes=string(DGETable.CGGCAGCAGGAAC);
%    end
%    if puck==5
%        IlluminaBarcodes=string(DGETable.CTGCTTTTGCACT);
%    end
%    if puck==6
%        IlluminaBarcodes=string(DGETable.CCTCACCCCTTTC);
%    end
    %NORMAL CODE:
    IlluminaBarcodes=DGETable.Properties.VariableNames;
    IlluminaBarcodes=IlluminaBarcodes(2:end);
    try
        GeneNames=string(DGETable.GENE);
    catch
        GeneNames=string(DGETable.Var1); %this is legacy        
    end
    DGE=table2array(DGETable(:,2:end));
%    DGE=csvread([IlluminaFolderRoot,PuckNames{puck},'\DGE.csv'],1,1);
    %Special code if analyzing the total # transcripts only DGE.
%    DGE=DGE';

    DGEallsums=sum(DGE,1);
    IlluminaBarcodesToAnalyze=min(length(find(DGEallsums>=IlluminaReadThreshold)),MaximumBarcodesToAnalyze);

    [sorted,sortindices]=sort(DGEallsums,'descend');
    DGESorted=DGE(:,sortindices(1:IlluminaBarcodesToAnalyze));
    IlluminaBarcodesDownsampled=IlluminaBarcodes(sortindices(1:IlluminaBarcodesToAnalyze));
    
    %Because we are only sequencing J bases after Truseq and UP, we cut out
    %the 8th J base after truseq.
    if BeadType=="180402"
        IlluminaBarcodesDownsampled=cellfun(@(x) [x(1:7),x(9:15)],IlluminaBarcodesDownsampled,'UniformOutput',false);
    end
    
    load([OutputFolders{puck},'AnalysisOutputs-selected'],'Bead');
    if BeadType=="180402"
        Illumina=MapLocationsFunTruseqUP14J(MappingOutputFolders{puck},Bead,IlluminaBarcodesDownsampled,PrimerNLigationSequence,PrimerUPLigationSequence,NumBases,SaveData,PuckNames{puck});
    end
    if BeadType=="ReversePhase"
        Illumina=MapLocationsFunTruseqUP(MappingOutputFolders{puck},Bead,IlluminaBarcodesDownsampled,PrimerNLigationSequence,PrimerUPLigationSequence,NumBases,SaveData,PuckNames{puck});    
    end
    NumNearestBarcodes=[Illumina.NumNearestBarcodes];
    hammingdistances=[Illumina.HammingDistance];
    MappedLocations=[Illumina.MappedLocation];




    %We now want to invert the Illumina => SOLiD mapping. We want to find,
    %for each SOLiD barcode, how many Illumina barcodes map onto it at the minimum
    %Hamming distance.
    %I think we could just do this using the Unique function. 
    
    %we want to do this only over the set of illumina barcodes that have
    %one SOLiD partner
    
    %% Find mapping from SOLiD to Illumina mapping  
    
    %Note that everything in the Illumina structure is referenced to the
    %*unique* bead barcodes, so we have to do:
    [UniqueBeadBarcodes,BBFirstRef,BBOccCounts]=unique([Bead.Barcodes]);
    Bead=Bead(BBFirstRef);
    %This typically cuts out almost no beads
    
    %In the MapLocations function, whenever an illumina barcode has more
    %than one SOLiD partner, we only store the index of one of the
    %partners. This is not a problem here because we are only using
    %illumina barcodes with a unique partner.
    
    NumNearestBarcodes=[Illumina.NumNearestBarcodes];
    IlluminaUnique=Illumina(hammingdistances<=2 & NumNearestBarcodes==1);
    IlluminaBarcodesUnique=IlluminaBarcodesDownsampled(hammingdistances<=2 & NumNearestBarcodes==1);
    DGEUnique=DGESorted(:,hammingdistances<=2 & NumNearestBarcodes==1); %this is the same as DGEGoodBeads
    %The illumina barcodes that go in here have all been matched to a unique SOLiD
    %barcode at HD<=2. So each illumina barcode must show up in exactly one
    %bead. sum([MappedBeads.NumIlluminaBarcodes]) can only be less than
    %length(IlluminaBarcodesUnique) if barcodes dropout because a different barcode
    %is closer in edit space to the same bead.
    droppedbarcodes=cell(1,length(Bead));
    IndexOfMatchingIlluminaBarcodes=cell(1,length(Bead));
    MatchingIlluminaHDs=cell(1,length(Bead));
    NumMatchingIlluminaBarcodes=cell(1,length(Bead));
    parfor jl=1:length(Bead)
%        if floor(jl/1000)==jl/1000
%            disp(num2str(jl))
%        end
        matchingbarcodes=find([IlluminaUnique.IndexofBeadBarcode]==jl);
        if length(matchingbarcodes)==0
            NumMatchingIlluminaBarcodes{jl}=0;
            MatchingIlluminaHDs{jl}=0;            
            continue
        end
        minval=min([IlluminaUnique(matchingbarcodes).HammingDistance]);
        minindex=find([IlluminaUnique.HammingDistance]==minval & [IlluminaUnique.IndexofBeadBarcode]==jl);
        droppedbarcodes{jl}=length(matchingbarcodes)-length(minindex);
        NumMatchingIlluminaBarcodes{jl}=length(minindex);
        MatchingIlluminaHDs{jl}=minval;
        IndexOfMatchingIlluminaBarcodes{jl}=minindex;
    end
    figure(36)
    histogram(cell2mat(droppedbarcodes))
    title('For each SOLiD barcode, the number of Illumina Barcodes at greater than min HD')
    set(gca,'yscale','log');
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    
    droppedbarcodes=sum(cell2mat(droppedbarcodes));
    NumMatchingIlluminaBarcodes=cell2mat(NumMatchingIlluminaBarcodes);
    MatchingIlluminaHDs=cell2mat(MatchingIlluminaHDs);
    %NOTE THAT THE INDEX HERE IS NOT THE INDEX IN THE DGE -- it's the index
    %in the DGE reordered by reads.
    MappedBeads=struct('Barcodes',{Bead.Barcodes},'Locations',{Bead.Locations},'Pixels',{Bead.Pixels},'HammingDistances',num2cell(MatchingIlluminaHDs),'NumIlluminaBarcodes',num2cell(NumMatchingIlluminaBarcodes),'IndexofIlluminaBarcode',IndexOfMatchingIlluminaBarcodes);

    
    %% Find bijective mapping from SOLiD to Illumina
    
    %In an ideal world, we would actually check how many reads each
    %MappedBead has. So if a given bead has two nearest neighbors, we can
    %check to make sure it's not a sequencing error.
    UniqueMappedBeads=MappedBeads([MappedBeads.NumIlluminaBarcodes]==1);
%    [UniqueMappedBeads(:).IlluminaBarcode]=deal(IlluminaBarcodesUnique([UniqueMappedBeads.IndexofIlluminaBarcode]));

    UniqueMappedDGE=DGEUnique(:,[UniqueMappedBeads.IndexofIlluminaBarcode]);
    UniqueMappedIlluminaBarcodes=IlluminaBarcodesUnique([UniqueMappedBeads.IndexofIlluminaBarcode]);
    %We are somehow losing 1000 illumina barcodes on this step, which is
    %weird because there are only 67 cases of a solid barcode mapping to
    %two illumina barcodes.
    if SaveData
        save([MappingOutputFolders{puck},'MappedBeads.mat'],'MappedBeads','-v7.3')
        save([MappingOutputFolders{puck},'BijectiveMapping.mat'],'UniqueMappedBeads','UniqueMappedDGE','UniqueMappedIlluminaBarcodes','GeneNames','-v7.3');
    end
    %% Plot figures

    DGEsums=sum(UniqueMappedDGE,1);
    figure(10)
    histogram(DGEsums,1:10:1000)
    title('Total Transcripts per Barcode for Bijectively Mapped Barcodes')
    set(gca,'yscale','log');
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    
    figure(7)
    clf
    for k=1:length(DGEsums)
        rectangle('Position',[UniqueMappedBeads(k).Locations(1),UniqueMappedBeads(k).Locations(2),5*DGEsums(k)/mean(DGEsums),5*DGEsums(k)/mean(DGEsums)],...
          'Curvature',[1,1], 'FaceColor','r')
    end
    title('Reads per bead for bijective barcodes')
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    
    figure(13)
    clf
    hold on
    for k=1:3    
        subplot(1,3,k)
        HD2=find([Illumina.HammingDistance]==k-1);
        NumNearestBarcodes=[Illumina(HD2).NumNearestBarcodes];
        histogram(NumNearestBarcodes,1:1:5)
        set(gca,'yscale','log');
        if k==2
            title(['Num SOLiD barcodes matching each Illumina barcode - for HD 0, 1, 2.',num2str(k-1)]);
        end
    end
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    
    
    
    
    figure(19)
    histogram([MappedBeads.NumIlluminaBarcodes],0:1:10);
    title('Number of Illumina Barcodes per SOLiD barcode')
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');

    
    pixelsperbead=cellfun(@(x) length(x),{MappedBeads.Pixels});
    pixelsperbeadbij=cellfun(@(x) length(x),{UniqueMappedBeads.Pixels});    
    figure(15)
    clf
    histogram(pixelsperbead)
    title('Pixels per Unique SOLiD barcode')
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');    
    
    figure(16)
    clf
    histogram(pixelsperbeadbij)
    title('Pixels per Bijective Pair')
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');    
    
    figure(20)
    histogram([UniqueMappedBeads.HammingDistances],0:1:NumBases);
    title('Distribution of HD between Bead and Illumina Barcodes for Bijective Pairs')
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    
    %plot the locations of SOLiD barcodes with matched illumina barcodes
    figure(22)
    clf
    hold on
    goodbeads=find([MappedBeads.NumIlluminaBarcodes]==1);
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
    info=imfinfo([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif']);
    ROI=[[1,info.Height];[1,info.Width]];
    ROIHeight=ROI(1,2)-ROI(1,1)+1;
    ROIWidth=ROI(2,2)-ROI(2,1)+1;
    BeadImageMatchedBeads=false(ROIHeight,ROIWidth);
    %export_fig([OutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    for qr=1:length(goodbeads)
        BeadImageMatchedBeads(MappedBeads(goodbeads(qr)).Pixels)=true;
    end
    title('Beads mapped bijectively (Green) or not mapped (Red) to Illumina Barcodes');

    badbeads=find([MappedBeads.NumIlluminaBarcodes]==0);
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
    info=imfinfo([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif']);
    ROI=[[1,info.Height];[1,info.Width]];
    ROIHeight=ROI(1,2)-ROI(1,1)+1;
    ROIWidth=ROI(2,2)-ROI(2,1)+1;
    BeadImageBadBeads=false(ROIHeight,ROIWidth);
    %export_fig([OutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    for qr=1:length(badbeads)
        BeadImageBadBeads(MappedBeads(badbeads(qr)).Pixels)=true;
    end
    BeadImageFinal=zeros(ROIHeight,ROIWidth,3);
    BeadImageFinal(:,:,2)=BeadImageMatchedBeads;
    BeadImageFinal(:,:,1)=BeadImageBadBeads;    
    imshow(BeadImageFinal)
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');

    
    %Plot the locations of SOLiD barcodes lacking matched illumina barcodes
    figure(23)
    clf
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
    info=imfinfo([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif']);
    ROI=[[1,info.Height];[1,info.Width]];
    ROIHeight=ROI(1,2)-ROI(1,1)+1;
    ROIWidth=ROI(2,2)-ROI(2,1)+1;
    BeadImageRandomBeads=false(ROIHeight,ROIWidth);
    %export_fig([OutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    for qr=1:min(2000,length(goodbeads))
        BeadImageRandomBeads(MappedBeads(goodbeads(qr)).Pixels)=true;
    end
    title('2000 Random Beads matched to Illumina Barcodes')
    imshow(BeadImageRandomBeads)

    
    ReadsPerBead=sum(UniqueMappedDGE,1);
    [ReadsSorted,SortingIndices]=sort(ReadsPerBead,'descend');
    SortedUniqueMappedBeads=UniqueMappedBeads(SortingIndices);

    figure(24)
    clf
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
    info=imfinfo([BaseName,pad(num2str(1),2,'left','0'),' channel ',int2str(1),suffix,' transform.tif']);
    ROI=[[1,info.Height];[1,info.Width]];
    ROIHeight=ROI(1,2)-ROI(1,1)+1;
    ROIWidth=ROI(2,2)-ROI(2,1)+1;
    BeadImageTopBeads=false(ROIHeight,ROIWidth);
    %export_fig([OutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    for qr=1:min(2000,length(SortedUniqueMappedBeads))
        BeadImageTopBeads(SortedUniqueMappedBeads(qr).Pixels)=true;
    end
    title('Top 2000 Beads matched to Illumina Barcodes versus 2000 random unmatched beads')
    BeadImageRandomBadBeads=false(ROIHeight,ROIWidth);
    %export_fig([OutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    for qr=1:min(2000,length(badbeads))
        BeadImageRandomBadBeads(MappedBeads(badbeads(qr)).Pixels)=true;
    end
    BeadImageFinal2=zeros(ROIHeight,ROIWidth,3);
    BeadImageFinal2(:,:,2)=BeadImageTopBeads;
    BeadImageFinal2(:,:,1)=BeadImageRandomBadBeads;    
    imshow(BeadImageFinal2)
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    
    %we now just want to plot
    
    if 1
    DentateMarkers=find(GeneNames=='Prox1'|GeneNames=='Npnt'|GeneNames=='C1ql2');
    DentateMarkerSum=sum(UniqueMappedDGE(DentateMarkers,:),1);
    CAMarkers=find(GeneNames=='Cck'|GeneNames=='Fibcd1'|GeneNames=='Pvrl3'|GeneNames=='Kcnq5');
    CAMarkerSum=sum(UniqueMappedDGE(CAMarkers,:),1);
    INMarkers=find(GeneNames=='Gad2'|GeneNames=='Gad1'|GeneNames=='Vip'|GeneNames=='Npy'|GeneNames=='Sst'|GeneNames=='Pvalb'|GeneNames=='Lhx6');
    INMarkerSum=sum(UniqueMappedDGE(INMarkers,:),1);
    figure(25)
    clf
    BeadImageDentate=false(ROIHeight,ROIWidth);
    BeadImageCA=false(ROIHeight,ROIWidth);
    BeadImageIN=false(ROIHeight,ROIWidth);
    for qr=1:length(UniqueMappedBeads)
        if DentateMarkerSum(qr)>=1
            BeadImageDentate(UniqueMappedBeads(qr).Pixels)=true;
        end
        if CAMarkerSum(qr)>=1
            BeadImageCA(UniqueMappedBeads(qr).Pixels)=true;
        end
        if INMarkerSum(qr)>=1
            BeadImageIN(UniqueMappedBeads(qr).Pixels)=true;
        end
    end
    imshow(BeadImageDentate)
    title('Red: Dentate; Green: CA fields; Blue: Interneuron')
    BeadImageFinal3=zeros(ROIHeight,ROIWidth,3);
    BeadImageFinal3(:,:,3)=BeadImageIN;
    BeadImageFinal3(:,:,2)=BeadImageCA;
    BeadImageFinal3(:,:,1)=BeadImageDentate;    
    imshow(BeadImageFinal3)
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');
    end
    

    %% Plots of the base balance and errors in the Illumina-SOLiD mapping.
    %ALL PLOTS DISPLAYED THIS WAY ARE IN LIGATION ORDER, NOT SEQUENCE ORDER
    disp('NOTE: The following analysis depends on the number of ligations, the ligation sequence, and the bead structure, and may require modification.')
    %We want to make a plot of the color space balance for the final
    %mapped barcodes
    BaseBalanceBarcodes=[UniqueMappedBeads.Barcodes];
    %The base 5 representations of the basecalls are:
    BaseBalanceBase5Barcodes=cellfun(@(x) reverse(string(x)),{dec2base(BaseBalanceBarcodes,5,NumBases)},'UniformOutput',false);
    BaseBalanceBase5Barcodes=BaseBalanceBase5Barcodes{1};

    BaseBalanceMatrix=zeros(5,NumBases);
    for jp=1:NumBases
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
    figure(26)
    b=bar(BaseBalanceMatrix');
    b(1).FaceColor='k';
    b(2).FaceColor='b';
    b(3).FaceColor='g';
    b(4).FaceColor='y';
    b(5).FaceColor='r';
    title('Base balance per ligation for bijectively mapped barcodes - SOLiD');
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');


    colorspacesequence={};
    for seqnum=1:length(UniqueMappedIlluminaBarcodes)
        seq=char(UniqueMappedIlluminaBarcodes(seqnum));
        if BeadType=="ReversePhase"
            [PrimerNcolorspace,badflagN]=bs2cs(seq(1:6),PrimerNLigationSequence,'T','T');
            [PrimerUPcolorspace,badflagUP]=bs2cs(seq(7:13),PrimerUPLigationSequence,'A','');
        end
        if BeadType=="180402"
            [PrimerNcolorspace,badflagN]=bs2cs(seq(1:7),PrimerNLigationSequence,'T','T');
            [PrimerUPcolorspace,badflagUP]=bs2cs(seq(8:14),PrimerUPLigationSequence,'A','');            
        end
        tmptest=num2str(cat(2,PrimerNcolorspace,PrimerUPcolorspace));
        colorspacesequence{seqnum}=tmptest(~isspace(tmptest));
    end
    colorspacesequence=string(colorspacesequence);
    %The base 5 representations of the basecalls are:
%    BaseBalanceBase5Barcodes=cellfun(@(x) reverse(string(x)),{dec2base(BaseBalanceBarcodes,5,NumBases)},'UniformOutput',false);
%    BaseBalanceBase5Barcodes=BaseBalanceBase5Barcodes{1};

    BaseBalanceMatrix2=zeros(5,NumBases);
    for jp=1:NumBases
        testcmp0(jp)='0';
        testcmp1(jp)='1';
        testcmp2(jp)='2';
        testcmp3(jp)='3';
        testcmp4(jp)='4';
    end
    BaseBalanceMatrix2(1,:)=sum(char(colorspacesequence')==testcmp0,1);
    BaseBalanceMatrix2(2,:)=sum(char(colorspacesequence')==testcmp1,1);
    BaseBalanceMatrix2(3,:)=sum(char(colorspacesequence')==testcmp2,1);
    BaseBalanceMatrix2(4,:)=sum(char(colorspacesequence')==testcmp3,1);
    BaseBalanceMatrix2(5,:)=sum(char(colorspacesequence')==testcmp4,1);
    figure(27)
    b=bar(BaseBalanceMatrix2');
    b(1).FaceColor='k';
    b(2).FaceColor='b';
    b(3).FaceColor='g';
    b(4).FaceColor='y';
    b(5).FaceColor='r';
    
    title('Base balance per ligation for bijective barcodes expected from Illumina');
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');    
    
    figure(28)
    b=bar(BaseBalanceMatrix2'-BaseBalanceMatrix');
    b(1).FaceColor='k';
    b(2).FaceColor='b';
    b(3).FaceColor='g';
    b(4).FaceColor='y';
    b(5).FaceColor='r';
    
    title('Base balance of Illumina minus base balance of SOLiD in color space');
    %This plot shows BIAS: it does not show number of errors. A 0 in this
    %plot indicates that there is no bias, not that there are no errors.
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');    
    
    %this plot shows the actual errors. If there is disagreement between
    %Illumina and SOLiD, we plot the SOLiD barcode.
    ErrorMatrix=zeros(5,NumBases);
    ErrorBases=char(colorspacesequence')~=char(BaseBalanceBase5Barcodes);
    ErrorMatrix(1,:)=sum(char(colorspacesequence')==testcmp0 & ErrorBases,1);
    ErrorMatrix(2,:)=sum(char(colorspacesequence')==testcmp1 & ErrorBases,1);
    ErrorMatrix(3,:)=sum(char(colorspacesequence')==testcmp2 & ErrorBases,1);
    ErrorMatrix(4,:)=sum(char(colorspacesequence')==testcmp3 & ErrorBases,1);
    ErrorMatrix(5,:)=sum(char(colorspacesequence')==testcmp4 & ErrorBases,1);
    figure(29)
    b=bar(ErrorMatrix');
    b(1).FaceColor='k';
    b(2).FaceColor='b';
    b(3).FaceColor='g';
    b(4).FaceColor='y';
    b(5).FaceColor='r';

    title('Solid Color Expected on Mismatched Ligations');
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');    

    figure(291)
    b=bar(sum(ErrorBases,1));
    title('Number of mismatched ligations, per base')
    export_fig([MappingOutputFolders{puck},'Report_',PuckNames{puck},'.pdf'],'-append');    

    %% Output barcodes to CSV in for quality checking afterwards
    %NOTE: All barcodes outputted this way are in SEQUENCE ORDER, NOT IN
    %ORDER BY LIGATION

    %Bijective Pairs: SOLiD side
    UniqueBeadBarcodesForExport=char(replace(BaseBalanceBase5Barcodes,{'0','1','2','3','4'},{'N','B','G','O','R'}));
    if NumBases<14 %This is to deal with the InverseLigationSequence -- the export barcodes have to be 14 bases long
        UniqueBeadBarcodesForExport(:,NumBases+1:14)='N';
        UniqueBeadBarcodesForExport=UniqueBeadBarcodesForExport(:,WhichLigationsAreMissing);
    end
    UniqueBeadBarcodesForExport=UniqueBeadBarcodesForExport(:,InverseLigationSequence);
    %Bijective pairs: Illumina side
    UniqueIlluminaBarcodesForExport=char(replace(colorspacesequence',{'0','1','2','3','4'},{'N','B','G','O','R'}));
    if NumBases<14 %This is to deal with the InverseLigationSequence -- the export barcodes have to be 14 bases long
        UniqueIlluminaBarcodesForExport(:,NumBases+1:14)='N';
        UniqueIlluminaBarcodesForExport=UniqueIlluminaBarcodesForExport(:,WhichLigationsAreMissing);
    end
    UniqueIlluminaBarcodesForExport=UniqueIlluminaBarcodesForExport(:,InverseLigationSequence);
	    
    %All SOLiD beads passing filter
    AllBaseBalanceBarcodes=[Bead.Barcodes];
    %The base 5 representations of the basecalls are:
    BeadBarcodesForExport=cellfun(@(x) reverse(string(x)),{dec2base(AllBaseBalanceBarcodes,5,NumBases)},'UniformOutput',false);
    BeadBarcodesForExport=BeadBarcodesForExport{1};
    BeadBarcodesForExport=char(replace(BeadBarcodesForExport,{'0','1','2','3','4'},{'N','B','G','O','R'}));
    if NumBases<14 %This is to deal with the InverseLigationSequence -- the export barcodes have to be 14 bases long
        BeadBarcodesForExport(:,NumBases+1:14)='N';
        BeadBarcodesForExport=BeadBarcodesForExport(:,WhichLigationsAreMissing);
    end

    BeadBarcodesForExport=BeadBarcodesForExport(:,InverseLigationSequence);
    
    %All Illumina Barcodes >>> these are in order of the DGE
    Nseq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';
    allilluminacolorspacesequence={};
    for seqnum=1:length(IlluminaBarcodesDownsampled)
        seq=char(IlluminaBarcodesDownsampled(seqnum));
        if length(seq)<13
            seq=Nseq(1:13);
        end
        if BeadType=="ReversePhase"
            [PrimerNcolorspace,badflagN]=bs2cs(seq(1:6),PrimerNLigationSequence,'T','T');
            [PrimerUPcolorspace,badflagUP]=bs2cs(seq(7:13),PrimerUPLigationSequence,'A','');
        end
        if BeadType=="180402"
            [PrimerNcolorspace,badflagN]=bs2cs(seq(1:7),PrimerNLigationSequence,'T','T');
            [PrimerUPcolorspace,badflagUP]=bs2cs(seq(8:14),PrimerUPLigationSequence,'A','');            
        end
        tmptest=num2str(cat(2,PrimerNcolorspace,PrimerUPcolorspace));
        allilluminacolorspacesequence{seqnum}=tmptest(~isspace(tmptest));
    end
%    isemptyfn=cell2mat(cellfun(@(x) ~isempty(x),allilluminacolorspacesequence,'UniformOutput',false));
    allilluminacolorspacesequence=string(allilluminacolorspacesequence);
    IlluminaBarcodesForExport=char(replace(allilluminacolorspacesequence',{'0','1','2','3','4'},{'N','B','G','O','R'}));
    if NumBases<14 %This is to deal with the InverseLigationSequence -- the export barcodes have to be 14 bases long
        IlluminaBarcodesForExport(:,NumBases+1:14)='N';
        IlluminaBarcodesForExport=IlluminaBarcodesForExport(:,WhichLigationsAreMissing);
    end
    IlluminaBarcodesForExport=IlluminaBarcodesForExport(:,InverseLigationSequence);
    
    save([MappingOutputFolders{puck},'ReadableSequences.mat'],'UniqueBeadBarcodesForExport','UniqueIlluminaBarcodesForExport','BeadBarcodesForExport','IlluminaBarcodesForExport','-v7.3')
    csvwrite([MappingOutputFolders{puck},'BijectiveBeadBarcodes.csv'],UniqueBeadBarcodesForExport);
    csvwrite([MappingOutputFolders{puck},'BijectiveIlluminaBarcodes.csv'],UniqueIlluminaBarcodesForExport);
    csvwrite([MappingOutputFolders{puck},'AllBeadBarcodes.csv'],BeadBarcodesForExport);
    csvwrite([MappingOutputFolders{puck},'AnalyzedIlluminaBarcodes.csv'],IlluminaBarcodesForExport);
    
    %% Print output
    fileid=fopen([MappingOutputFolders{puck},'Metrics.txt'],'a');
    fprintf(fileid,['\n\nThere were ',num2str(length(Bead)),' unique bead barcodes.\n',...
    'We analyzed ',num2str(IlluminaBarcodesToAnalyze),' Illumina barcodes with more than ',num2str(IlluminaReadThreshold),' reads.\n',...
    'There were ',num2str(length(IlluminaBarcodesUnique)),' Illumina barcodes uniquely matching a SOLiD barcode at HD<=2.\n',...
    'There were ',num2str(length(UniqueMappedIlluminaBarcodes)),' bijective pairings between SOLiD and Illumina barcodes.\n',...
    'A total of ',num2str(sum([MappedBeads([MappedBeads.NumIlluminaBarcodes]>1).NumIlluminaBarcodes])),' Illumina barcodes were dropped because multiple Illumina barcodes mapped onto a single SOLiD barcode.\n',...
    'A total of ',num2str(droppedbarcodes),' barcodes were dropped because the bead they mapped onto had another Illumina barcode that mapped onto it with lower Hamming distance.'
    ]);
    fclose(fileid);
    
    
    
    %
    %SOLiD barcodes, how many Illumina barcodes have been mapped onto it?
    %The Unique function already calculates this.
    

    %We can now add data to Bead. First, for each Bead, we want to take its
    %location, pixels, barcode, and its location in the DGE
%    UniqueBeads=

    
end
for puck=1:length(PuckNames)  
    mkdir(['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
    copyfile(MappingOutputFolders{puck},['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
end

%% Evaluate

%We can now use the Hough transform, and ask what percentage of Hough beads
%also have barcodes, as a way of analyzing what percentage of the beads
%were called. We call the Hough transform and find the centroids, and then
%ask about the fraction of Hough beads with centroids within a radius of a sequenced
%barcode ;

%We want to produce a plot with the fraction of illumina barcodes that are direct
%matches with a surface bead; one base away from a surface bead; two bases,
%etc., and also for the negative control