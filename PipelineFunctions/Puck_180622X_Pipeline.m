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

BalanceZScores=1; %This determines whether the basecalling should be balanced. See BeadSeqFun
SaveData=1;
IlluminaReadThreshold=10;
MaximumBarcodesToAnalyze=80000;
NumLigations=20; %We assume that the missing ligation is Truseq-4 if you are doing 13.
NumBases=14; %This is the number of bases sequenced
BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14];
%BarcodeSequence=[1,2,3,4,0,5,0,6,7,8,9,10,0,11,0,12,0,13]; %this determines how the numerical barcode is built from the ligations. Basically, the constant bases should be 0, and all other bases should be consecutive
FolderWithRawTiffs='D:\Slideseq\Raw\180626 - Pucks 180622-X\';
FolderWithProcessedTiffs='D:\Slideseq\Processed\';
tmpfolder='C:\Users\sgr\pucktmp\';
IndexFiles={'primer t through up-2','primer up-3 through up-4'}; %give the prefixes of each of the files. The next character after should be 't',
LigationToIndexFileMapping=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2];%for ligations 1:20, which file number are the ligations found in?
%for ligations 1:20, which value of t are they, within their file? This
%could also be deduced from the LigationToIndexFileMapping array
tnumMapping=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1,2,3,4];
PucksToAnalyze=[1,2,3,4,5,6,7,8,9,10,11,12,15,16,17,18,19,20,21,22,23,24,25,26]; %These are the values of xy from the .nd2 file that will be analyzed.
%Note that the order in PuckNames should match the order in the .nd2 file.
PuckNames={'Puck_180622_1','Puck_180622_2','Puck_180622_3','Puck_180622_4','Puck_180622_5','Puck_180622_6','Puck_180622_7','Puck_180622_8','Puck_180622_9','Puck_180622_10','Puck_180622_11','Puck_180622_12','Puck_180622_13','Puck_180622_14','Puck_180622_15','Puck_180622_16','Puck_180622_17','Puck_180622_18','Puck_180622_19','Puck_180622_20','Puck_180622_21','Puck_180622_22','Puck_180622_23','Puck_180622_24'}; %give the names of the pucks
OutputFolderRoot='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\';%\Puck_180106_1\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_2\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_3\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_4\'};
IlluminaFolderRoot='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\';
%SequencingFileName={'fei_S3_L001_R1_001.fastq'};
ImageSize=6030; %The image registration will output images that are not all the exact same size, because of stitching. So find_roi_stack_fun crops the images a bit. We can choose this value to be something like 0.95 * the size of the images. So e.g. for 3x3 it should be 0.95*2048*(2*(2/3)+1) = 4500. For 7x7 it is 0.95*2048*(1+6*2/3)=9728. I used 10400 previously though.
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
    OutputFolders{puck}=[OutputFolderRoot,PuckNames{puck},'\'];
    mkdir(OutputFolders{puck});    
end

display('Renaming Files');

%% Move files
for pucknum=1:length(PuckNames)
    mkdir([tmpfolder,PuckNames{pucknum}]);
    puck=PucksToAnalyze(pucknum);
    for ligation=1:NumLigations
        tnum=tnumMapping(ligation);
        if exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),'xy',num2str(puck),'.tif'],'file')
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),'xy',num2str(puck),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),'xy',pad(num2str(puck),2,'left','0'),'.tif'],'file')
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),'xy',pad(num2str(puck),2,'left','0'),'.tif'];        
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),'xy',num2str(puck),'.tif'],'file')
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),'xy',num2str(puck),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),'xy',pad(num2str(puck),2,'left','0'),'.tif'],'file')
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),'xy',pad(num2str(puck),2,'left','0'),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',num2str(tnum),'.tif'],'file')
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',num2str(tnum),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',pad(num2str(tnum),2,'left','0'),'.tif'],'file')
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'t',pad(num2str(tnum),2,'left','0'),'.tif'];
        elseif exist([FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'.tif'],'file')
            filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'xy',num2str(puck),'.tif'];
        else
            assert(1==0)
        end
            
            outputfilename=[ProcessedImageFolders{pucknum},PuckNames{pucknum},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'];
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

for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
%We now send the command to do the registration with find_roi
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
    disp(['Beginning basecalling for puck number ',num2str(puck)])
	[Bead BeadImage]=BeadSeqFun2(BaseName,suffix,OutputFolders{puck},BeadZeroThreshold,BarcodeSequence,20,NumLigations,PuckNames{puck},BalanceZScores);
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end

%% Match Illumina barcodes - to rerun the bead mapping only, start here

display('Matching Illumina Barcodes')
for puck=1:length(PuckNames)
    MapIlluminaToPuck(PuckNames{puck},BeadType,'PrimerNLigationSequence',PrimerNLigationSequence,'PrimerUPLigationSequence',PrimerUPLigationSequence,'InverseLigationSequence',InverseLigationSequence,'WhichLigationsAreMissing',WhichLigationsAreMissing,'IlluminaReadThreshold',IlluminaReadThreshold,'MaximumBarcodesToAnalyze',MaximumBarcodesToAnalyze,'OutputFolder',OutputFolders{puck},'IlluminaRootFolder',IlluminaFolderRoot,'NumBases',NumBases,'NumLigations',NumLigations,'BarcodeSequence',BarcodeSequence,'ProcessedImageFolder',ProcessedImageFolders{puck})    
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