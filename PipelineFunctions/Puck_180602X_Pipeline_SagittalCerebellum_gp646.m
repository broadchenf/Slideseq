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


%For future releases, consider setting the priority of the parpool
%processes to High:

%matlabpool open
%cmd_str = 'wmic process where name="MATLAB.exe" CALL setpriority 64';
%[~,~] = system(cmd_str);


%% Initialize
%Set up the paths to all of the various folders with functions in them that
%we will call:
clear all
close all

PythonPath='C:\Users\sgr\AppData\Local\Programs\Python\Python37\python.exe';%for gpdc2
DropboxPath='C:\Users\sgr\Dropbox\Projects\Project - SlideSeq\';
addpath('C:\Fiji.app\scripts',[DropboxPath,'BeadSeq Code'],[DropboxPath,'BeadSeq Code\find_roi'],[DropboxPath,'Pucks\PipelineFunctions']);
javaaddpath('C:\Program Files\MATLAB\R2017a\java\mij.jar');
javaaddpath(which('MatlabGarbageCollector.jar'))
%We assume that the nd2s have been exported to tiffs in the format:
%DescriptiveNametLYX, where DescriptiveName refers to one run of the microscope,  Y is a letter and X is a number, and L is the
%name of the ligation within that DescriptiveName file series.
%We convert the files to a final output format:
%Puck85 Ligation X Position AB
%BeadType="180402"; for 14bp barcodes from 180402 beads
BeadType="ReversePhase";

RenameFiles=1;

RunSignificanceAnalysis=0;
NumClusters=[11,17]; %Cerebellum is 1, hippocampus is 2

%Which pucks do we try to run analogizer on? 0 if don't run, otherwise, the
%number here refers to the element of DropseqDGEPaths and
%DropseqClusterPaths to use as a reference.
RunAnalogizer=2*ones(1,30);
AnalogizerBeadCutoff=10;
AnalogizerType="NMFReg";
%Cerebellum is 1. Hippocampus is 2.
DropseqDGEPaths={'\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT','\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Hippocampus'};
DropseqClusterPaths={'\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT\assign\F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS','\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Hippocampus\assign\F_GRCm38.81.P60Hippocampus.cluster.assign.RDS'};
DropseqMeanAndVariancePath=fullfile(DropboxPath,'BeadSeq Code\DGEMeansAndVariances');

NumPar=12;

CropImage=1;

EnforceBaseBalance=1;
BaseBalanceTolerance=0.05;

SaveData=1;
IlluminaReadThreshold=10;
MaximumBarcodesToAnalyze=160000;
NumLigations=20; %We assume that the missing ligation is Truseq-4 if you are doing 13.
NumBases=14; %This is the number of bases sequenced
BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14];
%BarcodeSequence=[1,2,3,4,0,5,0,6,7,8,9,10,0,11,0,12,0,13]; %this determines how the numerical barcode is built from the ligations. Basically, the constant bases should be 0, and all other bases should be consecutive
FolderWithRawTiffs='D:\Slideseq\Raw\180816 - Pucks 180728-X\';
FolderWithProcessedTiffs='D:\Slideseq\Processed\';
tmpfolder='D:\pucktmp\';
IndexFiles={'primers t through t-4','primers t-2 through t-4','primers up through up-4 then t-2 through t-4'}; %give the prefixes of each of the files. The next character after should be 't',
LigationToIndexFileMapping=[1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3];%for ligations 1:20, which file number are the ligations found in?
%for ligations 1:20, which value of t are they, within their file? This
%could also be deduced from the LigationToIndexFileMapping array
tnumMapping=[1,2,3,4,1,2,3,4,5,6,1,2,3,4,5,6,7,8,9,10];
PucksToAnalyze=15:24; %These are the values of xy from the .nd2 file that will be analyzed.
%Note that the order in PuckNames should match the order in the .nd2 file.
PuckNames={'Puck_180602_1','Puck_180602_2','Puck_180602_3','Puck_180602_4','Puck_180602_5','Puck_180602_6','Puck_180602_7','Puck_180602_8','Puck_180602_9','Puck_180602_10','Puck_180602_11','Puck_180602_12','Puck_180602_13','Puck_180602_14','Puck_180602_15','Puck_180602_16','Puck_180602_17','Puck_180602_18','Puck_180602_19','Puck_180602_20','Puck_180602_21','Puck_180602_22','Puck_180602_23','Puck_180602_24'}; %give the names of the pucks
OutputFolderRoot=[DropboxPath,'\Pucks\Barcodes\'];%\Puck_180106_1\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_2\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_3\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_4\'
IlluminaFolderRoot='C:\Slideseq\Illumina\';

%SequencingFileName={'fei_S3_L001_R1_001.fastq'};
ImageSize=6030; %The image registration will output images that are not all the exact same size, because of stitching. So find_roi_stack_fun crops the images a bit. We can choose this value to be something like 0.95 * the size of the images. So e.g. for 3x3 it should be 0.95*2048*(2*(2/3)+1) = 4500. For 7x7 it is 0.95*2048*(1+6*2/3)=9728. I used 10400 previously though.
XCorrBounds=[2800,3200,2800,3200]; %This is the ROI used in channel registration
RegisterColorChannels=1;
BeadZeroThreshold=1;
PixelCutoffRegistration=400;
PixelCutoffBasecalling=300;

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

if RenameFiles
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
end

%% Registration

display('Image Registration')

for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
    display(['Beginning registration on puck number ',num2str(puck)])
    %We now send the command to do the registration with find_roi
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    Suffix='_Stitched';
    find_roi_stack_fun(BaseName,Suffix,ImageSize,'PixelCutoff',PixelCutoffRegistration,'XCorrBounds',XCorrBounds,'RegisterColorChannels',1,'NumPar',NumPar);
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
	[Bead BeadImage]=BeadSeqFun6(BaseName,suffix,OutputFolders{puck},BeadZeroThreshold,BarcodeSequence,NumPar,NumLigations,PuckNames{puck},EnforceBaseBalance,BaseBalanceTolerance,'PixelCutoff',PixelCutoffBasecalling);
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end

%% Match Illumina barcodes - to rerun the bead mapping only, start here

display('Matching Illumina Barcodes')
%The illumina DGE should have the puck date and number, i.e. 180728_3,
%followed by an underscore, so 180728_3_
%The underscore is to distinguish 180728_2_ from 180728_23, for example
MappingOutputFolders={};
for puck=1:length(PuckNames)
    MappingOutputFolder=MapIlluminaToPuck(PuckNames{puck},BeadType,'PrimerNLigationSequence',PrimerNLigationSequence,'PrimerUPLigationSequence',PrimerUPLigationSequence,'InverseLigationSequence',InverseLigationSequence,'WhichLigationsAreMissing',WhichLigationsAreMissing,'IlluminaReadThreshold',IlluminaReadThreshold,'MaximumBarcodesToAnalyze',MaximumBarcodesToAnalyze,'OutputFolder',OutputFolders{puck},'IlluminaRootFolder',IlluminaFolderRoot,'NumBases',NumBases,'NumLigations',NumLigations,'BarcodeSequence',BarcodeSequence,'ProcessedImageFolder',ProcessedImageFolders{puck},'ImageSize',ImageSize,'NumPar',NumPar);    
    if (MappingOutputFolder=="NA")
        continue
    end
    MappingOutputFolders{puck}=MappingOutputFolder;
    %Record parameters into a file, "IlluminaParams.txt"
    fileid=fopen(fullfile(MappingOutputFolder,'IlluminaParams.txt'),'a');
    fprintf(fileid,['\n\nThe parameters for this Bead-To-Illumina matching run were:\n','\nBead Type: ',char(BeadType),'\nPrimerNLigationSequence: ',num2str(PrimerNLigationSequence),'\nPrimerUPLigationSequence: ',num2str(PrimerUPLigationSequence),'\nInverseLigationSequence: ',num2str(InverseLigationSequence),'\nWhichLigationsAreMissing: ',num2str(WhichLigationsAreMissing),'\nIlluminaReadThreshold: ',num2str(IlluminaReadThreshold),'\nMaximumBarcodesToAnalyze: ',num2str(MaximumBarcodesToAnalyze),'\nNumBases: ',num2str(NumBases),'\nNumLigations: ',num2str(NumLigations),'\nBarcodeSequence: ',num2str(BarcodeSequence)]);
    fclose(fileid);
    ConvertMatToR(OutputFolders{puck});
    tmp=strsplit(MappingOutputFolders{puck},'\');
    mappingstarttimereadable=tmp(end-1);
    mappingstarttimereadable=mappingstarttimereadable{1};
    mkdir(['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
    copyfile(MappingOutputFolders{puck},['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
end

if CropImage
    for puck=1:length(PuckNames)
        try
            BeadMappingFile=FindMostRecentMapping(OutputFolders{puck});
        catch
            continue
        end
        MappingOutputFolder=fullfile(OutputFolders{puck},BeadMappingFile);
        if ~exist(fullfile(MappingOutputFolder,'IlluminaParams.txt')) %Puck hasn't been mapped
            continue
        end
        load(fullfile(MappingOutputFolder,'BijectiveMapping.mat'))
        image=PlotGeneFromName('Malat1',GeneNames,UniqueMappedDGE,UniqueMappedBeads,'Overlay',1);
        h = imfreehand; %draw something 
        M = h.createMask();
        
        disp('Cropping complete')
        Locs=[UniqueMappedBeads.Locations];
%        GoodBeads=Locs(1,:) > croprectangle(1) & Locs(1,:) < (croprectangle(1)+croprectangle(3)) & Locs(2,:) > croprectangle(2) & Locs(2,:)<croprectangle(2)+croprectangle(4);
        GoodBeads=M(sub2ind(size(image),round(Locs(2,:)),round(Locs(1,:))))==1;
        UniqueMappedBeads=UniqueMappedBeads(GoodBeads);
        UniqueMappedDGE=UniqueMappedDGE(:,GoodBeads);
        UniqueMappedIlluminaBarcodes=UniqueMappedIlluminaBarcodes(GoodBeads);
        movefile(fullfile(MappingOutputFolder,'BeadLocationsForR.csv'),fullfile(MappingOutputFolder,'BeadLocationsForRUncropped.csv'));
        movefile(fullfile(MappingOutputFolder,'MappedDGEForR.csv'),fullfile(MappingOutputFolder,'MappedDGEForRUncropped.csv'));    
        movefile(fullfile(MappingOutputFolder,'BijectiveMapping.mat'),fullfile(MappingOutputFolder,'BijectiveMappingUncropped.mat'));    
        save(fullfile(MappingOutputFolder,'BijectiveMapping.mat'),'UniqueMappedBeads','UniqueMappedDGE','UniqueMappedIlluminaBarcodes','GeneNames','-v7.3');
        ConvertMatToR(OutputFolders{puck});    
    end
end

for puck=1:length(PuckNames)
    if ~RunAnalogizer(puck)
        continue
    end
    AtlasReference=RunAnalogizer(puck);
    DropseqDGEPath=DropseqDGEPaths{AtlasReference};
    DropseqClusterPath=DropseqClusterPaths{AtlasReference};
    display('Running Analogizer')
    %This will run the analogizer on the most recent Illumina Barcoding run.
    try
        BeadMappingFile=FindMostRecentMapping(OutputFolders{puck});
    catch
        continue
    end
    MappingOutputFolder=fullfile(OutputFolders{puck},BeadMappingFile);
    if ~exist(fullfile(MappingOutputFolder,'IlluminaParams.txt')) %Puck hasn't been mapped
        continue
    end
    if exist(fullfile(MappingOutputFolder,'AnalogizerParams.txt'))
        c=clock;
        analogizerstarttimereadable=[num2str(c(2)),'-',num2str(c(3)),'_',pad(num2str(c(4)),2,'left','0'),pad(num2str(c(5)),2,'left','0')];

        %If prior analogizer data exists, we save it and its parameter file
        %into a new labeled folder.
        mkdir(fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))
        try
            try
                movefile(fullfile(MappingOutputFolder,'Cluster*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))            
            catch
            end
            try
                movefile(fullfile(MappingOutputFolder,'Analogizer*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))            
            catch
            end
            %This is a failsafe
            try
                copyfile(fullfile(MappingOutputFolder,'Analogizer*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))
            catch
            end
            try
                copyfile(fullfile(MappingOutputFolder,'Cluster*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))
            catch
            end

            delete(fullfile(MappingOutputFolder,'Analogizer*'))
            delete(fullfile(MappingOutputFolder,'AnalogizerParams.txt'))
            delete(fullfile(MappingOutputFolder,'Cluster*'))
        catch

        end
    end
    thisbeadcutoff=AnalogizerBeadCutoff;
    
    while ~exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'))
        commandfile=fopen('C:\Analogizer.cmd','w');
        %fwrite(commandfile,['cd "C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code" & dir']);
        if AnalogizerType=="Analogizer"
            fwrite(commandfile,['cd "C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code" & "C:\Program Files\R\R-3.4.3\bin\Rscript.exe" AnalogizerScript.R ',PuckNames{puck},'\',BeadMappingFile,' ',DropseqDGEPath,' ',DropseqClusterPath,' ',num2str(thisbeadcutoff)]);
        elseif AnalogizerType=="NMFReg"
            %Bob says:
%           -da = data path atlas
%           -dp = data path puck
%           -t = tissue type (only cerebellum and hippo now
%           - c = UMI cutoff
%           - dge = DGE name minus file extension
%           - bl = Bead locations minus file extension
            switch RunAnalogizer(puck)
                case 1
                    tissuetype='cerebellum';
                case 2
                    tissuetype='hippocampus';
                otherwise
                    assert(1==0)
            end
            fwrite(commandfile,['C: & cd "',DropboxPath,'BeadSeq Code" & "',PythonPath,'" autoNMFreg_windows.py -da \\iodine-cifs\broad_macosko\data\NMFreg\data -dp "',fullfile(OutputFolders{puck},BeadMappingFile),'" -t ',tissuetype,' -c ',num2str(thisbeadcutoff),' -dge MappedDGEForR -bl BeadLocationsForR']);            
        end
        fclose(commandfile);
        !C:/Analogizer
        if ~exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'))
            disp('Analogizer Failed. Increasing bead cutoff and continuing.')
            thisbeadcutoff=thisbeadcutoff+5;
        end
        if AnalogizerType=="NMFReg" && exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv')) && ~exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignmentsOriginal.csv'))
            %NMFReg outputs the cluster assignments in the wrong format
            opts=detectImportOptions(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'));
            ClusterAssignments=readtable(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'),opts,'ReadVariableNames',true);
            newtable=table(ClusterAssignments.barcode,ClusterAssignments.atlas_cluster,'VariableNames',{'Var1','x'});
            movefile(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'),fullfile(MappingOutputFolder,'AnalogizerClusterAssignmentsOriginal.csv'))
            writetable(newtable,fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'));
        end
    end
    switch RunAnalogizer(puck)
        case 1
            copyfile(fullfile(DropboxPath,'BeadSeq Code\DGEMeansAndVariances\CerebellumVarianceByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerVarianceByDropseqCluster.csv'));
            copyfile(fullfile(DropboxPath,'BeadSeq Code\DGEMeansAndVariances\CerebellumExpressionByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerExpressionByDropseqCluster.csv'));
        case 2
            copyfile(fullfile(DropboxPath,'BeadSeq Code\DGEMeansAndVariances\HippocampusVarianceByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerVarianceByDropseqCluster.csv'));
            copyfile(fullfile(DropboxPath,'BeadSeq Code\DGEMeansAndVariances\HippocampusExpressionByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerExpressionByDropseqCluster.csv'));
    end
    %Record parameters into a file, "AnalogizerParams.txt"
    fileid=fopen([MappingOutputFolder,'\AnalogizerParams.txt'],'a');
    fprintf(fileid,['\n\nThe parameters for this Bead-To-Illumina matching run were:\n',...
    '\nAnalogizer Type: ',char(AnalogizerType),...
    '\nAnalogizer Bead Cutoff: ',num2str(thisbeadcutoff),...
    '\nDropseq DGE Path: ',DropseqDGEPath,...
    '\nDropseq Cluster Path: ',DropseqClusterPath,...
    ]);
    fclose(fileid);
    if RunSignificanceAnalysis
        for k=1:NumClusters(AtlasReference)
            PermutationTestByCluster(OutputFolders{puck},k)
        end
    end
end


for puck=1:length(PuckNames)  
    tmp=strsplit(MappingOutputFolders{puck},'\');
    mappingstarttimereadable=tmp(end-1);
    mappingstarttimereadable=mappingstarttimereadable{1};
    mkdir(['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
    copyfile(MappingOutputFolders{puck},['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
end

%% Evaluate

%We can now use the Hough transform, and ask what percentage of Hough beads
%also have barcodes, as a way of analyzing what percentage of the beads
%were called. We call the Hough transform and find the centroids, and then
%ask about the fraction of Hough beads with centroids within a radius of a sequenced
%barcode ;

%We want to produce a plot with the fraction of illumina barcodes that are diporect
%matches with a surface bead; one base away from a surface bead; two bases,
%etc., and also for the negative control