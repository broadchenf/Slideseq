%%%%SETUP:

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

%6) NOTE: SlideseqStitch.ijm is stored in C:\Fiji.app\macros, and the text
%is:
%params=getArgument()
%print("Running stitching with arguments:")
%print(params)
%run("Grid/Collection stitching", params);
%selectWindow('Fused');
%run("Save", "Tiff..., path=PuckOutputTmp.tif");
%print("Done.")
%eval("script","System.exit();");

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
IlluminaReadThreshold=5;
MaximumBarcodesToAnalyze=30000;
BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14]; %this determines how the numerical barcode is built from the ligations. Basically, the constant bases should be 0, and all other bases should be consecutive
FolderWithRawTiffs='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Raw\180123 - Pucks 180106-X\';
FolderWithProcessedTiffs='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Processed\';
tmpfolder='C:\Users\sgr\pucktmp\';
IndexFiles={'pucks180106_n_through_n-4_t','pucks180106_up_through_up-3_t','pucks180106_up-4_t'}; %give the prefixes of each of the files. The next character after should be 't',
LigationToIndexFileMapping=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3];%for ligations 1:20, which file number are the ligations found in?
%for ligations 1:20, which value of t are they, within their file? This
%could also be deduced from the LigationToIndexFileMapping array
tnumMapping=[1,2,3,4,5,6,7,8,9,10,1,2,5,6,7,8,9,10,1,2];
PuckNames={'Puck_180106_1','Puck_180106_2','Puck_180106_3','Puck_180106_4'}; %give the names of the pucks
OutputFolderRoot='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\';%\Puck_180106_1\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_2\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_3\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_180106_4\'};
IlluminaFolderRoot='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\';
%SequencingFileName={'fei_S3_L001_R1_001.fastq'};
ImageSize=4500; %The image registration will output images that are not all the exact same size, because of stitching. So find_roi_stack_fun crops the images a bit. We can choose this value to be something like 0.95 * the size of the images. So e.g. for 3x3 it should be 0.95*2048*(2*(2/3)+1) = 4500. For 7x7 it is 0.95*2048*(1+6*2/3)=9728. I used 10400 previously though.
TopBarcodeNum=30000; %when it comes to analyzing barcodes from the illumina data, we only consider the TopBarcodeNum barcodes with the most reads. This should roughly be the number of beads on the puck.
BeadZeroThreshold=1;

%TO DO: Propagate these variables and the number of bases through the code --
%currently 14 ligations is hardcoded in most places.
PrimerNLigationSequence = [2, 7, 1, 6, 5, 4, 3];
PrimerUPLigationSequence=[2, 7, 1, 6, 5, 4, 3];
InverseLigationSequence=[3,1,7,6,5,4,2,10,8,14,13,12,11,9];


OutputFolders={};
for puck=1:length(PuckNames)
    ProcessedImageFolders{puck}=[FolderWithProcessedTiffs,PuckNames{puck},'\'];
    mkdir([FolderWithProcessedTiffs,PuckNames{puck}]);
end
OutputXIndices=[[3,5];[9,11];[3,5];[9,11]]; %for each puck, give the X indices where it begins and ends
OutputYIndices=[[3,5];[2,4];[9,11];[8,10]];

display('Renaming Files');

%% Move files
for puck=1:length(PuckNames)
    mkdir([tmpfolder,PuckNames{puck}]);
    for ligation=1:20
        tnum=tnumMapping(ligation);
        for yposition=OutputYIndices(puck,1):OutputYIndices(puck,2)
            switch yposition
                case 1
                    yposstring='A';
                case 2
                    yposstring='B';
                case 3
                    yposstring='C';
                case 4
                    yposstring='D';
                case 5
                    yposstring='E';
                case 6
                    yposstring='F';
                case 7
                    yposstring='G';
                case 8
                    yposstring='H';
                case 9
                    yposstring='I';
                case 10
                    yposstring='J';
                case 11
                    yposstring='K';
                case 12
                    yposstring='L';
                case 13
                    yposstring='M';
            end
            for xposition=OutputXIndices(puck,1):OutputXIndices(puck,2)
                xpositionindex=find(OutputXIndices(puck,1):OutputXIndices(puck,2)==xposition);
                ypositionindex=find(OutputYIndices(puck,1):OutputYIndices(puck,2)==yposition);
                filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),yposstring,pad(num2str(xposition),2,'left','0'),'.tif'];
                if ~exist(filename)
                    filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),yposstring,pad(num2str(xposition),2,'left','0'),'.tif'];                    
                end
                outputfilename=[tmpfolder,PuckNames{puck},'\',PuckNames{puck},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Position_X_',pad(num2str(xpositionindex),2,'left','0'),'_Y_',pad(num2str(ypositionindex),2,'left','0'),'.tif'];
%                outputfilename=[OutputFolders{puck},PuckNames{puck},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Position_X_',pad(num2str(xpositionindex),2,'left','0'),'_Y_',pad(num2str(yposition),2,'left','0'),'.tif'];
                copyfile(filename,outputfilename)
            end
        end
    end
end
%We now open Miji instances and send the command, using parfor, to do the
%stitching
%% Stitching
display('Stitching Images')

for puck=1:length(PuckNames)
thispuckname=PuckNames{puck};
thisinputfolder=[tmpfolder,PuckNames{puck},'\'];
thisoutputfolder=ProcessedImageFolders{puck};

for ligation=1:20

    
    
    %THIS CODE WAS FOR USING MIJI. But there is some memory leak and Miji
    %slows down dramatically by about the 10th ligation.
%    Miji(false); %we don't want to start the FIJI gui
    %https://www.mathworks.com/matlabcentral/fileexchange/47545-mij--running-imagej-and-fiji-within-matlab
    %We got the following command by using the "record" function for macros in ImageJ.
    %We don't really need to do the replace replace replace here, we could
    %just concatenate strings. Unclear why we did it this way...
    command=replace(replace(replace(['type=[Filename defined position] order=[Defined by filename         ] grid_size_x=',num2str(OutputXIndices(puck,2)-OutputXIndices(puck,1)+1),' grid_size_y=',num2str(OutputYIndices(puck,2)-OutputYIndices(puck,1)+1),' tile_overlap=30 first_file_index_x=01 first_file_index_y=01 directory=DIRECTORYPATH file_names=PUCKNAME_Ligation_LL_Position_X_{xx}_Y_{yy}.tif output_textfile_name=Ligation_LL.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]'],'LL',pad(num2str(ligation),2,'left','0')),'PUCKNAME',thispuckname),'DIRECTORYPATH',thisinputfolder);
%    MIJ.run('Grid/Collection stitching',command);
    %We now need to make sure it saves the stitched image. We save to files
    %of the form 'Puck_85_Ligation_1_Stitched'
%    MIJ.selectWindow('Fused');
%    MIJ.run("Save",replace('Tiff..., path=[OUTPUTPATH]','OUTPUTPATH',replace([thisinputfolder,thispuckname,'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'],'\','\\')));
%    MIJ.run("Close");
%    MIJ.run("Quit");
    
%    movefile([thisinputfolder,thispuckname,'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'],[thisoutputfolder,thispuckname,'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'])
    commandfile=fopen('C:\FijiCommand.cmd','w');
    fwrite(commandfile,strcat('C:\Fiji.app\ImageJ-win64.exe --headless --console -macro SlideseqStitch.ijm "',command,'"'));
    fclose(commandfile);
    !C:/FijiCommand
    if exist('C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\PuckOutputTmp.tif')
        movefile('C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\PuckOutputTmp.tif',[thisoutputfolder,thispuckname,'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'])
    else
        movefile('C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\PipelineFunctions\PuckOutputTmp.tif',[thisoutputfolder,thispuckname,'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'])
    end
    %jheapcl
end
end

%Now, issue a command at the command line to change the smart sync status
%of the original input folder

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
%we now run Bead_Seq itself. Again, where parallelization is non-trivial in
%Bead_Seq, it is implemented naively using parfor.

display('Base Calling')

c=clock;
starttimereadable=[num2str(c(2)),'-',num2str(c(3)),'_',pad(num2str(c(4)),2,'left','0'),pad(num2str(c(5)),2,'left','0')];

for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
%We now send the command to do the registration with find_roi
    
    OutputFolders{puck}=[OutputFolderRoot,PuckNames{puck},'_',starttimereadable,'\'];
    mkdir(OutputFolders{puck});

    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
	[Bead BeadImage]=BeadSeqFun(BaseName,suffix,OutputFolders{puck},BeadZeroThreshold,BarcodeSequence,20);
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end

%% Match Illumina barcodes

display('Matching Illumina Barcodes')


for puck=1:length(PuckNames)
    %Load in the CSV:
    if ~(exist([IlluminaFolderRoot,PuckNames{puck},'\DGE.csv']))
        continue
    end
    DGETable=readtable([IlluminaFolderRoot,PuckNames{puck},'\DGE.csv'],'ReadVariableNames',true);

    %NORMAL CODE:
    IlluminaBarcodes=DGETable.Properties.VariableNames;
    IlluminaBarcodes=IlluminaBarcodes(2:end);
    GeneNames=string(DGETable.Var1);
    DGE=table2array(DGETable(:,2:end));
%    DGE=csvread([IlluminaFolderRoot,PuckNames{puck},'\DGE.csv'],1,1);

    DGEallsums=sum(DGE,1);
    IlluminaBarcodesToAnalyze=min(length(find(DGEallsums>=IlluminaReadThreshold)),MaximumBarcodesToAnalyze);

    [sorted,sortindices]=sort(DGEallsums,'descend');
    DGESorted=DGE(:,sortindices(1:IlluminaBarcodesToAnalyze));
    IlluminaBarcodesDownsampled=IlluminaBarcodes(sortindices(1:IlluminaBarcodesToAnalyze));
    
    
    load([OutputFolders{puck},'AnalysisOutputs-selected'],'Bead');
    Illumina=MapLocationsFunTruseq7UP7(OutputFolders{puck},Bead,IlluminaBarcodesDownsampled);
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
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    
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
    save([OutputFolders{puck},'MappedBeads.mat'],'MappedBeads','-v7.3')
    save([OutputFolders{puck},'BijectiveMapping.mat'],'UniqueMappedBeads','UniqueMappedDGE','UniqueMappedIlluminaBarcodes','-v7.3');
    
    %% Plot figures

    DGEsums=sum(UniqueMappedDGE,1);
    figure(10)
    histogram(DGEsums,1:10:1000)
    title('Total Reads per Barcode for Bijectively Mapped Barcodes')
    set(gca,'yscale','log');
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    
    figure(7)
    clf
    for k=1:length(DGEsums)
        rectangle('Position',[UniqueMappedBeads(k).Locations(1),UniqueMappedBeads(k).Locations(2),5*DGEsums(k)/mean(DGEsums),5*DGEsums(k)/mean(DGEsums)],...
          'Curvature',[1,1], 'FaceColor','r')
    end
    title('Reads per bead for bijective barcodes')
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    
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
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    
    
    
    
    figure(19)
    histogram([MappedBeads.NumIlluminaBarcodes],0:1:10);
    title('Number of Illumina Barcodes per SOLiD barcode')
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');

    
    pixelsperbead=cellfun(@(x) length(x),{MappedBeads.Pixels});
    pixelsperbeadbij=cellfun(@(x) length(x),{UniqueMappedBeads.Pixels});    
    figure(15)
    clf
    histogram(pixelsperbead)
    title('Pixels per Unique SOLiD barcode')
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');    
    
    figure(16)
    clf
    histogram(pixelsperbeadbij)
    title('Pixels per Bijective Pair')
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');    
    
    figure(20)
    histogram([UniqueMappedBeads.HammingDistances],0:1:14);
    title('Distribution of HD between Bead and Illumina Barcodes for Bijective Pairs')
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    
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
    %export_fig([OutputFolders{puck},'Report.pdf'],'-append');
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
    %export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    for qr=1:length(badbeads)
        BeadImageBadBeads(MappedBeads(badbeads(qr)).Pixels)=true;
    end
    BeadImageFinal=zeros(ROIHeight,ROIWidth,3);
    BeadImageFinal(:,:,2)=BeadImageMatchedBeads;
    BeadImageFinal(:,:,1)=BeadImageBadBeads;    
    imshow(BeadImageFinal)
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');

    
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
    %export_fig([OutputFolders{puck},'Report.pdf'],'-append');
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
    %export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    for qr=1:min(2000,length(SortedUniqueMappedBeads))
        BeadImageTopBeads(SortedUniqueMappedBeads(qr).Pixels)=true;
    end
    title('Top 2000 Beads matched to Illumina Barcodes versus 2000 random unmatched beads')
    BeadImageRandomBadBeads=false(ROIHeight,ROIWidth);
    %export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    for qr=1:min(2000,length(badbeads))
        BeadImageRandomBadBeads(MappedBeads(badbeads(qr)).Pixels)=true;
    end
    BeadImageFinal2=zeros(ROIHeight,ROIWidth,3);
    BeadImageFinal2(:,:,2)=BeadImageTopBeads;
    BeadImageFinal2(:,:,1)=BeadImageRandomBadBeads;    
    imshow(BeadImageFinal2)
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');
    
    %we now just want to plot
    
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
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');

    

    %% Plots of the base balance and errors in the Illumina-SOLiD mapping.
    disp('NOTE: The following analysis depends on the number of ligations, the ligation sequence, and the bead structure, and may require modification.')
    %We want to make a plot of the color space balance for the final
    %mapped barcodes
    BaseBalanceBarcodes=[UniqueMappedBeads.Barcodes];
    %The base 5 representations of the basecalls are:
    BaseBalanceBase5Barcodes=cellfun(@(x) reverse(string(x)),{dec2base(BaseBalanceBarcodes,5,14)},'UniformOutput',false);
    BaseBalanceBase5Barcodes=BaseBalanceBase5Barcodes{1};

    BaseBalanceMatrix=zeros(5,14);
    for jp=1:14
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
    bar(BaseBalanceMatrix')
    title('Base balance per ligation for bijectively mapped barcodes - SOLiD');
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');


    colorspacesequence={};
    for seqnum=1:length(UniqueMappedIlluminaBarcodes)
        seq=char(UniqueMappedIlluminaBarcodes(seqnum));
        [PrimerNcolorspace,badflagN]=bs2cs(seq(1:6),PrimerNLigationSequence,'T','T');
        [PrimerUPcolorspace,badflagUP]=bs2cs(seq(7:13),PrimerUPLigationSequence,'A','');
        tmptest=num2str(cat(2,PrimerNcolorspace,PrimerUPcolorspace));
        colorspacesequence{seqnum}=tmptest(~isspace(tmptest));
    end
    colorspacesequence=string(colorspacesequence);
    %The base 5 representations of the basecalls are:
%    BaseBalanceBase5Barcodes=cellfun(@(x) reverse(string(x)),{dec2base(BaseBalanceBarcodes,5,14)},'UniformOutput',false);
%    BaseBalanceBase5Barcodes=BaseBalanceBase5Barcodes{1};

    BaseBalanceMatrix2=zeros(5,14);
    for jp=1:14
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
    bar(BaseBalanceMatrix2')
    title('Base balance per ligation for bijective barcodes expected from Illumina');
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');    
    
    figure(28)
    bar(BaseBalanceMatrix2'-BaseBalanceMatrix')
    title('Base balance of Illumina minus base balance of SOLiD in color space');
    %This plot shows BIAS: it does not show number of errors. A 0 in this
    %plot indicates that there is no bias, not that there are no errors.
    export_fig([OutputFolders{puck},'Report.pdf'],'-append');    
    
    %this plot shows the actual errors. If there is disagreement between
    %Illumina and SOLiD, we plot the SOLiD barcode.
    ErrorMatrix=zeros(5,14);
    ErrorBases=char(colorspacesequence')~=char(BaseBalanceBase5Barcodes);
    ErrorMatrix(1,:)=sum(char(colorspacesequence')==testcmp0 & ErrorBases,1);
    ErrorMatrix(2,:)=sum(char(colorspacesequence')==testcmp1 & ErrorBases,1);
    ErrorMatrix(3,:)=sum(char(colorspacesequence')==testcmp2 & ErrorBases,1);
    ErrorMatrix(4,:)=sum(char(colorspacesequence')==testcmp3 & ErrorBases,1);
    ErrorMatrix(5,:)=sum(char(colorspacesequence')==testcmp4 & ErrorBases,1);
    figure(29)
    bar(ErrorMatrix')
    title('SOLiD Color Recorded on Mismatched Ligations');

    %% Output barcodes to CSV in for quality checking afterwards
    %NOTE: All barcodes outputted this way are in SEQUENCE ORDER, NOT IN
    %ORDER BY LIGATION
    %CHECK THAT THESE ARE NOT IN ORDER BY LIGATION AND ALL OTHER BARCODES
    %ARE!!

    %Bijective Pairs: SOLiD side
    UniqueBeadBarcodesForExport=char(replace(BaseBalanceBase5Barcodes,{'0','1','2','3','4'},{'N','B','G','O','R'}));
    UniqueBeadBarcodesForExport=UniqueBeadBarcodesForExport(:,InverseLigationSequence);
    %Bijective pairs: Illumina side
    UniqueIlluminaBarcodesForExport=char(replace(colorspacesequence',{'0','1','2','3','4'},{'N','B','G','O','R'}));
    UniqueIlluminaBarcodesForExport=UniqueIlluminaBarcodesForExport(:,InverseLigationSequence);
	    
    %All SOLiD beads passing filter
    AllBaseBalanceBarcodes=[Bead.Barcodes];
    %The base 5 representations of the basecalls are:
    BeadBarcodesForExport=cellfun(@(x) reverse(string(x)),{dec2base(AllBaseBalanceBarcodes,5,14)},'UniformOutput',false);
    BeadBarcodesForExport=BeadBarcodesForExport{1};
    BeadBarcodesForExport=char(replace(BeadBarcodesForExport,{'0','1','2','3','4'},{'N','B','G','O','R'}));
    BeadBarcodesForExport=BeadBarcodesForExport(:,InverseLigationSequence);
    
    %All Illumina Barcodes >>> these are in order of the DGE
    Nseq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';
    allilluminacolorspacesequence={};
    for seqnum=1:length(IlluminaBarcodes)
        seq=char(IlluminaBarcodes(seqnum));
        if length(seq)<13
            seq=Nseq(1:13);
        end
        [PrimerNcolorspace,badflagN]=bs2cs(seq(1:6),PrimerNLigationSequence,'T','T');
        [PrimerUPcolorspace,badflagUP]=bs2cs(seq(7:13),PrimerUPLigationSequence,'A','');
        tmptest=num2str(cat(2,PrimerNcolorspace,PrimerUPcolorspace));
        allilluminacolorspacesequence{seqnum}=tmptest(~isspace(tmptest));
    end
%    isemptyfn=cell2mat(cellfun(@(x) ~isempty(x),allilluminacolorspacesequence,'UniformOutput',false));
    allilluminacolorspacesequence=string(allilluminacolorspacesequence);
    IlluminaBarcodesForExport=char(replace(allilluminacolorspacesequence',{'0','1','2','3','4'},{'N','B','G','O','R'}));
    IlluminaBarcodesForExport=IlluminaBarcodesForExport(:,InverseLigationSequence);
    
    save([OutputFolders{puck},'ReadableSequences.mat'],'UniqueBeadBarcodesForExport','UniqueIlluminaBarcodesForExport','BeadBarcodesForExport','IlluminaBarcodesForExport','-v7.3')
    csvwrite([OutputFolders{puck},'BijectiveBeadBarcodes.csv'],UniqueBeadBarcodesForExport);
    csvwrite([OutputFolders{puck},'BijectiveIlluminaBarcodes.csv'],UniqueIlluminaBarcodesForExport);
    csvwrite([OutputFolders{puck},'AllBeadBarcodes.csv'],BeadBarcodesForExport);
    csvwrite([OutputFolders{puck},'AllIlluminaBarcodes.csv'],IlluminaBarcodesForExport);
    
    %% Print output
    fileid=fopen([OutputFolders{puck},'Metrics.txt'],'a');
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
    mkdir(['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',starttimereadable]);
    copyfile(OutputFolders{puck},['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',starttimereadable]);
end

%% Evaluate

%We can now use the Hough transform, and ask what percentage of Hough beads
%also have barcodes, as a way of analyzing what percentage of the beads
%were called. We call the Hough transform and find the centroids, and then
%ask about the fraction of Hough beads with centroids within a radius of a sequenced
%barcode centroid

%We want to produce a plot with the fraction of illumina barcodes that are direct
%matches with a surface bead; one base away from a surface bead; two bases,
%etc., and also for the negative control