%%%%SETUP:
%ONLY TRUE IF USING MIJI:
%1) Before using this function, you must set the MATLAB_JAVA environment
%variable to point to java 8. See also
%https://www.mathworks.com/matlabcentral/answers/130359-how-do-i-change-the-java-virtual-machine-jvm-that-matlab-is-using-on-windows
%For this to work, you must apparently be an administrator on your system,
%and you must run Matlab as an administrator.

%2) If your pucks are a size other than 7x7, you haveto make changes
%both in the yposition and xposition lines, and in the "command" line

%ONLY TRUE IF USING MIJI:
%3) To be able to stitch the pucks from matlab, you need to increase
%matlab's available Java heap space. Within matlab, go to
%Home>Preferences>Matlab>General>Java Heap Memory. 3gb should be fine for
%7x7 images.

%4) Currently, matlab does not make any of the required directories itself,
%although we could easily make them. You have to make the
%C:\users\sgr\pucktmp\PUCKNAME directory, and the C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\InputFolder directory
%And also the C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\PuckBarcodes directory to which BeadSeq will output

%ONLY TRUE IF USING MIJI:
%5) You must download the jheapcl matlab package and put it in the same
%directory as this function is in.

%6) When you are done with the code, you should make the raw data folder in
%Pucks and the InputFolder in find_roi online only via smart sync to save
%hard drive space. You should also delete the pucktmp directory.

%7) A number of paths are hardcoded in the code currently, e.g.
%"C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\helpers\vlfeat-0.9.20\toolbox\"
%in find_roi_stack_fun

%8) The size of the images is curretly hardcoded into find_roi_stack_fun,
%so if you change the size of the images you will have to change it in that
%function also. Look for 'PixelRegion'. The images we use currently should
%be in the range of 1:10660 pixels in each dimension.

%9) The sequencing data with the bead barcodes in it should be in the
%OutputDirectory


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

BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14]; %this determines how the numerical barcode is built from the ligations. Basically, the constant bases should be 0, and all other bases should be consecutive
FolderWithRawTiffs='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Raw\170822 - Pucks 8-6 and 8-7\';
FolderWithProcessedTiffs='C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Processed\';
tmpfolder='C:\Users\sgr\pucktmp\';
IndexFiles={'primers_n_n-1_n-3_up_up-1_','primers_n-2_n-4_','primers_up-2_up-3_up-4'}; %give the prefixes of each of the files. The next character after should be 't',
LigationToIndexFileMapping=[1,1,1,1,2,2,1,1,2,2,1,1,1,1,3,3,3,3,3,3];%for ligations 1:20, which file number are the ligations found in?
%for ligations 1:20, which value of t are they, within their file? This
%could also be deduced from the LigationToIndexFileMapping array
tnumMapping=[1,2,3,4,1,2,5,6,3,4,7,8,9,10,1,2,3,4,5,6];
PuckNames={'Puck_86','Puck_87'}; %give the names of the pucks
OutputFolders={'C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_86\','C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Barcodes\Puck_87\'};
SequencingFileName={'fei_S3_L001_R1_001.fastq'};
TopBarcodeNum=100000; %when it comes to analyzing barcodes from the illumina data, we only consider the TopBarcodeNum barcodes with the most reads. This should ruoghly be the number of beads on the puck.



for puck=1:length(PuckNames)
    ProcessedImageFolders{puck}=[FolderWithProcessedTiffs,PuckNames{puck},'\'];
end
OutputXIndices=[[11,17];[1,7]]; %for each puck, give the X indices where it begins and ends

display('Renaming Files');

%% Move files
for puck=1:length(PuckNames)
    for ligation=1:20
        tnum=tnumMapping(ligation);
        for yposition=1:7
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
            end
            for xposition=OutputXIndices(puck,1):OutputXIndices(puck,2)
                xpositionindex=find(OutputXIndices(puck,1):OutputXIndices(puck,2)==xposition);
                filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',pad(num2str(tnum),2,'left','0'),yposstring,pad(num2str(xposition),2,'left','0'),'.tif'];
                if ~exist(filename)
                    filename=[FolderWithRawTiffs,IndexFiles{LigationToIndexFileMapping(ligation)},'t',num2str(tnum),yposstring,pad(num2str(xposition),2,'left','0'),'.tif'];                    
                end
                outputfilename=[tmpfolder,PuckNames{puck},'\',PuckNames{puck},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Position_X_',pad(num2str(xpositionindex),2,'left','0'),'_Y_',pad(num2str(yposition),2,'left','0'),'.tif'];
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
    command=replace(replace(replace('type=[Filename defined position] order=[Defined by filename         ] grid_size_x=7 grid_size_y=7 tile_overlap=30 first_file_index_x=01 first_file_index_y=01 directory=DIRECTORYPATH file_names=PUCKNAME_Ligation_LL_Position_X_{xx}_Y_{yy}.tif output_textfile_name=Ligation_LL.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]','LL',pad(num2str(ligation),2,'left','0')),'PUCKNAME',thispuckname),'DIRECTORYPATH',thisinputfolder);
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
    find_roi_stack_fun(BaseName,Suffix,10400);
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end
%% Bead calling and sequencing
%we now run Bead_Seq itself. Again, where parallelization is non-trivial in
%Bead_Seq, it is implemented naively using parfor.

display('Base Calling')
    
for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
%We now send the command to do the registration with find_roi

    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
	[BeadBarcodes BeadLocations BeadImage]=BeadSeqFun(BaseName,suffix,OutputFolders{puck},2,BarcodeSequence,20)
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end

%% Match Illumina barcodes


%% Evaluate

%We can now use the Hough transform, and ask what percentage of Hough beads
%also have barcodes, as a way of analyzing what percentage of the beads
%were called. We call the Hough transform and find the centroids, and then
%ask about the fraction of Hough beads with centroids within a radius of a sequenced
%barcode centroid

%We want to produce a plot with the fraction of illumina barcodes that are direct
%matches with a surface bead; one base away from a surface bead; two bases,
%etc., and also for the negative control