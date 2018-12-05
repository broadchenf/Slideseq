%This script just identifies barcodes from Illumina sequencing data

%This is currently for the June beads, for which the sequence is JJJJJJTCTTCAGCGTTCCCGAGAJJJJJJNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
clear all
close all

%OutputDirectory='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\BigPuck-S00004-54250303';
%SequencingPath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\BigPuck-S00004-54250303\SamBigPuck_S4_L001_R1_001.fastq';

%OutputDirectory='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\SmallPuck-S00003-54260232';
%SequencingPath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\SmallPuck-S00003-54260232\SamSmallPuck_S3_L001_R1_001.fastq';

OutputDirectory='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\10k-S00005-54247365';
SequencingPath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\10k-S00005-54247365\Sam10k_S5_L001_R1_001.fastq';

%OutputDirectory='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\30k-S00006-54245371';
%SequencingPath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\30k-S00006-54245371\Sam30k_S6_L001_R1_001.fastq';

%OutputDirectory='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\100k-S00007-54236308';
%SequencingPath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\100k-S00007-54236308\Sam100k_S7_L001_R1_001.fastq';

%OutputDirectory='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\300k-S00008-54245370';
%SequencingPath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170811 - Extension oligo debugging\300k-S00008-54245370\Sam300k_S8_L001_R1_001.fastq';


reads=fastqread(SequencingPath);
UPSeq='TCTTCAGCGTTCCCGAGA';
AcceptableMismatches=4;
CountThreshold=2; %UMIs must occur MORE than this number of times to be counted
BarcodeFormat=[1:6,25:30];

%% What is the most common barcode?
if 0
observedBeadBarcodesx = cellfun(@(x) x(BarcodeFormat),...
        {reads.Sequence}, 'UniformOutput', false); % compare first 10 of each seq
    [uniqueBarcodesx, firstOccx, barcodeRefx] = unique(observedBeadBarcodesx, 'stable'); 
    [occCountsx, edges0x] = histcounts(barcodeRefx,1:length(firstOccx)); %histcounts(barcodeRef, 'BinMethod', 'integers'); %How often each barcode occurs
end
    %Note: the most common barcode is uniqueBarcodesx(find(occCountsx==max(occCountsx)));
%is the barcode that occurs the most. Here, it's just Ts
%% How many of the reads are >50% T?
if 0
Ts='TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT';
ReadsThatAreJustT = cell2mat(cellfun(@(x) sum(x==Ts(1:37))>18,{reads.Sequence},'UniformOutput',false));
ReadsThatAreJustTseq=reads(ReadsThatAreJustT);
%Note that most of the reads have >18Ts, and they are generally either
%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAA, or some random string of ~20
%nucleotides followed by Ts.
end
%% Identify the bead barcodes
ReadsThatHaveUP = cell2mat(cellfun(@(x) sum(x(7:24)==UPSeq)>(length(UPSeq)-AcceptableMismatches),{reads.Sequence},'UniformOutput',false));

ReadsThatHaveUPseq=reads(ReadsThatHaveUP);

    
observedBeadBarcodes = cellfun(@(x) x(BarcodeFormat),...
        {ReadsThatHaveUPseq.Sequence}, 'UniformOutput', false); % compare first 10 of each seq
    [uniqueBarcodes, firstOcc, barcodeRef] = unique(observedBeadBarcodes, 'stable'); 
    % uniqueUMIs = cell with unique UMI strings.
    % firstOcc = indices of first occurrence of unique UMI
    % umiRef = indices in uniqueUMIs list that correspond to observedUMIs list. 
    % Note: firstOcc matches indices of previous code's 'nonduplicates'
%    readsdedup = reads(firstOcc); % unique entries.
    
%We now cut out barcodes with less than say 10 reads
%Note that this works because each element occurs at least once. Otherwise we would get an offset because histcounts would delete some bins
    [occCounts, edges0] = histcounts(barcodeRef,1:length(firstOcc)); %histcounts(barcodeRef, 'BinMethod', 'integers'); %How often each barcode occurs
    
    barcode_num = 1:length(occCounts);

    barcode_g1 = barcode_num(occCounts > CountThreshold); %only count barcodes that occur more than 10 times
    readsdedup = ReadsThatHaveUPseq(firstOcc(barcode_g1));

%We now convert into color space

%% Look at Read 2
if 0
Read2SequencingPath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170701 - Puck 3 Extension Oligo\Undetermined_S0_R2_001.fastq';
read2=fastqread(Read2SequencingPath);
%This is mostly As 
end

%% Base space to color space conversion
%This is a super shitty thing we have to do because of the way we're doing
%the ligations

%for UP, we do primer UP ligation 1, then primer UP-1 ligation 1,
% Primer UP-1 ligation 2, primer UP-2 ligation 2, primer UP-3 ligation 2, and
% primer UP-4 ligation 2
PrimerUPLigationSequence=[2, 1, 6, 5, 4, 3];
PrimerNLigationSequence = [2, 7, 1, 6, 5, 4]; 
Base5Barcodes=zeros(1,length(readsdedup));
badflags=0;
for seqnum=1:length(readsdedup)
    [PrimerNcolorspace,badflagN]=bs2cs(readsdedup(seqnum).Sequence(BarcodeFormat(1:6)),PrimerNLigationSequence,'C','T');
    [PrimerUPcolorspace,badflagUP]=bs2cs(readsdedup(seqnum).Sequence(BarcodeFormat(7:12)),PrimerUPLigationSequence,'A','');
    colorspacesequence=cat(2,PrimerNcolorspace,PrimerUPcolorspace);
    if badflagN || badflagUP
        badflags=badflags+1;
    end
    for m=1:9%m=1:length(colorspacesequence) We just do 1:11 right now because puck 3 is missing the last ligation
        Base5Barcodes(seqnum)=Base5Barcodes(seqnum)+uint32(colorspacesequence(m))*5^(m-1);
    end
end
%For Primer N, we do Primer N ligation 1, Primer N ligation 2, Primer N-1 ligation 1, Primer N-1 ligation 2, Primer N-2 ligation 2, Primer N-3 ligation 2 
save([OutputDirectory,'Base5Barcodes'],'Base5Barcodes');


%% Negative control
%We're generating random barcodes to see what overlap we get:

NegBarcodes=zeros(1,length(readsdedup));
for negseq=1:length(readsdedup)
    for m=1:9
        NegBarcodes(negseq)=NegBarcodes(negseq)+uint32(ceil(4*rand()))*5^(m-1);
    end
end