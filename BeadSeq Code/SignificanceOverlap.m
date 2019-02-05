function [Cooccurrences,CorrelatingGenes,Anticooccurrences,AnticorrelatingGenes,GoodGeneList]=SignificanceOverlap(GeneList,PuckDirectories,varargin)
%NOTE THAT PuckDirectories is a cell of puck names, NOT a cell of puck
%directories.
%GeneList is a String Array.
%For each puck, we're going to analyze the six significant genes, and we're
%going to count genes that frequently (across pucks) correlate with multiple significant
%genes.

%Note that if any of these genes don't appear on the puck, it will crash,
%which is actually okay because it's an insurance policy against using
%super lowly expressed genes.
NumGenes=length(GeneList); %To test


Cooccurrences=false(22000,NumGenes,length(PuckDirectories));
CorrelatingGenes=string();
Anticooccurrences=false(22000,NumGenes,length(PuckDirectories));
AnticorrelatingGenes=string();


    ImageSizeX=6030;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="ImageSizeX"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        ImageSizeX=varargin{index+1};
    end
    ImageSizeY=6030;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="ImageSizeY"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        ImageSizeY=varargin{index+1};
    end


    NumSamples=100;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="NumSamples"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        NumSamples=varargin{index+1};
    end
    ZScoreCutoff=3;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="ZScoreCutoff"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        ZScoreCutoff=varargin{index+1};
    end
    BeadCutoff=4;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadCutoff"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadCutoff=varargin{index+1};
    end
    ClusterToAnalyze=0;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="ClusterToAnalyze"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        ClusterToAnalyze=varargin{index+1};
    end 
    FilterByCutoff=0;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="FilterByCutoff"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        FilterByCutoff=varargin{index+1};
    end 
    
    Spread=5;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Spread"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Spread=varargin{index+1};
    end
    
    DownsampleFactor=10;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="DownsampleFactor"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        DownsampleFactor=varargin{index+1};
    end
    
for puck=1:length(PuckDirectories)

PuckDirectory=GetPuckDirectory(PuckDirectories{puck});

BeadMappingFile=FindMostRecentMapping(PuckDirectory);

    if length(ClusterToAnalyze)==1 && ClusterToAnalyze==0
        [GeneNames,UniqueMappedBeads,UniqueMappedDGE,UniqueMappedBarcodes]=LoadBijectiveMapping(PuckName);
    elseif FilterByCutoff==0
        [UniqueMappedDGE,UniqueMappedBeads,UniqueMappedBarcodes,GeneNames]=DGEByCluster(PuckDirectory,ClusterToAnalyze,'BeadMappingFile',BeadMappingFile);
    elseif FilterByCutoff>0
        [UniqueMappedDGE,UniqueMappedBeads,UniqueMappedBarcodes,GeneNames]=DGEByCutoff(PuckDirectories{puck},ClusterToAnalyze,FilterByCutoff,'BeadMappingFile',BeadMappingFile);        
    end


%PDFOutputName=fullfile(PuckDirectory,BeadMappingFile,'J20_PlaqueGenes_OverlapOutput.pdf');
%First step is to make a sparse vector, where each column is a blurred
%bead
[GoodGeneList,GoodGeneIndices]=intersect(GeneNames,GeneList);
assert(length(GoodGeneList)==NumGenes) %we don't currently handle the case where this isn't true.

DGEBinary=UniqueMappedDGE>0;


BeadsPerGene=sum(DGEBinary,2);



%% Pixels by significant gene
currentbead=0;
DGEBinaryReduced=DGEBinary(GoodGeneIndices,:);
XIndices=zeros(1,nnz(DGEBinaryReduced)*(2*Spread+1)^2);
YIndices=zeros(1,nnz(DGEBinaryReduced)*(2*Spread+1)^2);

for jj=1:length(GoodGeneList)
    j=GoodGeneIndices(jj);
    listofbeads=find(DGEBinary(j,:));
    for i=1:length(listofbeads)
        currentbead=currentbead+1;
        %This doesn't work, we're hitting the edge still
        %tx=repmat((max(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)-Spread,1)):(min(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)+Spread,floor(ImageSizeX/DownsampleFactor))),2*Spread+1,1);
        %ty=repmat((max(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)-Spread,1)):(min(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)+Spread,floor(ImageSizeY/DownsampleFactor))),2*Spread+1,1)';

        tx=repmat((round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)-Spread):(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)+Spread),2*Spread+1,1);
        ty=repmat((round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)-Spread):(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)+Spread),2*Spread+1,1)';
        
        %We now restrict tx and ty to be on the image
        GoodIndices=tx<floor(ImageSizeX/DownsampleFactor) & tx >0 & ty > 0 & ty <floor(ImageSizeY/DownsampleFactor);
        tx=tx(GoodIndices);
        ty=ty(GoodIndices);
        
    %    a=toc(test1)
    %

        linearind=sub2ind([round(ImageSizeX/DownsampleFactor),round(ImageSizeY/DownsampleFactor)],tx(:),ty(:));
    
        XIndices(((currentbead-1)*(2*Spread+1)^2+1):((currentbead-1)*(2*Spread+1)^2+length(linearind)))=linearind;
        if length(linearind)<(2*Spread+1)^2
            XIndices((currentbead-1)*(2*Spread+1)^2+length(linearind):(currentbead*(2*Spread+1)^2))=-1;
        end
        YIndices(((currentbead-1)*(2*Spread+1)^2+1):(currentbead*(2*Spread+1)^2))=jj;
    end
    
    %    im=imgaussfilt(im,Spread); %Doing a Gaussian filter is too hard.
                                %Could also try rectangle function
%    test2=tic();
%    PixelsByBead(j,:)=reshape(im,[1,ImageSizeX*ImageSizeY]);%im(:); %this is slow why
%    b=toc(test2)
end
XIndices=XIndices(XIndices>0);
YIndices=YIndices(XIndices>0);

PixelsByGoodGene=sparse(XIndices,YIndices,ones(size(XIndices)),round(ImageSizeX/DownsampleFactor)*round(ImageSizeY/DownsampleFactor),size(DGEBinaryReduced,1));


%% Pixels by gene 
%In the V3 version of this function, we calculate the correlation between
%the significant genes and all other genes.
currentbead=0;
XIndices=zeros(1,nnz(DGEBinary)*(2*Spread+1)^2);
YIndices=zeros(1,nnz(DGEBinary)*(2*Spread+1)^2);

for j=1:size(DGEBinary,1)
    if j/100==floor(j/100)
        j
    end
    listofbeads=find(DGEBinary(j,:));
    for i=1:length(listofbeads)
        currentbead=currentbead+1;
        %This doesn't work, we're hitting the edge still
        %tx=repmat((max(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)-Spread,1)):(min(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)+Spread,floor(ImageSizeX/DownsampleFactor))),2*Spread+1,1);
        %ty=repmat((max(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)-Spread,1)):(min(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)+Spread,floor(ImageSizeY/DownsampleFactor))),2*Spread+1,1)';

        tx=repmat((round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)-Spread):(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)+Spread),2*Spread+1,1);
        ty=repmat((round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)-Spread):(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)+Spread),2*Spread+1,1)';
        
        %We now restrict tx and ty to be on the image
        GoodIndices=tx<floor(ImageSizeX/DownsampleFactor) & tx >0 & ty > 0 & ty <floor(ImageSizeY/DownsampleFactor);
        tx=tx(GoodIndices);
        ty=ty(GoodIndices);
        
    %    a=toc(test1)
    %

        linearind=sub2ind([round(ImageSizeX/DownsampleFactor),round(ImageSizeY/DownsampleFactor)],tx(:),ty(:));
    
        XIndices(((currentbead-1)*(2*Spread+1)^2+1):((currentbead-1)*(2*Spread+1)^2+length(linearind)))=linearind;
        if length(linearind)<(2*Spread+1)^2
            XIndices((currentbead-1)*(2*Spread+1)^2+length(linearind):(currentbead*(2*Spread+1)^2))=-1;
        end
        YIndices(((currentbead-1)*(2*Spread+1)^2+1):(currentbead*(2*Spread+1)^2))=j;
    end
    
    %    im=imgaussfilt(im,Spread); %Doing a Gaussian filter is too hard.
                                %Could also try rectangle function
%    test2=tic();
%    PixelsByBead(j,:)=reshape(im,[1,ImageSizeX*ImageSizeY]);%im(:); %this is slow why
%    b=toc(test2)
end
XIndices=XIndices(XIndices>0);
YIndices=YIndices(XIndices>0);

PixelsByGene=sparse(XIndices,YIndices,ones(size(XIndices)),round(ImageSizeX/DownsampleFactor)*round(ImageSizeY/DownsampleFactor),size(DGEBinary,1));

GeneInnerProducts=PixelsByGoodGene'*PixelsByGene;

%% Here we actually try to generate the null model. We want two A x A matrices, where A is the number of genes. The matrices are the expected overlap and standard devitaion 
%For consistency, we use identical code to what is above
%The insight here is that the time consuming part is making the sparse
%matrix, so we generate random samples for all genes at once and then
%collate them into a single sparse matrix. We could in principle generate
%many samples this way before making the sparse matrix. But this is okay
%for now.
NumReadsPerBead=sum(UniqueMappedDGE,1);
ProbabilityPerBead=NumReadsPerBead/sum(NumReadsPerBead);
MeanMatrix=cell(1,length(GoodGeneList));
StdMatrix=cell(1,length(GoodGeneList));
RandomSampleInnerProducts=zeros(length(GoodGeneList),size(DGEBinary,1),NumSamples);

for sample=1:NumSamples
    sample
    currentbead=0;
    XIndices=zeros(1,nnz(DGEBinary)*(2*Spread+1)^2);
    YIndices=zeros(1,nnz(DGEBinary)*(2*Spread+1)^2);

    for j=1:size(DGEBinary,1)
        NumBeads=sum(DGEBinary(j,:));
%        if j/100==floor(j/100)
%            j
%        end
        listofbeads=datasample(1:length(NumReadsPerBead),NumBeads,'Weights',ProbabilityPerBead,'Replace',false);
        for i=1:length(listofbeads)
            currentbead=currentbead+1;
            %This doesn't work, we're hitting the edge still
            %tx=repmat((max(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)-Spread,1)):(min(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)+Spread,floor(ImageSizeX/DownsampleFactor))),2*Spread+1,1);
            %ty=repmat((max(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)-Spread,1)):(min(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)+Spread,floor(ImageSizeY/DownsampleFactor))),2*Spread+1,1)';

            tx=repmat((round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)-Spread):(round(UniqueMappedBeads(listofbeads(i)).Locations(1)/DownsampleFactor)+Spread),2*Spread+1,1);
            ty=repmat((round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)-Spread):(round(UniqueMappedBeads(listofbeads(i)).Locations(2)/DownsampleFactor)+Spread),2*Spread+1,1)';

            %We now restrict tx and ty to be on the image
            GoodIndices=tx<floor(ImageSizeX/DownsampleFactor) & tx >0 & ty > 0 & ty <floor(ImageSizeY/DownsampleFactor);
            tx=tx(GoodIndices);
            ty=ty(GoodIndices);

        %    a=toc(test1)
        %

            linearind=sub2ind([round(ImageSizeX/DownsampleFactor),round(ImageSizeY/DownsampleFactor)],tx(:),ty(:));

            XIndices(((currentbead-1)*(2*Spread+1)^2+1):((currentbead-1)*(2*Spread+1)^2+length(linearind)))=linearind;
            if length(linearind)<(2*Spread+1)^2
                XIndices((currentbead-1)*(2*Spread+1)^2+length(linearind):(currentbead*(2*Spread+1)^2))=-1;
            end
            YIndices(((currentbead-1)*(2*Spread+1)^2+1):(currentbead*(2*Spread+1)^2))=j;
        end

        %    im=imgaussfilt(im,Spread); %Doing a Gaussian filter is too hard.
                                    %Could also try rectangle function
    %    test2=tic();
    %    PixelsByBead(j,:)=reshape(im,[1,ImageSizeX*ImageSizeY]);%im(:); %this is slow why
    %    b=toc(test2)
    end
    XIndices=XIndices(XIndices>0);
    YIndices=YIndices(XIndices>0);

    PixelsByRandomGene=sparse(XIndices,YIndices,ones(size(XIndices)),round(ImageSizeX/DownsampleFactor)*round(ImageSizeY/DownsampleFactor),size(DGEBinary,1));

    
    RandomSampleInnerProducts(:,:,sample)=PixelsByGoodGene'*PixelsByRandomGene;
end


RealMeanMatrix=mean(RandomSampleInnerProducts,3);
RealStdMatrix=std(RandomSampleInnerProducts,[],3);


ZScores=(GeneInnerProducts-RealMeanMatrix)./RealStdMatrix;
ZScores(isnan(ZScores))=0;
ZScores(isinf(ZScores))=0;

%% Clustering
%We only consider POSITIVE correlations, because inverse correlations will
%cluster together, which doesn't make sense
ZScoresBinary=double(ZScores>5);%-double(ZScores<-3);

imagesc(ZScoresBinary);

ZScorePassNumber=sum(ZScoresBinary,2)/size(ZScoresBinary,2);
ZScoresBinary(:,ZScorePassNumber>0.5)=0; %We cut out a few genes that correlate highly with everything.

%% We now want to form metagenes
if 0
for Gene=1:length(GoodGeneList)
    GenesForMetagene=find(ZScoresBinary(Gene,:));
    CountsPerBead=sum(DGEBinary(GenesForMetagene,:),1);
    PositiveBeadImage=zeros(ImageSizeX,ImageSizeY);
    for qr=1:length(UniqueMappedBeads)
        PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-10)):min(ImageSizeX,(round(UniqueMappedBeads(qr).Locations(1))+10)),max(1,(round(UniqueMappedBeads(qr).Locations(2)-10))):min(ImageSizeY,round(UniqueMappedBeads(qr).Locations(2)+10)))=CountsPerBead(qr);
    end
%    PositiveBeadImage=PositiveBeadImage/sum(sum(PositiveBeadImage));
    imagesc(PositiveBeadImage)
    title(['Metagene for gene ',GoodGeneList{Gene}]);
    export_fig(PDFOutputName,'-append')        
end
end


for j=1:NumGenes
    CorrGenesThisPuck=GeneNames(find(ZScores(j,:)>ZScoreCutoff & sum(DGEBinary,2)'>BeadCutoff));
    for p=1:length(CorrGenesThisPuck)
        [corgene,corindex]=intersect(CorrelatingGenes,CorrGenesThisPuck(p));
        if ~isempty(corgene)
            Cooccurrences(corindex,j,puck)=true;
        else
            if CorrelatingGenes(1)==""
                ind=1;
            else
                ind=length(CorrelatingGenes)+1;
            end
            CorrelatingGenes(ind)=CorrGenesThisPuck(p);
            Cooccurrences(ind,j,puck)=true;
        end 
    end
end

for j=1:NumGenes
    AntiCorrGenesThisPuck=GeneNames(find(ZScores(j,:)<-ZScoreCutoff & sum(DGEBinary,2)'>BeadCutoff));
    for p=1:length(AntiCorrGenesThisPuck)
        [corgene,corindex]=intersect(AnticorrelatingGenes,AntiCorrGenesThisPuck(p));
        if ~isempty(corgene)
            Anticooccurrences(corindex,j,puck)=true;
        else
            if AnticorrelatingGenes(1)==""
                ind=1;
            else
                ind=length(AnticorrelatingGenes)+1;
            end
            AnticorrelatingGenes(ind)=AntiCorrGenesThisPuck(p);
            Anticooccurrences(ind,j,puck)=true;
        end 
    end
end


end

Cooccurrences=Cooccurrences(1:length(CorrelatingGenes),:,:);
Anticooccurrences=Anticooccurrences(1:length(AnticorrelatingGenes),:,:);

%These are genes that correlate with at least two of the 6 test genes on at
%least 2 of the pucks.
%CorrelatingGenes(find(sum(sum(Cooccurrences,2)>=2,3)>=2)) %for each gene, we are counting how many times it correlates with at least 2 of the 6 significant genes on the puck. We find the genes that correlate with at least 2 significant genes on at least 2 pucks.

%AnticorrelatingGenes(find(sum(sum(Anticooccurrences,2)>=2,3)>=2))