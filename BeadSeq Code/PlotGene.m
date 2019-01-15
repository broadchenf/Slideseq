function BeadImage=PlotGene(GeneNum,UniqueMappedDGE,UniqueMappedBeads,varargin)

    PlotStyle="Circles";
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="PlotStyle"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        PlotStyle=varargin{index+1};
    end
    
    BeadCutoff=30; %this is only used for "Normalize" atm.
    Normalize=0;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Normalize"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Normalize=1;
    end
    
    
    BeadSizeFactor=1;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadSizeFactor"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadSizeFactor=varargin{index+1};
    end
    
    FigNum=25;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="FigNum"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        FigNum=varargin{index+1};
    end
    Overlay=1;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Overlay"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Overlay=varargin{index+1};
    end
    OverlayDownsample=5; %for the background beads in overlay, we display one bead out of every x, where this is x.
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="OverlayDownsample"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        OverlayDownsample=varargin{index+1};
    end
    
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Loadings"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Loadings=varargin{index+1};
    end

    
    Cluster=0;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Cluster"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Cluster=varargin{index+1};
    end
    %If cluster>0, then you also have to give a puck name.
    PuckName="None";
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="PuckName"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        PuckName=char(varargin{index+1});
    end
    if (length(Cluster)>1 || Cluster>0) && PuckName=="None"
        disp('If Cluster>0, you must specify a puck name.')
        assert(1==0)
    end

    
    %You are given a full DGE, and asked to plot on a subset.
    if length(Cluster)>1 || Cluster >0 %You now have to run DGEFromCluster, but without loading the DGE in.
        PuckDirectory=GetPuckDirectory(PuckName);
        BeadMappingFile=FindMostRecentMapping(PuckDirectory);
        ClusterPath=fullfile(PuckDirectory,BeadMappingFile,'AnalogizerClusterAssignments.csv');
        opts=detectImportOptions(ClusterPath);
        opts=setvartype(opts,'x','double');
        ClusterAssignments=readtable(ClusterPath,opts,'ReadVariableNames',true);
        ClusterAssignments.Properties.VariableNames={'Barcode','Cluster'};

        load(fullfile(PuckDirectory,BeadMappingFile,'BijectiveMapping.mat'),'UniqueMappedIlluminaBarcodes')
        %NOTE THAT WE ARE FILTERING HERE: we have thrown out slideseq beads with
        %<100 reads in the R program
        [C,ia,ib]=intersect(ClusterAssignments.Barcode,UniqueMappedIlluminaBarcodes);
        assert(length(C)==length(ClusterAssignments.Barcode))
        ClusterAssignments=ClusterAssignments(ia,:);
        UniqueMappedDGE=UniqueMappedDGE(:,ib);
        UniqueMappedBeads=UniqueMappedBeads(ib);
        %We now need to make sure the ClusterAssignments is in the same order as
        %the DGE
        BeadsInCluster=ismember(ClusterAssignments.Cluster,Cluster);
        UniqueMappedDGE=UniqueMappedDGE(:,BeadsInCluster);
        UniqueMappedBeads=UniqueMappedBeads(BeadsInCluster);
    end

    ImageSize=6030;

    MarkerSum=sum(UniqueMappedDGE(GeneNum,:),1);
    TotalSum=sum(UniqueMappedDGE,1);
    %for qr=1:length(UniqueMappedBeads)
    figure(FigNum)
    clf
    if PlotStyle=="Circles"
        if Overlay
            for qr=1:length(UniqueMappedBeads)
                if MarkerSum(qr)==0 && floor(qr/OverlayDownsample)==qr/OverlayDownsample
                    rectangle('Position',[UniqueMappedBeads(qr).Locations(1),UniqueMappedBeads(qr).Locations(2),35,35],...
                      'Curvature',[1,1], 'FaceColor','g','EdgeColor','None')
                end
            end
            
        end
        
        if ~Normalize
            for qr=1:length(UniqueMappedBeads)
                if MarkerSum(qr)>0
                    rectangle('Position',[UniqueMappedBeads(qr).Locations(1),UniqueMappedBeads(qr).Locations(2),min(300,50*MarkerSum(qr)*BeadSizeFactor),min(300,50*MarkerSum(qr)*BeadSizeFactor)],...
                      'Curvature',[1,1], 'FaceColor','b')
                end
            end
            
        else
            for qr=1:length(UniqueMappedBeads)
                if MarkerSum(qr)>0 && TotalSum(qr)>=BeadCutoff
                    rectangle('Position',[UniqueMappedBeads(qr).Locations(1),UniqueMappedBeads(qr).Locations(2),min(300,5000*MarkerSum(qr)*BeadSizeFactor/TotalSum(qr)),min(300,5000*MarkerSum(qr)*BeadSizeFactor/TotalSum(qr))],...
                      'Curvature',[1,1], 'FaceColor','b')
                end
            end
        end
        BeadImage=0;
    end
    
    if PlotStyle=="Imagesc"
        PositiveBeadImage=zeros(ImageSize,ImageSize);
        for qr=1:length(UniqueMappedBeads)
            if MarkerSum(qr)>0
                PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))+MarkerSum(qr);
%                PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))+MarkerSum(qr);
            end
        end
        BeadImage=PositiveBeadImage;
        imagesc(BeadImage);
        colormap(jet)
    end
    if PlotStyle=="ImagescFeiGene"
        PositiveBeadImage={}
        for qr=1:length(UniqueMappedBeads)
            if MarkerSum(qr)>0
                %PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))+MarkerSum(qr);
                PositiveBeadImage{qr} = [MarkerSum(qr) UniqueMappedBeads(qr).Locations(1) UniqueMappedBeads(qr).Locations(2)];
                %                PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))+MarkerSum(qr);
            
            
            
            
            end
        end
        
        Temp = PositiveBeadImage;
        FilteredPositions = cell2mat(Temp(cell2mat(cellfun(@(x) ~isempty(x),Temp,'UniformOutput',false))));

        FilteredPositions=round(reshape(FilteredPositions,3,length(FilteredPositions)/3));
        Xint = round((FilteredPositions(2,:))/4.02)+1;
        Yint = round((FilteredPositions(3,:))/4.02)+1;
        out = zeros(max([1500; 1500]));
        out(sub2ind(size(out),Xint,Yint)) = 1000*FilteredPositions(1,:);

        filterSizes = [7 7;11 11;17 17];
  
        maxFilterSize = max(filterSizes);
        padSize = (maxFilterSize - 1)/2;
        paddedImage = padarray(out,padSize,'replicate','both');
        intImage = integralImage(paddedImage);
        filteredImage1 = integralBoxFilter(intImage, filterSizes(2,:));

        filteredImage1 = imgaussfilt(filteredImage1, 3);
        
        BeadImage=filteredImage1;

    end
    
    if PlotStyle=="ImagescFeiGeneNormalized"
        PositiveBeadImage={}
        TotalSum=sum(sum(UniqueMappedDGE))/size(UniqueMappedDGE,2);
        for qr=1:length(UniqueMappedBeads)
            if MarkerSum(qr)>0
                MarkerSum(qr)=MarkerSum(qr)*TotalSum/sum(UniqueMappedDGE(:,qr));
                %PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))+MarkerSum(qr);
                PositiveBeadImage{qr} = [MarkerSum(qr) UniqueMappedBeads(qr).Locations(1) UniqueMappedBeads(qr).Locations(2)];
                %                PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))+MarkerSum(qr);
            
            
            
            
            end
        end
        
        Temp = PositiveBeadImage;
        FilteredPositions = cell2mat(Temp(cell2mat(cellfun(@(x) ~isempty(x),Temp,'UniformOutput',false))));

        FilteredPositions=round(reshape(FilteredPositions,3,length(FilteredPositions)/3));
        Xint = round((FilteredPositions(2,:))/4.02)+1;
        Yint = round((FilteredPositions(3,:))/4.02)+1;
        out = zeros(max([1500; 1500]));
        out(sub2ind(size(out),Xint,Yint)) = 1000*FilteredPositions(1,:);

        filterSizes = [7 7;11 11;17 17];
  
        maxFilterSize = max(filterSizes);
        padSize = (maxFilterSize - 1)/2;
        paddedImage = padarray(out,padSize,'replicate','both');
        intImage = integralImage(paddedImage);
        filteredImage1 = integralBoxFilter(intImage, filterSizes(2,:));

        filteredImage1 = imgaussfilt(filteredImage1, 3);
        
        BeadImage=filteredImage1;

    end

    %Plot Clusters as circles
     if PlotStyle=="ImagescFeiCluster"
        PositiveBeadImage={}
        for qr=1:length(UniqueMappedBeads)
            if MarkerSum(qr)>0
                %PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))+MarkerSum(qr);
                PositiveBeadImage{qr} = [MarkerSum(qr) UniqueMappedBeads(qr).Locations(1) UniqueMappedBeads(qr).Locations(2)];
                %                PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))+MarkerSum(qr);
            
            
            
            
            end
        end
        
        Temp = PositiveBeadImage;
        FilteredPositions = cell2mat(Temp(cell2mat(cellfun(@(x) ~isempty(x),Temp,'UniformOutput',false))));

        FilteredPositions=round(reshape(FilteredPositions,3,length(FilteredPositions)/3));
        Xint = round((FilteredPositions(2,:))/4.02)+1;
        Yint = round((FilteredPositions(3,:))/4.02)+1;
        out = zeros(max([1500; 1500]));
        out(sub2ind(size(out),Xint,Yint)) = 1000*FilteredPositions(1,:);

         se = offsetstrel('ball',5,5);
         %se = strel('sphere',20);
         dilatedI = imdilate(out,se);

       
        BeadImage=dilatedI;

     end
    
     if PlotStyle=="ImagescFei"
        PositiveBeadImage={}
        for qr=1:length(UniqueMappedBeads)
            if MarkerSum(qr)>0
                %PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))+MarkerSum(qr);
                PositiveBeadImage{qr} = [MarkerSum(qr) UniqueMappedBeads(qr).Locations(1) UniqueMappedBeads(qr).Locations(2)];
                %                PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))+MarkerSum(qr);
            
            
            
            
            end
        end
        
     
        BeadImage=PositiveBeadImage;

     end
     
     
     if PlotStyle=="ImagescClusterByLoading"
        PositiveBeadImage={}
        for qr=1:length(UniqueMappedBeads)
            if MarkerSum(qr)>0

                %PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+35)),max(1,(round(UniqueMappedBeads(qr).Locations(2))-35)):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+35)))+MarkerSum(qr);
                PositiveBeadImage{qr} = [MarkerSum(qr)*Loadings(qr) UniqueMappedBeads(qr).Locations(1) UniqueMappedBeads(qr).Locations(2)];
                %                PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))=PositiveBeadImage(max(1,(round(UniqueMappedBeads(qr).Locations(1))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(1))+25*min(MarkerSum(qr)*BeadSizeFactor,10))),max(1,(round(UniqueMappedBeads(qr).Locations(2))-25*min(MarkerSum(qr)*BeadSizeFactor,10))):min(ImageSize,(round(UniqueMappedBeads(qr).Locations(2))+25*min(MarkerSum(qr)*BeadSizeFactor,10))))+MarkerSum(qr);
            
            
            
            
            end
        end
        
        Temp = PositiveBeadImage;
        FilteredPositions = cell2mat(Temp(cell2mat(cellfun(@(x) ~isempty(x),Temp,'UniformOutput',false))));

        FilteredPositions=round(reshape(FilteredPositions,3,length(FilteredPositions)/3));
        Xint = round((FilteredPositions(2,:))/4.02)+1;
        Yint = round((FilteredPositions(3,:))/4.02)+1;
        out = zeros(max([1500; 1500]));
        out(sub2ind(size(out),Xint,Yint)) = 1000/max(FilteredPositions(1,:))*FilteredPositions(1,:);

         se = offsetstrel('ball',5,5);
         %se = strel('sphere',20);
         dilatedI = imdilate(out,se);

       
        BeadImage=dilatedI;

    end
    
    
    
    if PlotStyle=="Default"
    PositiveBeadImage=false(ImageSize,ImageSize);
    for qr=1:length(UniqueMappedBeads)
        if MarkerSum(qr)>0
            PositiveBeadImage(UniqueMappedBeads(qr).Pixels)=true;
        end
    end

    if Overlay
        NegBeadImage=false(ImageSize,ImageSize);
        for qr=1:length(UniqueMappedBeads)
            if ~MarkerSum(qr)>0 && floor(qr/OverlayDownsample)==qr/OverlayDownsample
                NegBeadImage(UniqueMappedBeads(qr).Pixels)=true;
            end
        end

%        BeadImage=ones(ImageSize,ImageSize,3);
%        BeadImage(:,:,1)=BeadImage(:,:,1)-PositiveBeadImage;
%        BeadImage(:,:,:)=BeadImage(:,:,:)-0.5*NegBeadImage;
 
%        BeadImage=zeros(ImageSize,ImageSize,3)+0.15;
%        BeadImage(:,:,1)=BeadImage(:,:,1)+PositiveBeadImage;
%        BeadImage(:,:,2)=BeadImage(:,:,2)+0.3*NegBeadImage;
%        BeadImage(:,:,3)=BeadImage(:,:,3)+0.3*NegBeadImage;

        BeadImage=zeros(ImageSize,ImageSize,3)+0.15;
        BeadImage(:,:,2)=BeadImage(:,:,2)+PositiveBeadImage;
        BeadImage(:,:,1)=BeadImage(:,:,1)+0.66*NegBeadImage;    
    else
        BeadImage=PositiveBeadImage;
    end
    imshow(BeadImage)
    end