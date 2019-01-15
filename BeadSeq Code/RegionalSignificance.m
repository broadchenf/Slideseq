function [GenesPassingSignificance,GeneNamesRestricted,pvals,MeanRegion1Fraction,MeanRegion2Fraction]=RegionalSignificance(PuckName,varargin)
    %This function can operate in two modes: asking is the gene enriched
    %relative to the puck as a whole, and asknig is it enriched relative to some
    %other region

    PuckDirectory=GetPuckDirectory(PuckName);

    BeadCountThreshold=4;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadCountThreshold"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadCountThreshold=varargin{index+1};
    end    
    FilterByCutoff=0;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="FilterByCutoff"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        FilterByCutoff=varargin{index+1};
    end
    ClusterToAnalyze=0;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Cluster"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        ClusterToAnalyze=varargin{index+1};
    end
    BeadMappingFile=FindMostRecentMapping(PuckDirectory);    
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadMappingFile"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadMappingFile=varargin{index+1};
    end
    Mode=1;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Mode"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Mode=varargin{index+1};
    end
    ImageSize=[6030,6030];
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="ImageSize"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        ImageSize=varargin{index+1};
    end
    Significance=0.001;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Significance"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Significance=varargin{index+1};
    end



    
    if length(ClusterToAnalyze)==1 && ClusterToAnalyze==0
        [GeneNames,Beads,DGE,Barcodes]=LoadBijectiveMapping(PuckName);
    elseif FilterByCutoff==0
        [DGE,Beads,Barcodes,GeneNames]=DGEByCluster(PuckDirectory,ClusterToAnalyze,'BeadMappingFile',BeadMappingFile);
    elseif FilterByCutoff>0
        [DGE,Beads,Barcodes,GeneNames]=DGEByCutoff(PuckName,ClusterToAnalyze,FilterByCutoff,'BeadMappingFile',BeadMappingFile);        
    end

    Malat1Image=PlotGeneFromName({'Malat1'},GeneNames,DGE,Beads,'PlotStyle','Default');
    
    if Mode==1
        %For Mode 1, the matrix is 
        %       Gene ~Gene
        %Region
        %~Region
        %So we are testing whether Gene is enriched in Region compared to
        %not region.
        
        disp('Draw mask')
        imagesc(Malat1Image)
        h = imfreehand; %draw something 
        mymask = h.createMask();        
        Locs=[Beads.Locations];
        BeadsInRegion=mymask(sub2ind(ImageSize,round(Locs(2,:)),round(Locs(1,:))))==1;
        
        TotalCounts=sum(DGE,2);
        GeneNamesRestricted=GeneNames(TotalCounts>BeadCountThreshold);
        
        CountsInRegion=sum(DGE(TotalCounts>BeadCountThreshold,BeadsInRegion),2);
        FractionInRegion=CountsInRegion./TotalCounts(TotalCounts>BeadCountThreshold);
        MeanRegionFraction=mean(FractionInRegion);
        CountsInRegionComplement=sum(DGE(TotalCounts>BeadCountThreshold,~BeadsInRegion),2);
        FractionInRegionComplement=CountsInRegionComplement./TotalCounts(TotalCounts>BeadCountThreshold);
        MeanRegionComplementFraction=mean(FractionInRegionComplement);

    
        for gene=1:length(GeneNamesRestricted)
            x=[CountsInRegion(gene) sum(CountsInRegion)-CountsInRegion(gene);CountsInRegionComplement(gene) sum(CountsInRegionComplement)-CountsInRegionComplement(gene)];
            [h,pval]=fishertest(x);
            hvals(gene)=h;
            pvals(gene)=pval;                
        end
        
        MeanRegion1Fraction=MeanRegionFraction;
        MeanRegion2Fraction="None";
    
        GenesPassingSignificance=GeneNamesRestricted(pvals<Significance);
        
    end    

    if Mode==2
        %For Mode 2, the matrix is
        %        Gene ~Gene
        %Region A
        %Region B
        
        disp('Draw mask for Region 1')
        imagesc(Malat1Image)
        h = imfreehand; %draw something 
        mymask = h.createMask();        
        Locs=[Beads.Locations];
        BeadsInRegion1=mymask(sub2ind(ImageSize,round(Locs(2,:)),round(Locs(1,:))))==1;

        disp('Draw mask for Region 2')
        imagesc(Malat1Image)
        h = imfreehand; %draw something 
        mymask = h.createMask();        
        Locs=[Beads.Locations];
        BeadsInRegion2=mymask(sub2ind(ImageSize,round(Locs(2,:)),round(Locs(1,:))))==1;
        
        
        TotalCounts=sum(DGE,2);
        GeneNamesRestricted=GeneNames(TotalCounts>BeadCountThreshold);
        
        CountsInRegion1=sum(DGE(TotalCounts>BeadCountThreshold,BeadsInRegion1),2);
        FractionInRegion1=CountsInRegion1./TotalCounts(TotalCounts>BeadCountThreshold);
        MeanRegion1Fraction=mean(FractionInRegion1);
        CountsInRegion1Complement=sum(DGE(TotalCounts>BeadCountThreshold,~BeadsInRegion1),2);
        FractionInRegion1Complement=CountsInRegion1Complement./TotalCounts(TotalCounts>BeadCountThreshold);
        MeanRegion1ComplementFraction=mean(FractionInRegion1Complement);

        CountsInRegion2=sum(DGE(TotalCounts>BeadCountThreshold,BeadsInRegion2),2);
        FractionInRegion2=CountsInRegion2./TotalCounts(TotalCounts>BeadCountThreshold);
        MeanRegion2Fraction=mean(FractionInRegion2);
        CountsInRegion2Complement=sum(DGE(TotalCounts>BeadCountThreshold,~BeadsInRegion2),2);
        FractionInRegion2Complement=CountsInRegion2Complement./TotalCounts(TotalCounts>BeadCountThreshold);
        MeanRegion2ComplementFraction=mean(FractionInRegion2Complement);

    
        for gene=1:length(GeneNamesRestricted)
            x=[CountsInRegion1(gene) sum(CountsInRegion1)-CountsInRegion1(gene);CountsInRegion2(gene) sum(CountsInRegion2)-CountsInRegion2(gene)];
            [h,pval]=fishertest(x);
            hvals(gene)=h;
            pvals(gene)=pval;                
        end
            
        GenesPassingSignificance=GeneNamesRestricted(pvals<Significance);
        
    end