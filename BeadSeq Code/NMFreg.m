function NMFreg(PuckName,TissueType,TissuePath,varargin)
    %PuckName is either the name of a puck or an input directory. If it is
    %an input BeadMapping directory, output directory must also be specified.
    UMICutoff=5;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="UMICutoff"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        UMICutoff=varargin{index+1};
    end
    AtlasFactors=30;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="AtlasFactors"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        AtlasFactors=varargin{index+1};
    end
    
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="OutputDirectory"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        OutputDirectory=varargin{index+1};
    else
        try
            OutputDirectory=fullfile(GetPuckDirectory(PuckName),FindMostRecentMapping(GetPuckDirectory(PuckName)));
        catch
            OutputDirectory=PuckName; %If PuckName not valid, we assume it is a dircetory, in which case OutputDirectory defaults to the same directory.
        end
    end

    try
        [GeneNames,UniqueMappedBeads,UniqueMappedDGE,UniqueMappedIlluminaBarcodes]=LoadBijectiveMapping(PuckName);
    catch
        disp('PuckName was not a valid puck name. Assuming it is a BeadMapping directory.')
        load(fullfile(PuckName,'BijectiveMapping.mat'));
    end


TissueDataPath=fullfile(TissuePath,TissueType);

%NOTE: for the matlab version, you have to add an extra comma
AtlasDGE=readtable(fullfile(TissueDataPath,'dge_hvgs_matlab.csv'));
AtlasNames=AtlasDGE.Properties.VariableNames;
AtlasNames=AtlasNames(2:end);
for p=1:length(AtlasNames) %cut out the x on the front of the rik genes
    tmp=char(AtlasNames{p});
    if tmp(1)=='x'
        AtlasNames{p}=tmp(2:end);
    end
end
    
AtlasDGE=table2array(AtlasDGE(:,2:end));

[GeneIntersection,ia,ib]=intersect(GeneNames,AtlasNames);
AtlasDGEReordered=AtlasDGE(:,ib)';
SSDGEReordered=UniqueMappedDGE(ia,:);

%Now subset the slideseq beads so they all have more than the specified
%number of reads:
SSDGEsubset=SSDGEReordered(:,sum(SSDGEReordered,1)>=UMICutoff);

%Now remove genes with no reads:
SSGenesToKeep=sum(SSDGEsubset,2)>0;
FinalGeneNames=GeneIntersection(SSGenesToKeep);
SSDGEfinal=SSDGEsubset(SSGenesToKeep,:);
AtlasDGEfinal=AtlasDGEReordered(SSGenesToKeep,:);

%Normalize each dataset so that cells sum to 1:

SSDGEfinal=SSDGEfinal./sum(SSDGEfinal,1);
AtlasDGEfinal=AtlasDGEfinal./sum(AtlasDGEfinal,1);

%we then scale the genes so that the genes have unit variance.
SSDGEfinal=SSDGEfinal./std(SSDGEfinal,[],2);
AtlasDGEfinal=AtlasDGEfinal./std(AtlasDGEfinal,[],2);

[Wa,Ha]=nnmf(AtlasDGEfinal,AtlasFactors);
HaNorm=Ha./std(Ha,[],2); %make the columns of Ha have unit variance
[~,MaxFactorIndex]=max(HaNorm,[],1);

%Map the max factors in NMF space onto the atlas clusters
AtlasClusterAssignments=readtable(fullfile(TissueDataPath,'cell_cluster_outcome.csv'));
AtlasClusters=table2array(AtlasClusterAssignments(:,2));
FactorHistogram=zeros(max(AtlasClusters),AtlasFactors);
for p=1:size(HaNorm,2)
    FactorHistogram(AtlasClusters(p),MaxFactorIndex(p))=FactorHistogram(AtlasClusters(p),MaxFactorIndex(p))+1;
end
FactorHistogram=FactorHistogram./sum(FactorHistogram,2);
[~,FactorToCelltypeMapping]=max(FactorHistogram,[],1);

%now we do the NNLS
Hs=zeros(AtlasFactors,size(SSDGEfinal,2));
for p=1:size(SSDGEfinal,2)
    Hs(:,p)=lsqnonneg(Wa,SSDGEfinal(:,p));
end
HsNorm=Hs./std(Hs,[],2); %make the columns of Hs have unit variance
[~,MaxBeadFactorIndex]=max(HsNorm,[],1);
AtlasClusters=FactorToCelltypeMapping(MaxBeadFactorIndex);
Barcodes=UniqueMappedIlluminaBarcodes(sum(SSDGEReordered,1)>=UMICutoff);
OutputTable=table((0:(length(Barcodes)-1))',Barcodes',MaxBeadFactorIndex',AtlasClusters');
OutputTable.Properties.VariableNames{2}='barcode';
OutputTable.Properties.VariableNames{3}='max_factor';
OutputTable.Properties.VariableNames{4}='atlas_cluster';
writetable(OutputTable,fullfile(OutputDirectory,'AnalogizerClusterAssignmentsOriginal.csv'));

OutputTable2=table(Barcodes',AtlasClusters');
OutputTable2.Properties.VariableNames{1}='Var1';
OutputTable2.Properties.VariableNames{2}='x';
writetable(OutputTable2,fullfile(OutputDirectory,'AnalogizerClusterAssignments.csv'));


mkdir(fullfile(OutputDirectory,'Analogizer_NMFreg_output'))

HsTable=table(Barcodes');
HsTable.Properties.VariableNames{1}='barcodes';
for p=1:size(Hs,1)
    HsTable(:,p+1)=table(Hs(p,:)');
end
writetable(HsTable,fullfile(OutputDirectory,'Analogizer_NMFreg_output',['Hs',num2str(AtlasFactors),'_0_0_17.csv']));

HaTable=table();
for p=1:size(Ha,1)
    HaTable(:,p)=table(Ha(p,:)');
end
writetable(HaTable,fullfile(OutputDirectory,'Analogizer_NMFreg_output',['Ha',num2str(AtlasFactors),'_0_0_17.csv']));
WaTable=table();
for p=1:size(Wa,2)
    WaTable(:,p)=table(Wa(:,p));
end
writetable(WaTable,fullfile(OutputDirectory,'Analogizer_NMFreg_output',['Wa',num2str(AtlasFactors),'_0_0_17.csv']));

figure(1)
imagesc(FactorHistogram)
saveas(figure(1),fullfile(OutputDirectory,'Analogizer_NMFreg_output','FactorHistogram.pdf'));