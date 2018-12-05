function Output=PermutationTestByCluster(PuckDirectory,ClusterToAnalyze,varargin)
%This takes about 20gb of ram, but you need about 75gb to calculate the
%pairwise distance matrix, after which you clear the other matrices.
%ClusterToAnalyze is a dropseq cluster, assuming the data has already been
%clustered using the AnalogizerScript.R function. If ClusterToAnalyze==0
%Variable arguments are:
%'PlotGenes', true or false. If true, will output a pdf with the
%significant genes at the 0.005 level plotted
%'NumSamples', number of samples for the null model. 1000 by default.
%Runtime scales linearly with numsamples.
%'BeadCutoff', minimum number of beads needed to assess significance of a
%gene. 15 by default.
%'BeadMappingFile', use the most recent bead mapping file by default, but
%you can supply a different one to use
%'FilterGenes', possible values are 0, in which case all genes are analyzed
%(not recommended, due to false positives); 1, in which case genes are
%filtered by within-dropseq-cluster expression; 2, in which case it's
%filtered by within-dropseq-cluster variance; or 3, in which case genes
%either match the expression cutoff or variance cutoff, and are labeled
%according to which they pass (or both).
%NOTE: In the output, passingnumber is an empty cell unless FilterGenes==3. If
%FilterGenes==3, then passingnumber is "10" if the gene passes expression
%but not variance, "01" if the gene passes variance but not expression, and
%"11" if it passes both.
%'EnforceReadNumbers' will enforce that all random samples have the same number of positive as the test sample at the cost of some computational
%time.
    EnforceReadNumbers=1;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="EnforceReadNumbers"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        EnforceReadNumbers=varargin{index+1};
    end
    FilterGenes=3;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="FilterGenes"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        FilterGenes=varargin{index+1};
    end
    PlotGenes=1;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="PlotGenes"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        PlotGenes=varargin{index+1};
    end
    NumSamples=1000;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="NumSamples"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        NumSamples=varargin{index+1};
    end
    BeadCutoff=15;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadCutoff"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadCutoff=varargin{index+1};
    end
    BeadMappingFile=FindMostRecentMapping(PuckDirectory);    
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="BeadMappingFile"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        BeadMappingFile=varargin{index+1};
    end

    if ClusterToAnalyze==0
        load(fullfile(PuckDirectory,BeadMappingFile,'BijectiveMapping.mat'))
    elseif exist(fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_UniqueMappedBeads.mat']))
        load(fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_UniqueMappedBeads.mat']));
        UniqueMappedDGE=ClusterUniqueMappedDGE;
        UniqueMappedBeads=ClusterUniqueMappedBeads;
        UniqueMappedIlluminaBarcodes=ClusterUniqueMappedIlluminaBarcodes;
        clear ClusterUniqueMappedDGE ClusterUniqueMappedBeads ClusterUniqueMappedIlluminaBarcodes
    elseif ~exist(fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_UniqueMappedBeads.mat']))
        [UniqueMappedDGE,UniqueMappedBeads,UniqueMappedIlluminaBarcodes,GeneNames]=DGEByCluster(PuckDirectory,ClusterToAnalyze,BeadMappingFile);
    end
    SignificanceOutputFile=fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_SignificanceOutput.mat']);
    if FilterGenes==1
        SignificanceOutputFile=fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_FilteredByExpression_SignificanceOutput.mat']);    
    elseif FilterGenes==2
        SignificanceOutputFile=fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_FilteredByVariance_SignificanceOutput.mat']);            
    elseif FilterGenes==3
        SignificanceOutputFile=fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(ClusterToAnalyze),'_FilteredByVarianceOrExpression_SignificanceOutput.mat']);            
    end
    
       
BeadCoords=[UniqueMappedBeads.Locations];
BeadXCoordMatrix=BeadCoords(1,:)'*ones(1,length(UniqueMappedBeads));
BeadYCoordMatrix=BeadCoords(2,:)'*ones(1,length(UniqueMappedBeads));

BeadPairwiseXValDifferences=BeadXCoordMatrix-BeadXCoordMatrix';
BeadPairwiseYValDifferences=BeadYCoordMatrix-BeadYCoordMatrix';

BeadPairwiseDistanceMat=sqrt(BeadPairwiseXValDifferences.*BeadPairwiseXValDifferences+BeadPairwiseYValDifferences.*BeadPairwiseYValDifferences);
clear BeadXCoordMatrix BeadYCoordMatrix BeadPairwiseXValDifferences BeadPairwiseYValDifferences
%To get the binedges, we do histcounts on the full matrix:
[vals,BinEdges]=histcounts(nonzeros(triu(BeadPairwiseDistanceMat)),100 );

NumReadsPerBead=sum(UniqueMappedDGE,1);
ProbabilityPerBead=NumReadsPerBead/sum(NumReadsPerBead); %The probability per bead is scaled by the number of reads per bead.

DGEBinary=UniqueMappedDGE>0;
%results=cell(1,size(UniqueMappedDGE,1));
pvals=zeros(1,size(UniqueMappedDGE,1));
effectsize=zeros(1,size(UniqueMappedDGE,1));
%delete(gcp('nocreate'));
%parpool(20)
passingnumber=cell(1,size(UniqueMappedDGE,1));

CountsPerGene=sum(DGEBinary,2);

ClusterCell={'C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','C29','C30','C31','C32','C33','C34','C35','C36','C37','C38','C39','C40','C41','C42','C43','C44','C45','C46','C47','C48','C49','C50'};
%Should be able to do both at the same time: just record in the file
%whether it passed filter by expression, variance, or both
if FilterGenes==1 || FilterGenes==3
    opts=detectImportOptions(fullfile(PuckDirectory,BeadMappingFile,'AnalogizerExpressionByDropseqCluster.csv'));
    opts=setvartype(opts,opts.VariableNames(2:end),'double');
    ExpressionTable=readtable(fullfile(PuckDirectory,BeadMappingFile,'AnalogizerExpressionByDropseqCluster.csv'),opts,'ReadVariableNames',true);
    ClusterNums={ClusterCell{1:size(ExpressionTable,2)}};
    ClusterNums{1}='GeneNames';
    ExpressionTable.Properties.VariableNames=ClusterNums;
    GoodGenes=table2array(ExpressionTable(:,ClusterToAnalyze+1))>1e-1;
    GoodGeneList=ExpressionTable.GeneNames(GoodGenes);
    if FilterGenes==3
        ExpressionGenes=GoodGeneList;
    end
end
if FilterGenes==2 || FilterGenes==3
    opts=detectImportOptions(fullfile(PuckDirectory,BeadMappingFile,'AnalogizerExpressionByDropseqCluster.csv'));
    opts=setvartype(opts,opts.VariableNames(2:end),'double');
    ExpressionTable=readtable(fullfile(PuckDirectory,BeadMappingFile,'AnalogizerExpressionByDropseqCluster.csv'),opts,'ReadVariableNames',true);    
    opts=detectImportOptions(fullfile(PuckDirectory,BeadMappingFile,'AnalogizerVarianceByDropseqCluster.csv'));
    opts=setvartype(opts,opts.VariableNames(2:end),'double');
    VarianceTable=readtable(fullfile(PuckDirectory,BeadMappingFile,'AnalogizerVarianceByDropseqCluster.csv'),opts,'ReadVariableNames',true);
    ClusterNums={ClusterCell{1:size(VarianceTable,2)}};
    ClusterNums{1}='GeneNames';
    VarianceTable.Properties.VariableNames=ClusterNums;
    GoodBeads=find(table2array(VarianceTable(:,ClusterToAnalyze+1))>0 & table2array(ExpressionTable(:,ClusterToAnalyze+1))>0); 
    GoodGenes=table2array(VarianceTable(GoodBeads,ClusterToAnalyze+1))./(table2array(ExpressionTable(GoodBeads,ClusterToAnalyze+1)).^2)>7.5 & table2array(ExpressionTable(GoodBeads,ClusterToAnalyze+1))>1e-2;
    GoodGeneList=VarianceTable.GeneNames(GoodBeads(GoodGenes));
    if FilterGenes==3
        VarianceGenes=GoodGeneList;
    end
%For debugging this, use:
    %scatter(table2array(ExpressionTable(:,ClusterToAnalyze+1)),table2array(VarianceTable(:,ClusterToAnalyze+1)))
    %set(gca,'xscale','log')
    %set(gca,'yscale','log')
    %hold on
%    plot(750*(0.1:.1:1e4).^2)
%we want to take the guys where variance > 750*expression^2
end




%% Do the test:
for geneval=1:size(UniqueMappedDGE,1)
    a=tic();
    if FilterGenes==3
        PassingVariance=0;
        PassingExpression=0;
        if any(GeneNames(geneval)==ExpressionGenes)
            PassingExpression=1;
        end
        if any(GeneNames(geneval)==VarianceGenes)
            PassingVariance=1;
        end
        if ~PassingVariance && ~PassingExpression
            pvals(geneval)=-1;
            continue
        end
        passingnumber{geneval}=[num2str(PassingExpression),num2str(PassingVariance)];
    elseif (FilterGenes>0 && ~any(GeneNames(geneval)==GoodGeneList))
        pvals(geneval)=-1;
        continue
    end

    NonzeroBeads=DGEBinary(geneval,:);
    NumBeads=nnz(NonzeroBeads);
    if NumBeads<BeadCutoff
        pvals(geneval)=-1;
%        disp(['Gene ',num2str(geneval),' completed in time ',num2str(toc(a))]);
        continue
    end
        
    %This is the true distribution
    PairwiseDistances=nonzeros(triu(BeadPairwiseDistanceMat(NonzeroBeads,NonzeroBeads))); %Triu takes the upper triangular part
    DistanceDist=histcounts(PairwiseDistances,BinEdges);
    DistSum=sum(DistanceDist);

    DistanceDist=DistanceDist/DistSum;

    %Now we generate a bunch of permuted distributions. There is a ton of duplication
    %here, because this calculation is the same regardless of geneval. It
    %only depends on the NUMBER of beads in which geneval appears.
    AverageDistribution=zeros(1,length(BinEdges)-1);
    RandomDists=cell(1,NumSamples);
    for p=1:NumSamples
        if EnforceReadNumbers           
            NonzeroBeadsRandom=datasample(1:length(NumReadsPerBead),NumBeads,'Weights',ProbabilityPerBead,'Replace',false);
            
            %NonzeroBeadsRandom=rand(3*NumBeads,length(NumReadsPerBead))<ProbabilityPerBead;            
            %num=find(cumsum(sum(NonzeroBeadsRandom,2)),1)>=NumBeads;
            %NonzeroBeadsRandom
%The following is too slow:
%            while true
%                NonzeroBeadsRandom=rand(100,length(NumReadsPerBead))/NumBeads<ProbabilityPerBead;            
%                if any(sum(NonzeroBeadsRandom,2)==NumBeads)
%                    NonzeroBeadsRandom=NonzeroBeadsRandom(find(sum(NonzeroBeadsRandom,2)==NumBeads,1),:);
%                    break
%                end
%            end                            
        else
            NonzeroBeadsRandom=rand(1,length(NumReadsPerBead))/NumBeads<ProbabilityPerBead;            
        end
%        NonzeroBeadsRandom=DGEBinary(geneval,randperm(size(DGEBinary,2))); %random integers. It would be faster to use randi probably and just generate random integers, but then we could get duplicates
        RandomDistTmp=nonzeros(triu(BeadPairwiseDistanceMat(NonzeroBeadsRandom,NonzeroBeadsRandom))); %Triu takes the upper triangular part
        RandomDists{p}=histcounts(RandomDistTmp,BinEdges); %There is maybe a way to avoid running histcounts so many times
        AverageDistribution=AverageDistribution+RandomDists{p};
        RandomDists{p}=RandomDists{p}/sum(RandomDists{p});
    end
    AverageDistribution=AverageDistribution/sum(AverageDistribution);
    %Now we calculate the distribution of norms between the random dists
    %and the average dist
    NormDist=zeros(1,NumSamples);
    for p=1:NumSamples
        NormDist(p)=sum(abs(RandomDists{p}-AverageDistribution));%This is the L1 norm, because it's fast
    end
    %Now we calculate the norm for the real distribution
    RealNorm=sum(abs(DistanceDist-AverageDistribution));
    %and a p val:
    pvals(geneval)=1-nnz(RealNorm>NormDist)/length(NormDist);
    effectsize(geneval)=RealNorm;
%    results{geneval}=[pval, RealNorm];
    disp(['Gene ',num2str(geneval),' completed in time ',num2str(toc(a)),' with p value ',num2str(pvals(geneval))]);

    if geneval/1000==floor(geneval/1000)
        Output=struct('Name',mat2cell(GeneNames,ones(length(CountsPerGene),1)),'p',mat2cell(pvals',ones(1,length(CountsPerGene))),'EffectSize',mat2cell(effectsize',ones(1,length(CountsPerGene))),'Counts',mat2cell(CountsPerGene,ones(length(CountsPerGene),1)),'PassingValue',passingnumber');
        save(SignificanceOutputFile,'Output')
    end
end

Output=struct('Name',mat2cell(GeneNames,ones(length(CountsPerGene),1)),'p',mat2cell(pvals',ones(1,length(CountsPerGene))),'EffectSize',mat2cell(effectsize',ones(1,length(CountsPerGene))),'Counts',mat2cell(CountsPerGene,ones(length(CountsPerGene),1)),'PassingValue',passingnumber');
save(SignificanceOutputFile,'Output')

if PlotGenes
    PDFOutputName=replace(SignificanceOutputFile,".mat","_significantgenes.pdf");
    SignificantGenes=find([Output.p]>=0 & [Output.p]<=0.005);
    for i=1:length(SignificantGenes)
        PlotGene(SignificantGenes(i),UniqueMappedDGE,UniqueMappedBeads,'Overlay',true);
        figure(25) %by default PlotGene uses figure(25).
        title(['Beads in Cluster ',num2str(ClusterToAnalyze),' - ',char(GeneNames(SignificantGenes(i))),' - p=',num2str(Output(SignificantGenes(i)).p),' - passingnum=',num2str(Output(SignificantGenes(i)).PassingValue)])
        export_fig(PDFOutputName,'-append')        
    end
    
    
    
end



%This is slow
%diagmat=diag(UniqueMappedDGE(geneval,:)>0);
%test4=diagmat*BeadPairwiseDistanceMat*diagmat;
%PairwiseDistanceDist=test4(find(test4>0));

%PairwiseDistanceDists{1}=diagmat*BeadPairwiseDistanceMat*diagmat