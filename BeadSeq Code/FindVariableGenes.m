function [VariableGenes,VariableGeneIndices,VMRZScoresOutput]=FindVariableGenes(DGE,GeneNames, varargin)

    Type='LogVMR';
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Type"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Type=varargin{index+1};
    end
    
    if string(Type)=="Dropseq"    
        
        
        
    elseif string(Type)=="LogVMR"

        DGERestricted=DGE(sum(DGE,2)>0,:);
        GeneNamesRestricted=GeneNames(sum(DGE,2)>0);
        expDGE=exp(DGERestricted);
        %We calculate the mean of the log values in non-log space.
        LogMean=log(mean(expDGE,2));
        %We then calculate the 
        LogVMR=log(var(expDGE,[],2)./mean(expDGE,2));
        LogVMR(isinf(LogVMR))=100;

        [~,SortIndices]=sort(LogMean,'descend');
        SortedVMR=LogVMR(SortIndices);
        SortedGeneNames=GeneNamesRestricted(SortIndices);
        interval=ceil(length(SortedGeneNames)/20);

        for j=1:20
            thisbin=((j-1)*interval+1):min((j*interval),length(SortedGeneNames));
            VMRZScores(thisbin)=zscore(SortedVMR(thisbin));
        end

        [~,InverseSortingIndices]=sort(SortIndices,'ascend'); %Note B=A(I), so B(J)=A(I(J)), so if we find J such that I(J)=1:len(A), then B(J)=A.
        VMRZScores=VMRZScores(InverseSortingIndices);

        IndicesPassingCutoff=find(VMRZScores>1);
        VariableGenes=GeneNamesRestricted(IndicesPassingCutoff);

        VMRZScoresOutput=zeros(1,size(DGE,2));
        VMRZScoresOutput(sum(DGE,2)>0)=VMRZScores;

        [~,VariableGeneIndices]=intersect(GeneNames,VariableGenes);
    end