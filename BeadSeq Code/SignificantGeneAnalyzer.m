function [SignificantGeneCount,SignificantGeneNames]=SignificantGeneAnalyzer(PuckList,Cluster)


SignificantGeneNames=string();
SignificantGeneCount=zeros(1,length(PuckList));
for p=1:length(PuckList)
    PuckDirectory=GetPuckDirectory(PuckList{p});
    BeadMappingFile=FindMostRecentMapping(PuckDirectory);
    try
        load(fullfile(PuckDirectory,BeadMappingFile,['Cluster_',num2str(Cluster),'_FilteredByVarianceOrExpression_SignificanceOutput.mat']))
    catch
        continue
    end
    SignificantGenes=find([Output.p]>=0 & [Output.p]<=0.005);
    for j=1:length(SignificantGenes)
        index=find(SignificantGeneNames==Output(SignificantGenes(j)).Name);
        if isempty(index)
            if SignificantGeneNames(1)==""
                index=1;
                SignificantGeneNames(1)=Output(SignificantGenes(j)).Name;
            else
                index=length(SignificantGeneNames)+1;
                SignificantGeneNames(index)=Output(SignificantGenes(j)).Name;
            end
        end
        SignificantGeneCount(index,p)=1;
    end
end