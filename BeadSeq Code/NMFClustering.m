function [GroupAssignmentsOutput,H,Y,GoodCells]=NMFClustering(DGE,GeneNames,varargin)

    k=30;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="k"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        k=varargin{index+1};
    end    

    Cutoff=10;
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="Cutoff"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        Cutoff=varargin{index+1};
    end    

    VariableGeneList="Default";
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="VariableGeneList"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        VariableGeneList=varargin{index+1};
    end  
    
    VarGenesType="Dropseq";
    index = find(cellfun(@(x) (all(ischar(x)) || isstring(x))&&(string(x)=="VarGenesType"), varargin, 'UniformOutput', 1));
    if ~isempty(index)
        VarGenesType=varargin{index+1};
    end    
    
    if VariableGeneList=="Default"
        [VariableGenes,VariableGeneIndices]=FindVariableGenes(DGE,GeneNames,'Type',VarGenesType);
    else
         vargenelist=readtable(VariableGeneList,'ReadVariableNames',0);
         vargenelist=vargenelist.Var1;
         [VariableGenes,VariableGeneIndices]=intersect(GeneNames,vargenelist);
    end
DGE=DGE(VariableGeneIndices,:);

GoodCells=sum(DGE,1)>=Cutoff;

SubSampledDGE=DGE(:,sum(DGE,1)>=Cutoff);

NewDGE=zeros(size(SubSampledDGE));

for p=1:size(SubSampledDGE,2)
    NewDGE(:,p)=log(SubSampledDGE(:,p)./sum(SubSampledDGE(:,p))*10000+1);
end
for q=1:size(SubSampledDGE,1)
    NewDGE(q,:)=NewDGE(q,:)./std(NewDGE(q,:));
end
%Because you filtered he cells after subsetting to variable genes, some
%variable genes will have no counts after removing cells.
NewDGE=NewDGE(~isnan(sum(NewDGE,2)),:);

%do NMF
[W,H]=nnmf(NewDGE,k);

for k=1:size(H,1)
    H(k,:)=H(k,:)/sum(H(k,:));
end
[~,GroupAssignments]=max(H,[],1);

GroupAssignmentsOutput=zeros(1,size(DGE,2));
GroupAssignmentsOutput(GoodCells)=GroupAssignments;

figure(1)
Y=tsne(H');
gscatter(Y(:,1),Y(:,2),GroupAssignments);