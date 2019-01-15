function =NMFClustering(DGE,VariableGenes)

NewDGE=zeros(size(DGE));

for p=1:size(DGE,2)
    NewDGE(:,p)=log(DGE(:,p)/sum(DGE(:,p))*10000+1);
end
for q=1:size(DGE,1)
    NewDGE(q,:)=NewDGE(q,:)/std(NewDGE(q,:));
end

%Need a list of variable genes. According to the Seurat package, the 
%FindVariableGenes function identifies genes that are outliers on a 'mean variability plot'. First, uses
%a function to calculate average expression (mean.function) and dispersion (dispersion.function)
%for each gene. Next, divides genes into num.bin (deafult 20) bins based on
%their average expression, and calculates z-scores for dispersion within
%each bin. The purpose of this is to identify variable genes while controlling for
%the strong relationship between variability and average expression.

%We calculate the mean of the log values in non-log space.
LogMean=log(mean(exp(DGE),2));
%We then calculate the 
LogVMR=log(var(exp(DGE,2))./mean(exp(DGE),2));
%do NMF
nnmf

Y=tsne(NewDGE);