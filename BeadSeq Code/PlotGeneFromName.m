function BeadImage=PlotGeneFromName(Name,GeneNames,UniqueMappedDGE,UniqueMappedBeads,varargin)
%Varargin: option 'Overlay', if true (default false), will plot all beads
%in 
%This can take a cell of gene names, i.e., {'Prox1','C1ql2',...}

try
    GeneNum=find(sum(GeneNames==Name,2));
catch
    [~,GeneNum]=intersect(GeneNames,Name);
end
if isempty(GeneNum)
    disp('No gene by that name.')
end

if ~isempty(varargin)
    BeadImage=PlotGene(GeneNum,UniqueMappedDGE,UniqueMappedBeads,varargin{:});
else
    BeadImage=PlotGene(GeneNum,UniqueMappedDGE,UniqueMappedBeads);
end