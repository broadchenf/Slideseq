function BeadImage=PlotCellFactorLoadings(PuckName,Cluster,mode,varargin)
%This is a hack on top of PlotGene. We get the factors and then pass them
%to PlotGene

FactorLoadings=GetFactorLoadings(PuckName);

if mode==1
    ThresholdedFactorLoadings=FactorLoadings(:,Cluster);
elseif mode==2
    ThresholdedFactorLoadings=FactorLoadings(:,Cluster).*sum(DGE,1);
end

%ThresholdedFactorLoadings=FactorLoadings(:,Cluster);
%ThresholdedFactorLoadings(ThresholdedFactorLoadings<0.05)=0;
%ThresholdedFactorLoadings(ThresholdedFactorLoadings>1)=1;


BeadImage=PlotGeneFromName('Factors',{'Factors'},ThresholdedFactorLoadings',Beads,varargin{:});