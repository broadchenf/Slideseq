%clear all
%close all
InSituBarcodePath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\InputFolder-Puck3-170701\Puck 3 Ligation AnalysisOutputs-selected-Ligation9.mat';
load(InSituBarcodePath);
ExSituBarcodePath='C:\Users\Sam\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\170719 - Puck 3 Extension Oligo 2\Base5Barcodes.mat';
load(ExSituqBarcodePath);

[IdentifiedBarcodes,IA,IB]=intersect(BeadBarcodes, Base5Barcodes);
[negIdentifiedBarcodes,negIA,negIB]=intersect(BeadBarcodes, NegBarcodes);