Most Slide-seq code is written in Matlab, with the exception of the NMFreg algorithm and some specific analysis scripts used for individual figures. The code is designed to run on Windows, but because Matlab is self-contained, the code should run on other operating systems with minor modifications.
The following readme is intended to assist interested users in working with the analysis code developed for the Slideseq paper.

****It is not intended to be a distributable, user-friendly software package.****
We provide no warranty. Please reach out to us if you need assistance, and we will help you however we can.

We intend to release a full Slide-seq software package at some point in the coming months.

*********Directory Structure***********
Processed slide-seq data is stored in a directory entitled "Barcodes". The base calling data for each puck is stored in a folder entitled "Puck_YYMMDD_NN", where YYMMDD is the day it was made, and NN is a number assigned to that puck.

Within the Puck_YYMMDD directory, there is a BeadMapping directory, of the form BeadMapping_M-D_hhmm, where M and D are the month and date (either 1 or 2 digits), hh is an hour and mm is a minute marker. This denotes the time when the mapping between solid barcodes and illumina barcodes was run, to prevent future mappings from overwriting old mappings.
For pucks that were mapped to illumina data after December 2018, the form is BeadMapping_YYYY_M-D_hhmm

If you are looking to analyze your data with matlab, the core data for the puck, including the DGE, gene names, and bead locations, is stored in BijectiveMapping.mat within the BeadMapping folder. Bead locations are retrieved as [UniqueMappedBeads.Locations] after loading that .mat file. We strongly recommend using the library of functions we have built for the purpose. Read the Other Analysis and Examples sections. Make sure to pay attention to the "Setting Up Your Analysis" section.
If you are looking to analyze your data with R or another language, bead locations and the DGE are available as MappedDGEForR.csv and MappedLocationsForR.csv, also in the BeadMappingFile. Otherwise, you are largely on your own, although the analyses in Fig. 1C,D; S2; and 2B were done using specialized code in R, and code for that analysis will be made available.

If you are making your own pucks, read the Puck Processing section.

*********PUCK PROCESSING************
Pucks are processed using the Puck_X_Pipeline.mat script. Parameters in the preamble must be altered, as follows:

***INPUTS***
PythonPath: If you are using NMFReg, this should point to your python executable. If you are not running NMFReg with the pipeline, this variable has no effect.
BeadseqCodePath: This variable should point to your Beadseq Code directory.
PipelineFunctionPath: This variable should point to your PipelineFunctions directory.
BeadType: If using beads with "Sequence 1" in the supplement, this should be "Reverse Phase." If using beads with "Sequence 2," this should be "180402."
RenameFiles: Should be 1 by default.
RunSignificanceAnalysis: If 1, then after the Illumina reads are mapped onto beads in space and after NMFReg is run, the algorithm described in the supplement under "Significant Gene Calling" will automatically be run for each NMFreg cluster.
NumClusters: This is the number of clusters present in the scRNAseq or snRNAseq dataset being used for NMFReg.
RunAnalogizer: This is a list of N numbers, where N is the number of pucks being analyzed. If non-zero, NMFReg will be run on the corresponding puck. The integer used corresponds to a particular tissue type, with 1=Cerebellum and 2=Hippocampus. Code can be modified to allow more other tissue types. The length of RunAnalogizer should be the same as the length of PuckNames.
AnalogizerBeadCutoff: The NMFReg algorithm will exclude beads with fewer than this number of variable genes.
AnalogizerType="NMFReg" by default.
NMFregFactorNum: Number of factors to use for NMF by NMFreg. 30 by default.
NMFregPath: The path to the references for NMFreg
DropseqDGEPaths: Deprecated.
DropseqClusterPaths: Deprecated.
DropseqMeanAndVariancePath: This should point to a file that includes the mean and variance for each gene in each cluster for the dataset being used as a reference for NMFReg.
NumPar: This is the maximum number of threads available for ParPool.
CropImage: If 1, the program will prompt the user to crop the puck freehand after mapping illumina barcodes. The image used for cropping shows expression of Malat1.
CropSuffix: A suffix to be appended to the BijectiveMapping.mat, MappedDGEForR.csv, and MappedLocationsForR.csv files after cropping.
EnforceBaseBalance: If set to 1, then the base caller will adjust the gain on the various channels until the number of pixels called as the most common base minus the number of pixels called as the least common base is only BaseBalanceTolerance times the total number of pixels. This is only computed on the central 1000x1000 pixels of the image, and is problematic if the beads are sparse/if background pixels are called as a base.
SaveData: 1 by default.
IlluminaReadThreshold and MaximumBarcodesToAnalyze: For mapping, we consider at most MaximumBarcodesToAnalyze beads, requiring that each bead have at least IlluminaReadThreshold trancripts.
NumLigations: 20 by default
NumBases: 14 by default
BarcodeSequence: Indicates which ligations should be used in the barcode.
FolderWithRawTiffs: Points to the folder that the raw tiffs are stored in.
FolderWithProcessedTiffs: Points to the folder that the renamed and registered tiffs are stored in.

IndexFiles: The pipeline is designed assuming that some the ligations are imaged into nd2 files, and are then exported using nis elements to TIFFs. After exporting, they have the format DescriptiveNametTTxyMM, where TT refers to the timepoint within the nd2, and MM is the index of the ROI corresponding to that puck in the nd2. IndexFiles is a list of the "DescriptiveName"s for each of the files.
LigationToIndexFileMapping: For each ligation, this is the IndexFile within which the ligation can be found.
tnumMapping: For each ligation, this is the timepoint corresponding to that ligation within its IndexFile.
PucksToAnalyze: This indicates which ROIs in the IndexFiles correspond to pucks. All IndexFiles should have the same number of ROIs.
PuckNames: Descriptive names given to the pucks, used throughout when creating the file structure. The ith entry of PuckNames is the name that will be assigned to the ROI at position PucksToAnalyze(i). By default, the puck name should be of the form Puck_YYMMDD_NN, where NN is an index.
OutputFolderRoot: This is the "Barcodes" folder, to which base calls will be outputted.
IlluminaFolderRoot: This is the folder within which one can expect to find DGEs. The DGE must be in that folder, and must contain "YYMMDD_NN." or "YYMMDD_NN_" If there is not a unique DGE in the IlluminaFolderRoot directory with one of those two constructs in its name, the puck will be skipped during alignment of Illumina and Solid barcodes.
ImageSize: This is the size of the images. It should be slightly smaller than the actual size. During image registration, the program will crop the images to take pixels 1:ImageSize by 1:ImageSize, to ensure all the images are the same size. (Sometimes they might not be, in cases of stitching errors.)
XCorrBounds: This gives the corners of the box used for cross-correlation channel coregistration, in the form x1, x2, y1, y2.
RegisterColorChannels: Register the colors of each image to each other before image registration.
BeadZeroThreshold: During basecalling, only accept barcodes if they have at most BeadZeroThreshold "zeros" (i.e., if they failed to be called as a base on at most BeadZeroThreshold bases).
PixelCutoffRegistration: Pixel cutoff to use during color channel coregistration. Should be set manually.
PixelCutoffBasecalling: Pixel cutoff to use during basecalling. Should be set manually. Probably should be the same or slightly lower than PixelCutoffRegistration.
DropBases: If set to 1, then during basecalling, the program will count how many barcodes it has passing filter when it uses all 14 meaningful ligations, as well as how many barcodes it has passing filter when it uses only 13 ligations, for each of the 14 subsets of 13 ligations. It will then drop one ligation if dropping that ligation increases the number of barcodes passing filter by at least 4000 barcodes. The point here is that if some bases are very low quality, or failed image restration, they will dramatically reduce the number of barcodes passing filter for that puck.
BeadSizeCutoff: This is the minimum number of pixels belonging to a connected component in Beadseqfun, in order for that connected component to be called as a barcode. For most of the duration of the project, it was set to 50. If set to 50, only high-quality beads pass the filter, and the number of barcodes passing filter following base calling is then a good predictor of the quality of the puck (for a 3mm puck, 50k-60k is acceptable, 60k-70k is good, and 70k+ is outstanding). However, if image registration is poor, it tends to decrease the size of the connected components, and one can lose beads as a result. We found that setting BeadSizeCutoff to 30 makes the basecalling more robust against minor issues with image registration, but then the number of barcodes passing filter no longer appears to be a good predictor of the quality of the puck. To understand this, one should look in particular at the first page of the "Report_BaseCalling" pdf. In particular, if the image registration was good, there should be a clear hump visible in the histogram on the first page of "Report_Bascalling" (see e.g. Report_BaseCalling_Puck_180819_9), with a peak around 100-200 and a trough around 50. If the image registration was bad, leading to division of beads, that hump will disappear. In that case, having a lower BeadSizeCutoff can greatly assist in recovering beads.

PrimerNLigationSequence, PrimerUPLigationSequence: This is the order in which the ligations are performed.
InverseLigationSequence: This is the inverse of PrimerNLigationSequence, and is used when writing out the color space barcodes in human readable format.
WhichLigationsAreMissing: In case a ligation is missing, this variable can be used to get the program to ignore the fact that that ligation is missing. However, in that case, it's usually easier just to substitute an existing ligation for the missing ligation. For example, if Truseq ligation 2 is missing, simply using Truseq Ligation 1 as a stand-in for Truseq Ligation 2 (using the tnumMapping variable) is very effective.

NOTE: Some parameters are also hard-coded into the functions. For example, there are several parameters hard-coded into Beadseqfun6, in particular regarding the treatment of phasing, which may need to be altered for different microscopy setups.

***OUTPUTS***
The pipeline outputs a large number of files. The most important are as follows:

Basecalling outputs (outputted to Barcodes/PuckName):
1) BeadImage: All pixels that belong to a connected component with a barcode that passed filter are plotted in this image.
2) ChannelCalls1 through 4: For each channel, ChannelCallsN shows all pixels that were called as the Nth channel.
3) Report_Basecalling: The first page shows the number of pixels having each barcode, for *all* barcodes identified in the basecalling. Note that this report is written out *before* the connected components segment of the basecalling. Thus, for example, the mode number of pixels per barcode is 1. The second page shows the base balance across all barcodes with more than BeadSizeCutoff pixels. The third image is a histogram of the total number of zeros (i.e., bases on which no channel passed the PixelCutoffBasecalling threshold) per barcode.
4) Metrics: This records the parameters used during basecalling, and the number of barcodes passing filter. Note that the number of barcodes passing filter, as recorded by metrics, is the same as the number of white clusters displayed in BeadImage.

Barcode Matching:
Upon matching SOLiD barcodes to Illumina barcodes, the pipeline generates a function within Barcodes\PuckName entitled "BeadMapping_YYYY_MM-DD_hhmm." This is intended to avoid accidentally overwriting the bead mapping data if the bead mapping is later rerun with a different DGE.
1) Metrics.txt gives a variety of useful metrics regarding the bead mapping. The most important of these is the total number of bijective pairings, which is the key metric. If you are using BeadSizeCutoff=50, then the ratio of Bijective Pairings to Barcodes passing filter is a good measure of what fraction of beads were mapped in Illumina. However, with BeadSizeCutoff=30, because some beads can give rise to two barcodes passing filter, and because the number of barcodes depends on image registration (with worse image registration generating more barcodes due to bead splitting) that ratio is no longer clearly a good metric.
2) BijectiveMapping.mat: This contains the DGE for bijectively mapped beads, the Beads structure (which contains in it for each bijectively mapped bead, the bead barcode and the location of the bead), the Barcodes cell which contains the Illumina barcode for each bead, and the GeneNames cell, which contains a list of all gene names in the DGE. Note that the jth column of the DGE is Beads(j), and the ith row of the DGE is GeneNames{i}. Note that a SOLiD barcode B and an Illumina barcode I are said to be bijectively mapped if there are no other Illumina barcodes with hamming distance equal or closer to B than I, and if there are no other SOLiD barcodes with hamming distance equal or closer to I than B.
3) Report_Puck: This is your one-stop-shop for puck quality metrics. The important pages are as follows. Page 1 shows the hamming distance distribution mapping SOLiD barcodes to Illumina barcodes, for the SOLiD barcodes found during basecalling (in blue) and for SOLiD barcodes obtained by scrambling the barcodes found by basecalling (in orange). Note that if using DropBases=1, and if a base is dropped, then one expects the blue distribution to be peaked at 1. Otherwise, one expects it to be peaked at 0. 
4) MappedDGEForR: The bijectively mapped DGE, exported as a csv for use in R.
5) MappedLocationsForR: The locations of the beads in MappedDGEForR, exported as a CSV for use in R.

NMFReg:
In order to use NMFReg, DropseqDGEPaths and DropseqClusterPaths must be set, and appropriately formatted DGEs and Cluster mappings must be present at the indicated paths. All output is within the BeadMapping directory.

1) AnalogizerClusterAssignmentsOriginal: For each unique mapped barcode, this contains the maximum factor and the cluster to which the barcode was assigned.
2) AnalogizerClusterAssignments: This is a subset of AnalogizerClusterAssignmentsOriginal.
3) Analogizer_NMFreg_output\Hs30_0_0_17: This is the H matrix for slideseq data (as opposed to Ha30_0_0_17, which is the H matrix for the atlas data). Using this Hs matrix and AnalogizerClusterAssignmentsOriginal to determine the mapping between factors and clusters, one can determine the weight associated with each cluster for each bead. That is done by the GetFactorLoading function.

Significance Analysis:
The function that implements the significant gene analysis is entitled PermutationTestByCluster. In order to use the significance analysis, DropseqMeanAndVariancePath must be set, and appropriate formatted matrices containing the means and variances for each gene in each cluster must be present at the target of that path variable. This is necessary because the significance analysis will only consider genes with high expression, high variance, or both, in the indicated cell type. All output is within the BeadMapping directory.

When run on a particular cluster, PermutationTestByCluster outputs a pdf of all genes passing the significance threshold at the p=0.005 threshold within that cluster. It also outputs the SignificanceOutput variable (e.g. "Cluster_1_FilteredByVarianceOrExpression_SignificanceOutput.mat") which contains, for each, gene, the associated p value within that cluster.

Use the SignificantGeneAnalyzer function to aggregate the output of PermutationTestByCluster across multiple pucks.

Note that as currently implemented, and as used throughout the paper, PermutationTestByCluster analyzes all beads assigned to a particular cluster by NMFreg. We recommend adjusting the code (using GetFactorLoadings) to restrict instead to the set of beads that have more than a particular factor loading in the cluster of interest. This is a more principled way of determining which beads to perform the analysis on. 

*********OTHER ANALYSIS***********
>>>>>>>Setting up your analysis:
After installing Matlab, we recommend that users edit the contents of the GetPuckDirectory function in Beadseq Code to point to their Barcodes directory. The Barcodes directory is the directory labeled "Barcodes" in the Geo repository if you are using the data from the Manuscript, or is the target of the OutputFolderRoot variable in the pipeline function if you are using your own pucks.
After editing that line in GetPuckDirectory, you will be able to conveniently access the data from the paper using LoadBijectiveMapping, as described below.

In addition, in order to be able to run NMFreg on the tissue types used in the paper, you must take the following steps:


In order to run 

As with most Matlab functions, optional arguments are provided as Name,Value pairs.

>>>>>>>>LoadBijectiveMapping
Syntax: [GeneNames,Beads,DGE,Barcodes]=LoadBijectiveMapping(PuckName)
PuckName is the descriptive name of the puck used in the pipeline, and is also the name of the folder corresponding to that puck within the Barcodes directory.

Optional inputs:
--"CropSuffix": loads the variables associated with the cropped puck with a given CropSuffix.
--"BeadMappingFile": The bead mapping file to use. By default, it will use the most recently dated file.

Outputs:
This function loads the variables in the BijectiveMapping.mat file for the most recently dated BeadMapping directory associated with the given puck.
---GeneNames: A list of genes
---Beads: A struct of data associated with each bijectively mapped bead. The most important field is Locations, which gives the centroid of the bead location. Pixels gives a list of the pixels (in linear indexing) that belonged to the bead.
---DGE: A DGE, genes by beads, where the genes in the DGE are in the same order as GeneNames and the beads are in the same order as Beads.
---Barcodes: Illumina barcodes associated with beads.

>>>>>>>>GetPuckDirectory
Syntax: directory=GetPuckDirectory(PuckName)

This function simply returns the string BarcodesDirectory\PuckName, where BarcodesDirectory is (hard-coded) a link to the directory Barcodes. Users are advised to alter this function prior to using the code package.

>>>>>>>>DGEByCutoff
Syntax: [DGE,Beads,Barcodes,GeneNames, FactorWeight]=DGEByCutoff(PuckName,Cluster,Cutoff)

DGEByCutoff subsets the DGE on the basis of an NMFReg factor loading. It keeps all beads with an NMFReg factor loading greater than a specified cutoff. Note that in this approach, a given bead can belong to more than one cluster. This is the method used in Fig. 2C.

Inputs:
PuckName is the descriptive name of the puck used in the pipeline function. It is also the name of the folder containing the data for that puck in the Barcodes directory. Alternatively, PuckName can be a puck directory, of the form GetPuckDirectory(PuckName)
Cluster is an integer or a list of integers indicating the clusters to keep.
Cutoff is the NMFReg factor loading cutoff.

Outputs:
Same as LoadBijectiveMapping, but all variables are subsetted to the beads passing the cutoff filter.
FactorWeight is a list of the factor loadings for each bead in the subsetted DGE, in each of the specified clusters. Note that the order of the columns of FactorWeight is given by the Cluster variable. Hence, if we have Cluster=[6,3,4], then the first column of FactorWeight is the factor weightings for cluster 6.

Optional inputs:
CropSuffix: Passes this argument to GetFactorLoadings
BeadMappingFile: Passes this argument to GetFactorLoadings


>>>>>>>>DGEByCluster
Syntax: [DGE,Beads,Barcodes,GeneNames]=DGEByCluster(PuckDirectory,ClusterToAnalyze)
DGEByCluster subsets the DGE by assigning each bead to the cluster with the greatest factor loading. Note that in this approach, each bead belongs to exactly one cluster. This is the method used in figures 3 and 4. However, it was used for figures 3 and 4 primarily because the DGEByCutoff method was developed later, and in general we recommend starting with the DGEByCutoff method.

Inputs:
PuckDirectory: Unlike DGEByCutoff, DGEByCluster takes a directory as input, for legacy reasons. You can supply it GetPuckDirectory(PuckName).
ClusterToAnalyze: Same as the "Cluster" variable for DGEByCutoff.

Optional inputs:
CropSuffix: Passes this argument to GetFactorLoadings
BeadMappingFile: Passes this argument to GetFactorLoadings

>>>>>>>>GetFactorLoadings
Syntax: [FactorWeight,ClusterUniqueMappedDGE,ClusterUniqueMappedBeads,ClusterUniqueMappedIlluminaBarcodes,GeneNames]=GetFactorLoadings(PuckName,varargin)

Note that DGEByCutoff is a wrapper on GetFactorLoadings.

Inputs:
--"PuckName" - the name used for the puck in the Barcodes directory.

Outputs:
FactorWeight: A matrix of the factor weights for each atlas cluster, for each bead. Technically, the i,jth entry is the L2 norm of all factors mapping to Atlas cluster j, for the ith bead.
[ClusterUniqueMappedDGE,ClusterUniqueMappedBeads,ClusterUniqueMappedIlluminaBarcodes,GeneNames]: same as the outputs for LoadBijectiveMapping, but restricted to the cluster(s) specified by the Cluster argument. These are provided as optional outputs for convenience.

Optional inputs:
--BeadMappingFile: See above
--CropSuffix: See above
--SaveOutput: 0 by default. If set to 1, the outputs will be saved to fullfile(PuckDirectory,BeadMappingFile,['FactorWeights',CropSuffix,'.mat']);
--Cluster: If one or more clusters are provided (multiple can be provided as a list), the DGE and outputs are restricted to beads that were called as one of the specified clusters by NMFReg, using the maximum factor loading method. Note that GetFactorLoadings doesn't accept "FilterByCutoff" as an argument.

Note that even if Cluster is not specified, the DGE you get with GetFactorLoadings is smaller than the DGE you get with LoadBijectiveMapping, because GetFactorLoadings only returns the beads that passed the threshold used for NMFReg (at least 5 variable transcripts by default).

>>>>>>>>NMFreg
Syntax: NMFreg(PuckName,TissueType,TissuePath,varargin)
This function runs NMFreg.

In order to run NMFreg, you must have the files dge_hvgs_matlab.csv and cell_cluster_outcome.csv in a directory (see TissueType, below). dge_hvgs_matlab.csv is the atlas DGE (beads by genes), with barcodes in the first column and genes in the first row. cell_cluster_outcome.csv is a csv with two columns, the first having barcodes and the second having atlas cluster assignments. Note that the first entry of dge_hvgs_matlab.csv must be a comma, i.e. the 0,0th position of the array must be empty for the gene names to be read correctly.

Inputs:
--PuckName: A valid PuckName. If the Puckname is not valid, then it is interpreted as a BeadMapping directory, which is then also the default OutputDirectory.
--TissueType: A character string. The directory TissuePath\TissueType must exist, and must contain the files dge_hvgs_matlab.csv and cell_cluster_outcome.csv. Acceptable tissue types at time of deposition are 'hippocampus', 'cerebellum', 'kidney' and 'liver', if TissuePath points to the NMFreg Dependencies folder deposited in Geo.
--TissuePath: See TissueType.

Optional inputs:
--UMICutoff: The minimum number of variable genes a bead must have to be included in NMFreg.
--AtlasFactors: The number of factors to use in the atlas NMF decomposition. 30 by default.
--OutputDirectory: The directory to write the output files to.

Outputs: The function has no outputs. However, it creates the following files within the OutputDirectory:
--The Analogizer_NMFreg_output directory, which contains the H matrix for the atlas decomposition (Ha_30_0_0_17), the W matrix for the atlas decomposition (Wa_30_0_0_17), and the H matrix for the slide-seq decomposition (Hs_30_0_0_17). In addition, it contains a heatmap showing the mapping from NMFreg factors to Atlas cell types.
--AnalogizerClusterAssignmentsOriginal.csv, which gives the barcode for each bead called by NMFreg, its maximum factor loading, and the corresponding atlas cluster as determined using the maximum factor loading method.
--AnalogizerClusterAssignments.csv, which just contains barcodes and atlas cluster assignment.

>>>>>>>>RegionalSignificance
Syntax: [GenesPassingSignificance,GeneNamesRestricted,pvals,MeanRegion1Fraction,MeanRegion2Fraction]=RegionalSignificance(PuckName,varargin)
This function allows the user to determine which genes have significantly different expression in two user-defined regions. The comparison is a Fisher exact test, performed on the matrix: [X1, T1-X1; X2, T2-X2], where X1 is the counts of the gene in region 1, and T1 is the total counts of all genes in region 1.

Output:
---GenesPassingSignificance is a cell of all of the names of genes passing the significance threshold (0.001 by default).
---GeneNamesRestricted is the set of all genes passing the BeadCountThreshold. GenesPassingSignificance=GeneNamesRestricted(pvals<Significance).
---pvals is the matrix of raw p values. It is the same size as GeneNamesRestricted.
---MeanRegion1Fraction: For all genes in GeneNamesRestricted, MeanRegion1Fraction is the fraction of gene expression in region 1.
---MeanRegion2Fraction: For all genes in GeneNamesRestricted, MeanRegion2Fraction is the fraction of gene expression in region 2.

Optional Inputs:
---"Mode": This is either 1 or 2. If 1 (as by default), then it prompts the user to draw a single region, and sets the second region to be the complement of the first region. Thus, it evaluates whether genes are significantly enriched or depleted relative to the rest of the puck.
---"BeadCountThreshold": Genes must have more than this many counts (4 by default) on the whole puck to be evaluated. Note that the false discovery rate is the number of genes passing BeadCountThreshold times the Significance threshold.
---"BeadMappingFile": The bead mapping file to use. By default, it will use the most recently dated file.
---"ClusterToAnalyze" and "FilterByCutoff: If ClusterToAnalyze is 0 (by default), all beads will be considered. Otherwise, if ClusterToAnalyze is greater than 0 and FilterByCutoff is 0, then analysis will be restricted to beads that belong to the cluster numbered ClusterToAnalyze by NMFReg, according to the maximum factor loading method. If ClusterToAnalyze is greater than 0 and FilterByCutoff is greater than 0, then beads will be analyzed if the factor loading corresponding to the cluster numbered ClusterToAnalyze is greater than FilterByCutoff. Note that ClusterToAnalyze can be a list, in which case all clusters in the list will be included.
---"ImageSize": The size of the puck images. [6030, 6030] by default.
---"Significance": The significance threshold for reporting GenesPassingSignificance.


>>>>>>>>GetBeadsInRegion
Syntax: BeadIndices=GetBeadsInRegion(DGE,Beads,GeneNames)
This function displays an image of all beads expressing Malat1, which is typically helpful for visualizing cell bodies. It then prompts the user to crop the image freehand.

Output: A binary vector with the same size as Beads. Beads(BeadIndices) gives all beads that are within the convex hull of the region drawn by the user.

>>>>>>>>FindVariableGenes
Syntax: [VariableGenes,VariableGeneIndices,VMRZScoresOutput]=FindVariableGenes(DGE,GeneNames)
Calculates a list of variable genes from the given gene, using the logVMR method similar to that in the Seurat package.

Inputs:
--DGE
--GeneNames

Outputs:
--VariableGenes: a list of variable genes
--VariableGeneIndices: VariableGenes=GeneNames(VariableGeneIndices)
--VMRZScoresOutput: By default, the variable genes are those with variance Z scores greater than 1 in their bin. This variable gives the raw Z scores to facilitate alternative cutoffs.

>>>>>>>>NMFClustering
Syntax: [GroupAssignments,H,Y,GoodCells]=NMFClustering(DGE,GeneNames,varargin)
This function performs NMF clustering on the DGE. Finally, it plots a tSNE of the resulting clustering.

Inputs:
DGE: The slide-seq DGE.
GeneNames: A list of GeneNames, used to intersect the variable genes with the DGE gene names.

Outputs:
--GroupAssignments: A list of NMF cluster assignments for each bead. The clustering is done by maximum factor.
--H: The H matrix from the NMF decomposition.
--Y: The coordinates of each bead in the tSNE. The tSNE can be plotted with "gscatter(Y(:,1),Y(:,2),GroupAssignments)";
--GoodCells: A logical array indicating which beads were retained for the NMF. Beads are only retained if they have more than some number of variable genes (10 by default, see Cutoff, below). 

Optional Inputs:
--"k": the number of components to be used in NMF.
--"Cutoff": The minimum number of variable genes required for a bead to be retained for NMF.
--"VariableGeneList": A path to a list of variable genes, which shoud be an N by 1 csv. Lists for the cerebellum, hippocampus, and kidney are included for convenience in the BeadseqCode repository. If VariableGenes="Default" (the default value), then the variable genes will be calculated using FindVariableGenes, which may have poor performance, and which has not been validated.

This function is intended to be used with PlotGroupsInSpace, see "Examples" below.

>>>>>>>>PlotGeneFromName:
Syntax: BeadImage=PlotGeneFromName(Name,GeneNames,DGE,Beads)
This is the most useful plotting function. Name is a name of a gene (i.e. a member of GeneNames), or a cell of gene names, for plotting. Note that if multiple gene names are supplied then the counts of the genes supplied are summed for plotting. Supplying GeneNames as the first argument (as well as the second argument) will plot the sum of all counts.

NOTE: PlotGene and PlotGeneFromName do not currently accept the FilterByCutoff argument for subsetting beads by a factor loading cutoff, unlike RegionalSignificance and PermatutionTestByCluster. Instead, users are encouraged to use DGEByCutoff to subset the DGE and Bead variables appropriately, and then supply them directly to the PlotGeneFromName function. Likewise, although one is free to use the "Cluster" argument to subset beads based on maximum NMFreg factor loading, it would be faster to subset the DGE manually using DGEByCluster if one wants to plot many different genes.

Optional inputs:
--"PlotStyle":
------"Default" - Plots all beads that have the gene in green, and all other beads in red.
------"Circles" - plots all beads having the gene as a blue circle with radius proportional to the number of counts of the gene on that puck. Plots a subset of other beads as green circles of constant radius.
------"Imagesc" - plots the counts as a heatmap.
------"Imagescfei" - Plots the counts as a heatmap, with smoothing, as described under "Density Plots" in the supplementary material. If using ImagescFei, the image will not automatically be displayed: one must subsequently run imagesc(BeadImage) to view the image.
------"ImagescClusterByLoading" - REPLACE THIS WITH FILTERBYCUTOFF
--"Overlay": If false (true by default) for "Circles" or Default", beads that do not express the gene will not be displayed.
--"OverlayDownsample": 5 by default. The number of beads without any counts of the gene plotted by "Overlay" is equal to the total number of such beads divided by OverlayDownsample.
--"Cluster": 0 by default (i.e. no subsetting of beads). Limits display to beads that were determined to belong to the numbered cluster by the maximum factor loading method in NMFReg. Note that a list may be supplied, in which case all beads belonging to clusters in the list will be included.
--"FilterByCutoff": Not currently supported.
--"PuckName": If Cluster is non-zero, then "PuckName" must be specified so the program can go and retrieve the list of beads belonging to each cluster.
--"Normalize": If "Normalize" is set to 1, then the counts of the indicated transcripts as a fraction of the total number of counts on the bead is plotted. If 0, the absolute number of transcripts is plotted. "Normalize" is currently only implemented for PlotStyle="Circles".
--"BeadSizeFactor": For "Circles", the bead radius is scaled by BeadSizeFactor.
--"FigNum": Determines which figure number to use when plotting.

Examples:
PlotGeneFromName('Mybpc1',GeneNames,DGE,Beads,'PlotStyle','Circles','Cluster',2,'PuckName','Puck_180819_12');
Note that if you are using the 'Cluster' argument for PlotGeneFromName (which only plots genes restricted to beads called as a particular 

PlotGeneFromName will also accept cells containing multiple gene names, and will plot them together. For example:
PlotGeneFromName({'Pcp2','Pcp4'},GeneNames,DGE,Beads,'PlotStyle','Imagesc');


>>>>>>>>PlotCellFactorLoadings
Syntax: BeadImage=PlotCellFactorLoadings(PuckName,Cluster,mode,varargin)

For the given Cluster, this function plots the factor loading for the given Cluster for each bead on the puck. The plotting is the same as if it were plotting a gene with a number of counts equal to the factor loading. Use of BeadSizeFactor may be required for PlotStyle="Circles".

If mode=1, it plots the raw factor loading for each bead.
If mode=2, it plots the factor loading multiplied by the number of transcripts appearing on the bead.

All variable arguments are passed directly to PlotGeneFromName.

>>>>>>>>PlotGroupsInSpace
Syntax: PlotGroupsInSpace(Beads,GroupAssignment)
This function returns an image of all the beads in Beads, colored according to a scalar in GroupAssignment. It is intended to be used with NMFClustering.

>>>>>>>>SignificanceOverlap
Syntax: [Cooccurrences,CorrelatingGenes,Anticooccurrences,AnticorrelatingGenes,GoodGeneList]=SignificanceOverlap(GeneList,PuckDirectories,varargin)

Inputs:
--GeneList: A string array of gene names. All the genes in the genome will be compared to the genes in GeneList to identify correlation or anticorrelation. Adding additional genes to GeneList increases the runtime, but not linearly, at least when GeneList is small.
----NOTE: GeneList will be alphabetized and intersected with the set of genes occurring on the pucks by the function. We refer to the resulting list as GoodGeneList below. If GeneList contains some genes that don't occur on all the pucks, it will crash.
--PuckNames: A cell of puck names to analyze.

Outputs:
--Correlations: A binary matrix of size A x length(GeneList) x length(PuckNames). A is the number of genes that correlated with at least one of the genes in GeneList on at least one puck. If the i,j,kth entry is 1, then CorrelatingGenes{i} correlated with the jth entry of GoodGeneList on the kth puck, where the two genes are said to be correlated if the overlap of the two genes is at least ZScoreCutoff z scores greater than the mean overlap between the shuffled gene images for the same genes.
--CorrelatingGenes: A list of the genes on the first axis of Correlations.
--Anticorrelations: Same as correlations, but for genes with overlaps at least ZScoreCutoff z scores below the mean.
--AnticorrelatingGenes: A list of the genes on the first axis of Anticorrelations.
--GoodGeneList: A list of the genes in GeneList that also appear in the DGEs of the labeled pucks.

Optional Inputs:
--Spread: 5 by default. Beads are represented by squares with side length equal to 2*Spread after downsampling.
--DownsampleFactor: 10 by default. Images of the beads are downsampled by this factor prior to the analysis, to accelerate things. For the purpose of calculating the overlap, beads are effectively represented by squares with side length equal to Spread*DownsampleFactor.
--ZScoreCutoff: 3 by default. After calculating Z scores, genes with Z greater than ZScoreCutoff are considered to be correlated on that puck, while those with Z less than negative ZScoreCutoff are considered to be anticorrelated on that puck. The false positive rate is extremely high, so it is necessary to have multiple pucks to eliminate false positives.
--BeadCutoff: 4 by default. We only consider genes that are represented on more than BeadCutoff beads on the puck.
--ImageSizeX and ImageSizeY: 6030 by default.
--NumSamples: 100 by default. This is the number of samples calculated in order to determine a mean and variance for calculating Z scores. The runtime is roughly linear in this value.
--ClusterToAnalyze and FilterByCutoff: 0 by default. See description of these variables under RegionalSignificance.

We do not currently output the raw matrices of Z scores, although in principle it might be better to do clustering on that matrix, rather than doing clustering on the Cooccurrences and Anticooccurrences matrix. Note that currently all images are represented as binary, so there is no accounting for cases in which a bead has two counts of the same gene, or when two adjacent beads express the same gene (ideally, you might sum the expression counts of overlapping beads, but this is not currently done).

>>>>>>>>PermutationTestByCluster
function Output=PermutationTestByCluster(PuckDirectory,ClusterToAnalyze,varargin)

Inputs:
--PuckDirectory, of the form PuckDirectory=GetPuckDirectory(PuckName)
--ClusterToAnalyze, a cluster from the NMF factorization.

Outputs:
--An Output struct, which includes the p value associated with each gene. P=-1 if the gene did not pass filter, either due to insufficient reads on the puck or due to insufficient expression or variance..
--Optionally, a PDF. See "PlotGenes" below.

Optional inputs:
--FilterByCutoff: See description under "RegionalSignificance"
--EnforceReadNumbers: If set to 1, then all random permutations generated by the function (for the purpose of generating the null model) have the same number of transcripts as are found on the puck. If set to 0, then the number of transcripts in each random sample will be a pseudobinomial distribution with mean equal to the number of transcripts found on the puck. If EnforceReadNumbers=0, the function will run ~10x faster, but will generate a non-uniform P distribution peaked at 0.5.
--FilterGenes: If set to 0, the function will run on all genes in the DGE. If set to 1, genes will be filtered according to their expression levels in dropseq. If set to 2, they will be filtered by their variance in dropseq. If set to 3, they will be filtered by either expression or variance. If FilterGenes is set to 3, then the output variable and PDFs contain a binary string indicating whether a gene passed the variance ('01') or expression ('10') filter, or both ('11'). Note that if dropseq expression and variance levels are not available, FilterGenes must be set to 0, or else there will be an error. In order for the function to run with FilterGenes>0, there must be a file titled AnalogizerExpressionByDropseqCluster.csv and/or a file called AnalogizerVarianceByDropseqCluster.csv in the beadmapping file. If RunAnalogizer is non-zero for the puck in the pipeline function, this file is automatically transferred by the pipeline function from a directory within BeadseqCode into the BeadMappingFile. Otherwise, it is incumbent upon the user to generate the appropriate file and transfer it manually. For each gene and each cluster from the Drop-seq data, the file gives the expression or variance (respectively) in that cluster in the dropseq data.
--PlotGenes: If set to 1, a PDF will be generated for each gene identified as significant at the 0.005 level. The PDF will be labeled according to the cluster number and the value of FilterGenes and will be put into the appropriate BeadMapping folder. If a very large number of genes are found to be significant, this can take a large amount of time.
--NumSamples (1000 by default): This is the number of random samples that will be generated for the purpose of generating the null model. The smallest measurable P value is 1/NumSamples. The runtime is proportional to NumSamples.
--BeadCutoff (15 by default): Genes with fewer than this many counts across the entire puck will not be analyzed.
--BeadMappingFile: See descriptions of this option above.

***************EXAMPLES****************
%Load data for puck 180618_4:
>>[GeneNames,Beads,DGE,Barcodes]=LoadBijectiveMapping('Puck_180618_4');

%Plot the counts of Pcp4 on the puck:
>>PlotGeneFromName('Pcp4',GeneNames,DGE,Beads);
%Plot as a heatmap instead:
>>PlotGeneFromName('Pcp4',GeneNames,DGE,Beads,'PlotStyle','Imagesc');

%Plot the combined counts of Pcp4 and Cux2:
>>PlotGeneFromName({'Pcp4','Cux2'},GeneNames,DGE,Beads);

%Do a tSNE and plot the spatial bead assignments using pre-existing variable genes:
>>Assignments=NMFClustering(DGE,Beads,GeneNames,'k',20,'VariableGenes','YOURPATH\BeadSeq Code\VariableGeneLists\F_GRCm38.81.P60Hippocampus.selected_genes.txt');
>>PlotGroupsInSpace(Beads,Assignments)

%Discover variable genes de novo and do NMF:
>>Assignments=NMFClustering(DGE,Beads,GeneNames,'k',20);
>>PlotGroupsInSpace(Beads,Assignments)

%Find spatially significant genes
>>output=PermutationTestByCluster(GetPuckDirectory('Puck_180620_4'),5,'FilterByCutoff',0.5)

%Get the subset of genes assigned to the CA2/3 cluster:
>>[ClusterDGE,ClusterBeads,ClusterBarcodes,ClusterGeneNames]=DGEByCluster(GetPuckDirectory('Puck_180618_4'),6);
%Plot the expression of Pcp4 using the subsetted DGE
>>PlotGeneFromName('Pcp4',ClusterGeneNames,ClusterDGE,ClusterBeads);
%Plot as a heatmap:
>>PlotGeneFromName('Pcp4',ClusterGeneNames,ClusterDGE,ClusterBeads,'PlotStyle','Imagesc');
%Plot the expression of Pcp4 from the full DGE:
>>PlotGeneFromName('Pcp4',GeneNames,DGE,Beads,'Cluster',6,'PuckName','Puck_180618_4');

%Find all genes that overlap significantly with Kctd12, Aldoc, and Prkcd, separately, on each of four hippocampal pucks. This will take some time to run:
>>[Cooccurrences,CooccuringGenes,Anticooccurrences, AnticooccurringGenes,GoodGeneList]=SignificanceOverlap(["Kctd12","Aldoc","Prkcd"],{"Puck_180819_10","Puck_180819_11","Puck_180819_12"});
%Find the genes that correlate with at least two of the three genes on at least two of the three pucks:
>>CorrelatingGenes(find(sum(sum(Cooccurrences,2)>=2,3)>=2))