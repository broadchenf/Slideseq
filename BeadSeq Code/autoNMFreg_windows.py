
# coding: utf-8

# In[1]:


import argparse 
parser = argparse.ArgumentParser(description = "handle inputs from sam matlab script to run NMFreg on puck data")
parser.add_argument("-da", type=str,
                   help = "Pass the data path for the atlas here")
parser.add_argument("-dp", type=str,
                   help = "Pass the datapath for the puck here")
parser.add_argument("-t", type=str,
                   help = "string of tissue type as is in data reference directory eg hippocampus, cerebellum, etc")
parser.add_argument("-c", type=int,
                   help = "cutoff for the UMI filtering of the DGE for the puck")
parser.add_argument("-dge", type=str,
                   help = "name of the dge file. No extension assums csv")
parser.add_argument("-bl", type=str,
                   help = "name of the bead location file passed here. No extension assumes csv")
print(parser)



# In[2]:


args = parser.parse_args()

#args = parser.parse_args('-da /Volumes/broad_macosko/bstickels/data/slideseq/NMFreg -dp /Volumes/broad_macosko/bstickels/data/slideseq/NMFreg/autoseq_trial -t cerebellum -c 100'.split())
#vars(args)


# In[3]:


import pandas as pd
import numpy as np
from IPython.display import display
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
from sklearn.preprocessing import StandardScaler
from IPython.display import display
import datetime

now = datetime.datetime.now()
import scipy.optimize
import scipy.stats



####Enter names of path###
#atlas data_path
data_path = args.da

#puck data path
data_path_puck = args.dp
#data_path = /Users/tkamath/Macosko_Lab/NMFreg/NMFreg/data

tissue_name = args.t
#datapath for the puck

#puck_number = '180430_1'
#puck_data_path = '/Volumes/broad_macosko/aleks/NMFreg/data/seqpuck{}_data'.format(puck_number)

#UMI threshold for the cutoff
UMI_threshold = args.c

puck_dge_name = args.dge

bead_locations = args.bl

tissue_data_path = "{}\\{}".format(data_path,tissue_name)
#if not os.path.exists(tissue_data_path):
#	os.makedirs(tissue_data_path)
print(data_path)
print(data_path_puck)
print(tissue_name)
print(UMI_threshold)
print(tissue_data_path)
print(sys.version)


# In[4]:


#### Define saving function #### 
def save_result(name):
    plt.savefig("{}{}.png".format(NMFreg_output, name),
                bbox_inches='tight', transparent=True, dpi=None)


# In[5]:

NMFreg_output = '{}\\Analogizer_NMFreg_output\\'.format(data_path_puck)
if not os.path.exists(NMFreg_output):
    os.makedirs(NMFreg_output)
print(NMFreg_output)
f=open("{}Log.txt".format(NMFreg_output),'w')
f.write("Beginning to read in Slideseq DGE\n")
f.close()
    
#### Read in count and coordinate data ####
dge_path = "{}\\{}.csv".format(data_path_puck,puck_dge_name)
dge = pd.read_csv(dge_path, header = 0, index_col = 0)
dge = dge.T
dge = dge.reset_index()
dge = dge.rename(columns={'index':'barcode'})

f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Slideseq DGE Read-in Complete\n")
f.close()
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Beginning to read in Slideseq Locations\n")
f.close()
coords = pd.read_csv("{}\\{}.csv".format(data_path_puck,bead_locations), header = 0)
coords = coords.rename(columns={'Barcodes':'barcode'})
coords = coords.rename(columns={'barcodes':'barcode'})
df_merged = dge.merge(coords, right_on='barcode', left_on='barcode')
counts = df_merged.drop(['xcoord', 'ycoord'], axis=1)
#coords.shape
#dge = []
#df_merged = []
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Slideseq Locations read-in complete\n")
f.close()


### Read in atlas data ###
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Beginning to read in Atlas data\n")
f.close()

atlas_dge = pd.read_csv("{}\\dge_hvgs.csv".format(tissue_data_path))
print(atlas_dge.shape)
#atlas_dge.head(5)

cell_clusters = pd.read_csv("{}\\cell_cluster_outcome.csv".format(tissue_data_path))
print(cell_clusters.shape)
#cell_clusters.head(5)
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Atlas data read-in complete\n")
f.close()

##Take only the genes in both the atlas and the slideseq
atlas_genes = atlas_dge.columns.tolist()
puck_genes = counts.columns.tolist()[1:]

gene_intersection = list(set(atlas_genes) & set(puck_genes))

atlasdge = atlas_dge[gene_intersection]

counts = counts[['barcode'] + gene_intersection]

##Now that we have taken the intersection of atlas and slideseq genes, take only those beads with more than the specified number of transcripts in the intersection.

sample_info = df_merged[['barcode','xcoord', 'ycoord']]

counts_barcodestotals = np.sum(counts.drop('barcode', axis=1), axis=1)
print(counts_barcodestotals.shape)
sample_info['total_counts'] = counts_barcodestotals

sample_info_grthreshold = sample_info.loc[sample_info['total_counts'] > UMI_threshold]
coords = sample_info_grthreshold
df_merged_grthreshold = counts.merge(sample_info_grthreshold, right_on='barcode', left_on='barcode')
counts_grthreshold = df_merged_grthreshold.drop(['xcoord', 'ycoord', 'total_counts'], axis=1)

counts = counts_grthreshold
df_merged = []
df_merged_grthreshold = []
counts_grthreshold = []
print(counts.shape)
print(counts[0:4].sum(axis=1))

##Because we dropped beads, we now have some genes with no reads. We have to drop them as well

counts_nob = counts.drop("barcode", axis=1)
print(counts_nob.shape)
#compute the number of UMIs per gene
UMIspergene = np.sum(counts_nob, axis=0)
#genes with more than 0 UMIs
print(UMIspergene.shape)
#get the genes with more than 0 UMI
gr0_UMIspergene = UMIspergene[1:][UMIspergene[1:] > 0]
print(gr0_UMIspergene.shape)
print(gr0_UMIspergene.head(5))
counts = counts[['barcode'] + gr0_UMIspergene.index.tolist()]
print(counts.shape)

#### Plot count/coordinate data ####
plt.figure(figsize = (10, 10))
plt.scatter(coords['xcoord'], coords['ycoord'], c='k', s=4, alpha=0.6);
plt.axis('equal');
save_result("filtered_tissue_coverage")
plt.close()

#plt.figure(1)
plt.figure(figsize = (12, 12))
plt.set_cmap('viridis_r')
plt.scatter(coords['xcoord'], coords['ycoord'], c=coords['total_counts'], s=4, alpha=0.6);
plt.axis('equal');
plt.colorbar();
save_result("bead_counts")
plt.close()


#Because we dropped some puck genes, we have to subset the atlas DGE again.
atlas_genes = atlas_dge.columns.tolist()
puck_genes = counts.columns.tolist()[1:]

gene_intersection = list(set(atlas_genes) & set(puck_genes))

atlasdge = atlas_dge[gene_intersection]



puckcounts = counts.set_index(counts['barcode'])
puckcounts = puckcounts.drop('barcode', axis=1)

cell_totalUMI = np.sum(puckcounts, axis = 1)
puckcounts_cellnorm = np.true_divide(puckcounts, cell_totalUMI[:,None])

puckcounts_scaled = StandardScaler(with_mean=False).fit_transform(puckcounts_cellnorm)
print(puckcounts_scaled.shape)

# In[19]:



cell_totalUMIa = np.sum(atlasdge, axis = 1)
atlasdge_cellnorm = np.true_divide(atlasdge, cell_totalUMIa[:,None])
atlasdge_scaled = StandardScaler(with_mean=False).fit_transform(atlasdge_cellnorm)


# In[20]:


### Perform NMF on the atlas data to find basis for projection #####
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Performing Atlas NMF\n")
f.close()

from sklearn.decomposition import NMF
K = 30
alpha = 0
l1_ratio = 0
random_state = 17
model = NMF(n_components=K, init='random', random_state = random_state, alpha = alpha, l1_ratio = l1_ratio)
Ha = model.fit_transform(atlasdge_scaled)
Wa = model.components_

f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Atlas NMF Complete\n")
f.close()

## Save output

np.savetxt("{}Ha{}_{}_{}_{}.csv".format(NMFreg_output, K, alpha, l1_ratio, random_state), Ha, delimiter=",")
np.savetxt("{}Wa{}_{}_{}_{}.csv".format(NMFreg_output, K, alpha, l1_ratio, random_state), Wa, delimiter=",")

Ha_norm = StandardScaler(with_mean=False).fit_transform(Ha)
Ha_norm = pd.DataFrame(Ha_norm)
Ha_norm['barcode'] = atlasdge.index.tolist()
maxloc = Ha_norm.drop('barcode', axis=1).values.argmax(axis=1)
cell_clusters['maxloc'] = maxloc
cell_clusters.head()


# In[21]:


#num_atlas_clusters = np.unique(cell_clusters['cluster']).size
#### Interpret these factors to cell type assignments #####
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Mapping Atlas Cells to Clusters\n")
f.close()

num_atlas_clusters = max(cell_clusters['cluster'])
bins = num_atlas_clusters
factor_to_celltype_df = pd.DataFrame(0, index=range(1, num_atlas_clusters+1), columns=range(K))
print(factor_to_celltype_df)
for k in range(K):
    n, bins, patches = plt.hist(cell_clusters['cluster'][cell_clusters['maxloc'] == k],
            bins, range = (0.5,(bins+0.5)), facecolor='green', alpha=0.75)
    print(n)
    factor_to_celltype_df[k] = n.astype(int)

factor_to_celltype_df = factor_to_celltype_df.T

factor_total = np.sum(factor_to_celltype_df, axis = 1)
factor_to_celltype_df_norm = np.true_divide(factor_to_celltype_df, factor_total[:,None])

f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Atlas Cluster Mapping Complete\n")
f.close()


plt.figure()

cx = sns.clustermap(factor_to_celltype_df_norm, fmt = 'd',
                cmap="magma_r", linewidth=0.5, col_cluster = False,
                   figsize=(7, 9))
ax = sns.clustermap(factor_to_celltype_df_norm, fmt = 'd',
                cmap="magma_r", linewidth=0.5, col_cluster = False,
                   annot = factor_to_celltype_df.loc[cx.dendrogram_row.reordered_ind],
                   figsize=(7, 9))
save_result("AtlasClusterMapping")
plt.close()


# In[23]:


from tqdm import tqdm

###### Get a factor to cell type and cell type to factor dictionaries #######
maxloc_fc = factor_to_celltype_df.values.argmax(axis=1)
factor_to_celltype_dict = {factor : ctype + 1 for factor, ctype in enumerate(maxloc_fc)}

celltype_to_factor_dict = {}
for c in range(1, num_atlas_clusters + 1):
    celltype_to_factor_dict[c] = [k for k, v in factor_to_celltype_dict.items() if v == c]



##### Perform NNLS with the atlas basis #######
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Performing NNLS\n")
f.close()

WaT = Wa.T
XsT = puckcounts_scaled.T

Hs_hat = []
for b in tqdm(range(XsT.shape[1])):
    h_hat = scipy.optimize.nnls(WaT, XsT[:, b])[0]
    if b == 0:
        Hs_hat = h_hat
    else:
        Hs_hat = np.vstack((Hs_hat, h_hat))

Ha = pd.DataFrame(Ha)
Ha['cellname'] = atlasdge.index.tolist()
Ha_indexed = Ha.set_index('cellname')

Hs = pd.DataFrame(Hs_hat)
Hs['barcode'] = puckcounts.index.tolist()
Hs_indexed = Hs.set_index('barcode')

Hs_indexed.to_csv("{}Hs{}_{}_{}_{}.csv".format(NMFreg_output, K, alpha, l1_ratio, random_state), index=True)

f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("NNLS complete\n")
f.close()

## Scale the Ha and Hs matrices to unit variance

Ha_norm = pd.DataFrame(StandardScaler(with_mean=False).fit_transform(Ha_indexed))
Hs_norm = pd.DataFrame(StandardScaler(with_mean=False).fit_transform(Hs_indexed))

## Re-add indices after scaling

Ha_norm['cellname'] = atlasdge.index.tolist()
Ha_norm_indexed = Ha_norm.set_index('cellname')

Hs_norm['barcode'] = puckcounts.index.tolist()
Hs_norm_indexed = Hs_norm.set_index('barcode')



##### Assign barcodes a cluster assignment from the atlas data set based on max factors
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Assigning slideseq beads to clusters\n")
f.close()

maxloc_s = Hs_norm_indexed.values.argmax(axis=1)
barcode_clusters = pd.DataFrame()
barcode_clusters['barcode'] = Hs_norm_indexed.index.tolist()
barcode_clusters['max_factor'] = maxloc_s

barcode_clusters['atlas_cluster'] = barcode_clusters['barcode']

for c in range(1, num_atlas_clusters + 1):
    condition = np.isin(barcode_clusters['max_factor'], celltype_to_factor_dict[c])
    barcode_clusters['atlas_cluster'][condition] = c
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Slideseq bead cluster assignment complete\n")
f.close()

#########################################################

barcode_clusters.to_csv("{}\\AnalogizerClusterAssignments.csv".format(data_path_puck), index=True)

####### Plot the assignments onto the puck ########
f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Plotting assignments onto the puck\n")
f.close()

from matplotlib import colors

for i in range(1,1+max(barcode_clusters['atlas_cluster'])):
    boolcol = (barcode_clusters['atlas_cluster']==i)
    sub_df = barcode_clusters.copy()
    sub_df['bool'] = boolcol

    plt.figure(figsize =(10, 10))
    plt.set_cmap('copper_r')
    plt.scatter(coords['xcoord'], coords['ycoord'], c=sub_df['bool'], s=4, alpha=0.6)
    plt.title(tissue_name + 'atlas cluster' + str(i))
    plt.axis('equal')
    save_result(tissue_name + str(i) + '_puckplots')
    plt.close()

f=open("{}Log.txt".format(NMFreg_output),'a')
f.write("Script complete\n")
f.close()

exit()
