'''
This script performs standard quality control on the Cell Ranger output data. 
'''

#import modules that will be needed
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import pickle
from sklearn.cluster import KMeans


#########################################################Main Thread####################################

#Set file locations
qc_dir = '/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/quality_control/'
cellranger_loc = '/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs/'
cellbender_loc = '/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellbender_outputs/learning_rate_0.00005/'
grna_ref_loc = cellranger_loc + 'aggr/crispr_analysis/protospacer_calls_per_cell.csv'
#Set the pool list
pools = ['Pool_' + str(x) for x in range(1,9,1)]

#Set scanpy parameters
sc.settings.set_figure_params(dpi=80, facecolor="white")

#Define a function to do basic qc plotting
def plot_basic_qc(cellranger_raw, pool, file_prefix):
    
    #Set output directory for figures
    sc.settings.figdir = qc_dir + pool + '/'
    
    #Make names of gene features unique 
    cellranger_raw.var_names_make_unique()
    #Basic filtering
    hfob_filt = cellranger_raw
    sc.pp.filter_cells(hfob_filt, min_genes=500) #Every cell has to have at least 500 genes (I can bump this up later)
    sc.pp.filter_genes(hfob_filt, min_cells=3) #Every gene has to be in at least 3 cells
    
    #Make plot of most highly expressed genes
    sc.pl.highest_expr_genes(hfob_filt, n_top=20, show=False, save='.' + file_prefix)
    
    #Calculate data on mitochondrial reads
    hfob_filt.var['mt'] = hfob_filt.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(hfob_filt, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    #Make qc plots
    sc.pl.violin(hfob_filt, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                 stripplot=False, multi_panel=True, show=False, save='.' + file_prefix)
    sc.pl.scatter(hfob_filt, x='total_counts', y='pct_counts_mt', show=False, save=".total_UMIs.mito_pct." + file_prefix)
    sc.pl.scatter(hfob_filt, x='total_counts', y='n_genes_by_counts', show=False, save=".total_UMIs.genes_per_cell." + file_prefix)
    
    #Plot cells without removing any additional cells
    sc.pp.normalize_total(hfob_filt, target_sum=1e4)
    sc.pp.log1p(hfob_filt)
    sc.pp.highly_variable_genes(hfob_filt, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(hfob_filt, show=False, save='.' + file_prefix)
    hfob_filt.raw = hfob_filt
    sc.tl.pca(hfob_filt, svd_solver='arpack')
    sc.pl.pca_variance_ratio(hfob_filt, log=True,show=False, save=".no_filt." + file_prefix)
    sc.pp.neighbors(hfob_filt, n_neighbors=10, n_pcs=40)
    sc.tl.umap(hfob_filt)
    sc.pl.umap(hfob_filt, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], show=False, save=".no_filt." + file_prefix)
    
    #Binarize qc values and store in observations, then replot
    hfob_filt.obs['genes_by_counts_lt_2k'] = (hfob_filt.obs['n_genes_by_counts'] < 2000) * 1
    hfob_filt.obs['genes_by_counts_lt_2.5k'] = (hfob_filt.obs['n_genes_by_counts'] < 2500) * 1
    hfob_filt.obs['genes_by_counts_lt_3k'] = (hfob_filt.obs['n_genes_by_counts'] < 3000) * 1
    hfob_filt.obs['genes_by_counts_lt_3.5k'] = (hfob_filt.obs['n_genes_by_counts'] < 3500) * 1
    hfob_filt.obs['genes_by_counts_lt_4k'] = (hfob_filt.obs['n_genes_by_counts'] < 4000) * 1
    hfob_filt.obs['total_counts_lt_8k'] = (hfob_filt.obs['total_counts'] < 8000) * 1
    hfob_filt.obs['total_counts_lt_7k'] = (hfob_filt.obs['total_counts'] < 7000) * 1
    hfob_filt.obs['total_counts_lt_6k'] = (hfob_filt.obs['total_counts'] < 6000) * 1
    hfob_filt.obs['total_counts_lt_5k'] = (hfob_filt.obs['total_counts'] < 5000) * 1
    hfob_filt.obs['total_counts_lt_4k'] = (hfob_filt.obs['total_counts'] < 4000) * 1
    hfob_filt.obs['pct_counts_mt_gt_5'] = (hfob_filt.obs['pct_counts_mt'] > 5) * 1
    hfob_filt.obs['pct_counts_mt_gt_8'] = (hfob_filt.obs['pct_counts_mt'] > 8) * 1
    hfob_filt.obs['pct_counts_mt_gt_10'] = (hfob_filt.obs['pct_counts_mt'] > 10) * 1
    hfob_filt.obs['pct_counts_mt_gt_12'] = (hfob_filt.obs['pct_counts_mt'] > 12) * 1
    hfob_filt.obs['pct_counts_mt_gt_15'] = (hfob_filt.obs['pct_counts_mt'] > 15) * 1
    sc.pl.umap(hfob_filt, color=['total_counts_lt_8k', 'total_counts_lt_7k', 'total_counts_lt_6k', 'total_counts_lt_5k', 'total_counts_lt_4k'], show=False, save=".no_filt.binarized.total_counts_descending." + file_prefix)
    sc.pl.umap(hfob_filt, color=['genes_by_counts_lt_4k', 'genes_by_counts_lt_3.5k', 'genes_by_counts_lt_3k', 'genes_by_counts_lt_2.5k', 'genes_by_counts_lt_2k'], show=False, save=".no_filt.binarized.gene_counts_descending." + file_prefix)
    sc.pl.umap(hfob_filt, color=['pct_counts_mt_gt_5', 'pct_counts_mt_gt_8', 'pct_counts_mt_gt_10', 'pct_counts_mt_gt_12', 'pct_counts_mt_gt_15'], show=False, save=".no_filt.binarized.mito_pcts_ascending." + file_prefix)

#Read in files and run qc tests
cellranger_results = {}
cellbender_results = {}
for pool in pools:
    cellranger_results[pool] = sc.read_10x_h5(cellranger_loc + pool + '/filtered_feature_bc_matrix.h5')
    cellbender_results[pool] = sc.read_10x_h5(cellbender_loc + pool + '/cell_bender.output_filtered.h5')
    plot_basic_qc(cellranger_results[pool], pool, pool + '_pre_cellbender')
    plot_basic_qc(cellbender_results[pool], pool, pool + '_post_cellbender')

#Get IDs of the good pools (This is being done manually based on the QC data)
good_pools = pools[:]
good_pools.remove('Pool_2')

#Make K-Means Clustering Plots for a K of 2 for each Pool with the CellBender Results
for pool in pools:
    temp = cellbender_results[pool][:]
    #Set output directory for figures
    sc.settings.figdir = qc_dir + pool + '/'
    kmeans = KMeans(n_clusters=2, random_state=0).fit(temp.obsm['X_pca']) 
    temp.obs['kmeans2'] = kmeans.labels_.astype(str)
    sc.pl.umap(temp, color=['kmeans2'], show=False, save=".no_filt.2-means_cluster")

#Make Leiden Clustering Plot with the CellBender Results
for pool in pools:
    temp = cellbender_results[pool]
    #Set output directory for figures
    sc.settings.figdir = qc_dir + pool + '/'
    sc.tl.leiden(temp, resolution=0.05, random_state=0)
    sc.pl.umap(temp, color=['leiden'], show=False, save=".no_filt.leiden_cluster")
    
#Create a dictionary of filtered cellbender results filtered on the Leiden results
filtered_cellbender_results = cellbender_results.copy()
del filtered_cellbender_results['Pool_2']
for pool in good_pools:
    #Do the filtering
    filtered_cellbender_results[pool] = filtered_cellbender_results[pool][filtered_cellbender_results[pool].obs['leiden'] == '0']
    #Rename the cellular barcodes in each pool to match the pool name
    filtered_cellbender_results[pool].obs_names = [x.replace('-1','-' + pool.replace('Pool_', '')) for x in filtered_cellbender_results[pool].obs_names]
    #Plot qc metrics
    sc.settings.figdir = qc_dir + pool + '/'
    sc.pl.violin(filtered_cellbender_results[pool], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                 stripplot=False, multi_panel=True, show=False, save='.post_leiden')

####################################### Development Checkpoint ########################################
#outfile = open(qc_dir + 'aggr/leiden_filtered_results.pkl','wb')
#pickle.dump(filtered_cellbender_results,outfile)
#outfile.close()
#Read in pickle file
with open(qc_dir + 'aggr/leiden_filtered_results.pkl', 'rb') as f:
    filtered_cellbender_results = pickle.load(f)
######################################################################################################

#Read in cellbender output files again and filter them for the good cells
cellbender_good_results = {}
for pool in good_pools:
    #Read in file again
    cellbender_good_results[pool] = sc.read_10x_h5(cellbender_loc + pool + '/cell_bender.output_filtered.h5')
    #Rename droplet ids and store pool info
    cellbender_good_results[pool].obs_names = [x.replace('-1','-' + pool.replace('Pool_', '')) for x in cellbender_good_results[pool].obs_names]
    cellbender_good_results[pool].obs['Pool'] = pool.replace('Pool_', '')
    #Use ENSG codes for gene names
    cellbender_good_results[pool].var_names=cellbender_good_results[pool].var['gene_ids']
    #Filter for the likely cell-containing droplets
    cellbender_good_results[pool] = cellbender_good_results[pool][filtered_cellbender_results[pool].obs_names]
    #Make variable names unique
    cellbender_good_results[pool].var_names_make_unique()

#Merge the good cells into a single scanpy object
cellbender_merged_good = ad.concat(cellbender_good_results.values())

#Run basic QC metrics for the merged good cells
sc.settings.figdir = qc_dir + 'aggr/'
#Calculate data on mitochondrial reads
cellbender_merged_good.var['mt'] = cellbender_merged_good.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(cellbender_merged_good, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#Make qc plots
sc.pl.violin(cellbender_merged_good, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
             stripplot=False, multi_panel=True, show=False, save=True)

#Read in the csv reference file of gRNAs as a pandas dataframe
grna_ref_raw = pd.read_table(grna_ref_loc, sep=",")
#Append on the information about whether the cell has been retained thus far
grna_ref_raw['passed_init_qc'] = [True if x in cellbender_merged_good.obs_names.to_list() else False for x in grna_ref_raw['cell_barcode']]
#Create lists of guide RNAs present in cells with 1 sgRNA and those with multiples
single_sgRNAs = grna_ref_raw[(grna_ref_raw['num_features'] == 1) & (grna_ref_raw['passed_init_qc'])]['feature_call']
single_sgRNAs_unique, single_sgRNAs_counts = np.unique(single_sgRNAs, return_counts=True)
multi_sgRNAs = grna_ref_raw[(grna_ref_raw['passed_init_qc'])]['feature_call']
multi_sgRNAs = [x.split('|') for x in multi_sgRNAs]
multi_sgRNAs = [item for sublist in multi_sgRNAs for item in sublist]
multi_sgRNAs_unique, multi_sgRNAs_counts = np.unique(multi_sgRNAs, return_counts=True)
#Make histograms of the sgRNA counts
fig, ax = plt.subplots()
ax.hist(single_sgRNAs_counts, bins=20)
ax.axvline(np.median(single_sgRNAs_counts), color='k', linestyle='dashed', linewidth=1)
fig.tight_layout()
plt.savefig(qc_dir + 'aggr/single_sgRNA_counts.hist.png')
plt.close()
fig, ax = plt.subplots()
ax.hist(multi_sgRNAs_counts, bins=20)
ax.axvline(np.median(multi_sgRNAs_counts), color='k', linestyle='dashed', linewidth=1)
fig.tight_layout()
plt.savefig(qc_dir + 'aggr/multi_sgRNA_counts.hist.png')
plt.close()
#Calculate basic statistics
np.median(single_sgRNAs_counts)
np.min(single_sgRNAs_counts)
np.median(multi_sgRNAs_counts)
np.min(multi_sgRNAs_counts)

#Run qc process on aggregated data
cellbender_merged_good_copy = cellbender_merged_good[:]
sc.pp.filter_cells(cellbender_merged_good_copy, min_genes=500)
sc.pp.normalize_total(cellbender_merged_good_copy, target_sum=1e4)
sc.pp.log1p(cellbender_merged_good_copy)
sc.pp.highly_variable_genes(cellbender_merged_good_copy, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(cellbender_merged_good_copy, svd_solver='arpack')
sc.pl.pca_variance_ratio(cellbender_merged_good_copy, log=True,show=False, save=".no_filt.aggr.filtered")
sc.pp.neighbors(cellbender_merged_good_copy, n_neighbors=10, n_pcs=40)
sc.tl.umap(cellbender_merged_good_copy)
sc.pl.umap(cellbender_merged_good_copy, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'Pool'], show=False, save=".no_filt.aggr.filtered")

#Filter the data to remove likely doublets and dying cells
cellbender_merged_good_filtered = cellbender_merged_good[cellbender_merged_good.obs.total_counts < 90000, :]
cellbender_merged_good_filtered = cellbender_merged_good_filtered[cellbender_merged_good_filtered.obs.pct_counts_mt < 10, :]
#Create a version to save for export
cellbender_export = cellbender_merged_good_filtered[:]

#Append on the information about whether the cell has been retained thus far
grna_ref_raw['passed_secondary_qc'] = [True if x in cellbender_merged_good_filtered.obs_names.to_list() else False for x in grna_ref_raw['cell_barcode']]
#Create lists of guide RNAs present in cells with 1 sgRNA and those with multiples
secondary_single_sgRNAs = grna_ref_raw[(grna_ref_raw['num_features'] == 1) & (grna_ref_raw['passed_secondary_qc'])]['feature_call']
secondary_single_sgRNAs_unique, secondary_single_sgRNAs_counts = np.unique(secondary_single_sgRNAs, return_counts=True)
secondary_multi_sgRNAs = grna_ref_raw[(grna_ref_raw['passed_secondary_qc'])]['feature_call']
secondary_multi_sgRNAs = [x.split('|') for x in secondary_multi_sgRNAs]
secondary_multi_sgRNAs = [item for sublist in secondary_multi_sgRNAs for item in sublist]
secondary_multi_sgRNAs_unique, secondary_multi_sgRNAs_counts = np.unique(secondary_multi_sgRNAs, return_counts=True)
#Make histograms of the sgRNA counts
fig, ax = plt.subplots()
ax.hist(secondary_single_sgRNAs_counts, bins=20)
ax.axvline(np.median(secondary_single_sgRNAs_counts), color='k', linestyle='dashed', linewidth=1)
fig.tight_layout()
plt.savefig(qc_dir + 'aggr/single_sgRNA_counts.post_secondary_qc.hist.png')
plt.close()
fig, ax = plt.subplots()
ax.hist(secondary_multi_sgRNAs_counts, bins=20)
ax.axvline(np.median(secondary_multi_sgRNAs_counts), color='k', linestyle='dashed', linewidth=1)
fig.tight_layout()
plt.savefig(qc_dir + 'aggr/multi_sgRNA_counts.post_secondary_qc.hist.png')
plt.close()
#Calculate basic statistics
np.median(secondary_single_sgRNAs_counts)
np.min(secondary_single_sgRNAs_counts)
np.median(secondary_multi_sgRNAs_counts)
np.min(secondary_multi_sgRNAs_counts)

#Run qc process on post-secondary qc data
sc.pp.filter_cells(cellbender_merged_good_filtered, min_genes=500)
sc.pp.normalize_total(cellbender_merged_good_filtered, target_sum=1e4)
sc.pp.log1p(cellbender_merged_good_filtered)
sc.pp.highly_variable_genes(cellbender_merged_good_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(cellbender_merged_good_filtered, svd_solver='arpack')
sc.pl.pca_variance_ratio(cellbender_merged_good_filtered, log=True,show=False, save=".no_filt.aggr.filtered")
sc.pp.neighbors(cellbender_merged_good_filtered, n_neighbors=10, n_pcs=40)
sc.tl.umap(cellbender_merged_good_filtered)
sc.pl.umap(cellbender_merged_good_filtered, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'Pool'], show=False, save=".no_filt.aggr.post_secondary_qc")

#Write exportable annData object to file for use in Seurat/R
cellbender_export.write_h5ad(qc_dir + 'aggr/aggr.post_qc.h5ad')