'''
This script performs standard quality control on the Cell Ranger output data. 
'''

#import modules that will be needed
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pickle
import pybiomart

#########################################################Main Thread####################################

#Set crispri screen file locations
qc_dir = '/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/quality_control/'
filtered_cell_loc = os.path.join(qc_dir, "aggr", "all_cells.fully_qced.pkl")
grna_ref_loc = '/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/protospacer_calls_per_cell.csv'

#Set bulk hFOB tpm location
bulk_loc = "/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/DE/hFOB/hFOB_tpm.txt"
hmsc_loc = "/mnt/isilon/sfgi/chesia/analyses/grant/rnaSeq/BMP2/BMP2_vs_Contr.txt"
#Set gene length file
gene_len_loc = "/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.length.txt"
#Set file of osteoblast marker genes
hfob_loc = qc_dir + 'hFOB_marker_genes.csv'

#Set scanpy parameters
sc.settings.set_figure_params(dpi=80, facecolor="white")

#Read in bulk hFOB RNA-seq results
bulk_raw = pd.read_csv(bulk_loc, sep="\t")
#Filter for ENSG genes
bulk_filt = bulk_raw.loc[[True if "ENSG" in x else False for x in bulk_raw.gene_id.to_list()],:]
bulk_filt['cut_id'] = [x.split('.')[0] for x in bulk_filt.gene_id.to_list()]

#Read in bullk hMSC results
hmsc_raw = pd.read_csv(hmsc_loc, sep="\t")
hmsc_filt = hmsc_raw.loc[[True if "ENSG" in x else False for x in hmsc_raw.GENE.to_list()],:]
hmsc_filt['cut_id'] = [x.split('.')[0] for x in hmsc_filt.GENE.to_list()]
hmsc_filt.index = hmsc_filt['cut_id']
#Convert FPKM to TPM
for fpkm_col in [x for x in hmsc_filt.columns if "FPKM" in x]:
    scale_factor = hmsc_filt[fpkm_col].sum()
    hmsc_filt[fpkm_col.replace("FPKM", "TPM")] = (hmsc_filt[fpkm_col] / scale_factor) * 1e6

#Read in hfob gene reference
hfob_ref_raw = pd.read_table(hfob_loc, sep=",")
#Make dictionary mapping and get maturation and mineralization genes
osteo_mapping = {hfob_ref_raw['ENSG'][x] : hfob_ref_raw['Symbol'][x] for x in range(len(hfob_ref_raw))}
proliferation_genes = [hfob_ref_raw['Symbol'][x] for x in range(len(hfob_ref_raw)) if hfob_ref_raw['Category'][x] == "Proliferation"]
maturation_genes = [hfob_ref_raw['Symbol'][x] for x in range(len(hfob_ref_raw)) if hfob_ref_raw['Category'][x] == "Maturation"]
mineralization_genes = [hfob_ref_raw['Symbol'][x] for x in range(len(hfob_ref_raw)) if hfob_ref_raw['Category'][x] == "Mineralization"]

#Read in pickle file of screen cells
with open(filtered_cell_loc, 'rb') as f:
    screen_adata = pickle.load(f)

#Read in file of grna information
grna_ref_raw = pd.read_table(grna_ref_loc, sep=",")
grna_ref_raw.index = grna_ref_raw.cell_barcode
#Append on grna info o the anndata object
temp = pd.concat([pd.Series([grna_ref_raw.loc[x,"feature_call"] if x in grna_ref_raw.index else "" for x in screen_adata.obs.index]),pd.Series([grna_ref_raw.loc[x,"num_features"] if x in grna_ref_raw.index else 0 for x in screen_adata.obs.index])], axis=1)
temp.columns = ["guide_name", "guide_count"]
temp.index = screen_adata.obs.index
screen_adata.obs = pd.concat([screen_adata.obs, temp], axis=1)

#Filter for the no guide and non-targeting guide cells
no_guide_cells = screen_adata[screen_adata.obs['guide_count'] == 0]
nt_cells = screen_adata[[True if screen_adata.obs.guide_count[x] != 0 and 'rs' not in screen_adata.obs.guide_name[x] else False for x in range(len(screen_adata.obs))]]
#Create pseudo bulk vectors of expression for each set and combine into dataframe
no_guide_pseudo = no_guide_cells.X.sum(axis=0).A1 if isinstance(no_guide_cells.X, np.matrix) else no_guide_cells.X.sum(axis=0)
no_guide_pseudo = np.array(no_guide_pseudo).flatten() 
nt_pseudo = nt_cells.X.sum(axis=0).A1 if isinstance(nt_cells.X, np.matrix) else nt_cells.X.sum(axis=0)
nt_pseudo = np.array(nt_pseudo).flatten() 
pseudo_bulk_df = pd.DataFrame({
    'No_Guide_Expr': no_guide_pseudo,
    'NT_Expr': nt_pseudo
})
pseudo_bulk_df.index = screen_adata.var_names

#Read in file of gene lengths to convert pseudobulk counts to TPMs
gene_len_raw = pd.read_csv(gene_len_loc, sep="\t", header=None)
gene_len_raw.columns = ['raw_id', 'symbol', 'bp']
gene_len_raw['cut_id'] = [x.split('.')[0] for x in gene_len_raw.raw_id.to_list()]
gene_len_raw.index = gene_len_raw['cut_id']
#Extract gene lengths for pseudobulk genes
gene_lens = [gene_len_raw.loc[x, 'bp'] if x in gene_len_raw.index else None for x in pseudo_bulk_df.index]
#Filter the pseudobulk_df for genes with non-missing expression
pseudo_bulk_df_filt = pseudo_bulk_df.loc[[True if x != None else False for x in gene_lens],:]
gene_lens_filt = pd.Series([x for x in gene_lens if x != None])
gene_lens_filt.index = pseudo_bulk_df_filt.index

##### Calculate TPM for pseudobulk expression #####
#Convert gene lengths to kb
gene_lengths_kb = gene_lens_filt / 1000
# Normalize the counts by gene length for both conditions
pseudo_bulk_df_filt['No_Guide_Expr_norm'] = pseudo_bulk_df_filt['No_Guide_Expr'] / gene_lengths_kb
pseudo_bulk_df_filt['NT_Expr_norm'] = pseudo_bulk_df_filt['NT_Expr'] / gene_lengths_kb
# Calculate scaling factors (sum of all length-normalized counts)
no_guide_scaling_factor = pseudo_bulk_df_filt['No_Guide_Expr_norm'].sum()
nt_scaling_factor = pseudo_bulk_df_filt['NT_Expr_norm'].sum()
# Calculate TPM
pseudo_bulk_df_filt['No_Guide_TPM'] = (pseudo_bulk_df_filt['No_Guide_Expr_norm'] / no_guide_scaling_factor) * 1e6
pseudo_bulk_df_filt['NT_TPM'] = (pseudo_bulk_df_filt['NT_Expr_norm'] / nt_scaling_factor) * 1e6
# Drop the intermediate normalized columns
pseudo_bulk_df_filt = pseudo_bulk_df_filt.drop(columns=['No_Guide_Expr_norm', 'NT_Expr_norm'])
#Make a dataframe of tpms
tpm_df = pseudo_bulk_df_filt[['No_Guide_TPM', "NT_TPM"]]
temp = bulk_filt.cut_id.to_list()
tpm_df = tpm_df.loc[[True if x in temp else False for x in tpm_df.index],:]

#Add on the bulk hFOB data to the tpm_df
for repl in list(bulk_filt['sample'].unique()):
    #Get the data for the current bulk replicate
    temp = bulk_filt.query("sample == @repl")
    temp.index = temp.cut_id
    #Append on the new data
    tpm_df[repl] = temp.loc[list(tpm_df.index), 'tpm']
    
#Add on the bulk hMSC data to the tpm_df
tpm_df = tpm_df.loc[[x for x in tpm_df.index if x in hmsc_filt.index],:]
tpm_df = pd.concat([tpm_df, hmsc_filt.loc[tpm_df.index,[x for x in hmsc_filt.columns if "TPM" in x]]], axis = 1)
    
#### Make the PCA Plot #####
# Standardize the data
scaler = StandardScaler()
tpm_scaled = scaler.fit_transform(tpm_df.T)  # Transpose to have experiments as rows
# Apply PCA
pca = PCA(n_components=2)  # We are interested in the first two components
principal_components = pca.fit_transform(tpm_scaled)
# Create a DataFrame for the principal components
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=tpm_df.columns)
#Categorize the experiments
categories = []
for experiment in pca_df.index:
    if experiment.startswith("hFOBdiff_"):
        categories.append("Diff. hFOB")
    elif experiment.startswith("hFOBundiff_"):
        categories.append("Undiff. hFOB")
    elif "BMP2" in experiment:
        categories.append("hMSC-Osteo")
    elif "Contr" in experiment:
        categories.append("hMSC")
    else:
        categories.append("Screen")
        
pca_df['Category'] = categories 
# Define colors for each category
colors = {
    "Diff. hFOB": "#026602",
    "Undiff. hFOB": "#03fa03",
    "hMSC-Osteo": "#004778",
    "hMSC": "#0396fc",
    "Screen": "red"
}

# Plot the results
plt.figure(figsize=(8, 6))
for category in pca_df['Category'].unique():
    subset = pca_df[pca_df['Category'] == category]
    plt.scatter(subset['PC1'], subset['PC2'], label=category, color=colors[category], s=100)

plt.title('PCA of RNA-seq TPM Data')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title="Category")
plt.grid(True)

# Save the plot to the specified location
plt.savefig(os.path.join(qc_dir, 'aggr', 'bulk_pseudo-bulk-compare.pca.png'), dpi=300, bbox_inches='tight')

#Filter for marker genes
tpm_marker_df = tpm_df.loc[osteo_mapping.keys(),:]
tpm_marker_df.index = [osteo_mapping[x] for x in tpm_marker_df.index]
tpm_marker_df.columns = ["No Guide Screen", "NT Guide Screen", "Diff. Bulk 1", "Diff. Bulk 2", "Diff. Bulk 3", "Undiff. Bulk 1", "Undiff. Bulk 2", "Undiff. Bulk 3"]
#Melt the DataFrame to long format for plotting
tpm_marker_long = tpm_marker_df.reset_index().melt(id_vars='index', var_name='Experiment', value_name='TPM') 
# Rename 'index' to 'Gene' for clarity
tpm_marker_long = tpm_marker_long.rename(columns={'index': 'Gene'})
# Add a 'Category' column based on the gene lists
def assign_category(gene):
    if gene in proliferation_genes:
        return "Proliferation"
    elif gene in maturation_genes:
        return "Maturation"
    elif gene in mineralization_genes:
        return "Mineralization"
    else:
        return "Unknown"

tpm_marker_long['Category'] = tpm_marker_long['Gene'].apply(assign_category)
#Reset experiment names
tpm_marker_long['Experiment'] = tpm_marker_long.Experiment.str.replace("No_Guide_TPM", 'Screen No Guide')
tpm_marker_long['Experiment'] = tpm_marker_long.Experiment.str.replace("NT_TPM", 'Screen NT Guide')
tpm_marker_long['Experiment'] = tpm_marker_long.Experiment.str.replace("hFOBdiff_rep", 'Diff. hFOB ')
tpm_marker_long['Experiment'] = tpm_marker_long.Experiment.str.replace("hFOBundiff_rep", 'Undiff. hFOB ')
tpm_marker_long['Experiment'] = tpm_marker_long.Experiment.str.replace("TPM.BMP2", "hMSC-Osteo ")
tpm_marker_long['Experiment'] = tpm_marker_long.Experiment.str.replace("TPM.Contr", "hMSC ")

#Create a column for sorting on the category
def assign_category_number(cat):
    if cat == "Proliferation":
        return 1
    elif cat  == "Maturation":
        return 2
    elif cat == "Mineralization":
        return 3
    else:
        return 4
tpm_marker_long['Cat_Sort'] = tpm_marker_long['Category'].apply(assign_category_number)
#Execute sort
tpm_marker_long.sort_values(['Cat_Sort', 'Gene', 'Experiment'], inplace=True)

#Create the faceted bar plot
plt.figure(figsize=(8, 6))
g = sns.FacetGrid(tpm_marker_long, col='Gene', hue='Category', col_wrap=3, sharey=False, height=4, aspect=1.5, legend_out=False)
g.map(sns.barplot, 'Experiment', 'TPM', order=tpm_marker_long['Experiment'].unique(), ci=None)

for ax in g.axes[6:]:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

# Adjusting the layout for readability
g.set_axis_labels("", "TPM")
g.set_titles("{col_name}")

# Save the plot to the specified location
plt.savefig(os.path.join(qc_dir, 'aggr', 'marker_genes_faceted_barplot.png'), dpi=300, bbox_inches='tight')
plt.clf()

