from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.io as pio
from scipy import stats
from skbio.diversity import alpha_diversity
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from umap import UMAP

def save_fig(fig_, savepath):
    fig_.update_layout(dragmode='pan', margin=dict(l=30, r=30, t=30, b=30))
    fig_.write_html(savepath, config={'scrollZoom': True, 'displaylogo': False})


pio.templates.default = 'plotly_white'

# %% codecell
# load data
tax_table = pd.read_csv('data/phylo_taxtable.csv', index_col=0)
otutable = pd.read_csv('data/phylo_otutable.csv', index_col=0)
metadata = pd.read_csv('data/phylo_sample-metadata.txt', sep='\t', index_col=0)

assert (tax_table.index == otutable.index).all(), "ASV in taxonomy and otutable don't match!"
metadata = metadata.loc[otutable.columns]
assert all(otutable.columns ==
           metadata.index), "Columns in otutable and index in metadata do not match!"

# %% markdown
# # Preprocessing

# %% codecell
# add column with all taxon levels concatenated
def concatenate_levels(row):
    return '-'.join(filter(pd.notna, row))


tax_table['full'] = tax_table.apply(concatenate_levels, axis=1)

# %% codecell
# create tables with only known phylum
tax_table_phylum = tax_table.dropna(subset='Phylum')
otutable_phylum = otutable.reindex(tax_table_phylum.index)
assert (tax_table_phylum.index == otutable_phylum.index).all(
), "ASV in taxonomy and otutable don't match!"

# %% codecell
# create tables with only known genus
tax_table_genus = tax_table.dropna(subset='Genus')
otutable_genus = otutable.reindex(tax_table_genus.index)
assert (tax_table_genus.index == otutable_genus.index).all(
), "ASV in taxonomy and otutable don't match!"

# %% codecell
mean_asv = otutable_phylum.mean(axis=1)

# plot distribution of mean counts of ASVs
fig = ff.create_distplot([mean_asv[mean_asv < 100]], group_labels=['ASV presence'])
fig.update_layout(title='Average ASV presence (known Phylum)',
                  xaxis_title='Mean of Counts',
                  yaxis_title='Density',
                  showlegend=False)
fig.show()
# %%
figures_dir = Path('results/figures')
figures_dir.mkdir(parents=True, exist_ok=True)
save_fig(fig, figures_dir / 'average_asv_presence_phylum_all.html')

# %% md
# ## Removing low-abundance ASVs and singletons


# %% md
# Removing low-abundance ASVs as they might be associated with sequencing errors
# and our study is not focused on detecting rare species or understanding low-abundance dynamics.
# Also removing singletons (ASVs appearing in only one sample)

# %% codecell
otutable_phylum_filtered = otutable_phylum[mean_asv > 0.8]
otutable_phylum_filtered = otutable_phylum_filtered[(otutable_phylum_filtered > 0).sum(axis=1) > 1]

# %% codecell
print('ASVs with known phylum before filtering: ', len(otutable_phylum))
print('ASVs with known phylum after filtering: ', len(otutable_phylum_filtered))

# %% codecell
# plot distributions of mean counts of ASVs (filtered)
mean_asv = otutable_phylum_filtered.mean(axis=1)
fig = ff.create_distplot([mean_asv[mean_asv < 100]], group_labels=['ASV presence'])
fig.update_layout(title='Average ASV presence (known Phylum)',
                  xaxis_title='Mean of Counts',
                  yaxis_title='Density',
                  showlegend=False)
fig.show()
# %%
save_fig(fig, figures_dir / 'average_asv_presence_filtered.html')

# %% md
# ### Checking that filtering the dataset didn't affect alpha-diversity and composition on phylum level

# %% codecell
# check that alpha diversity is not affected after filtering
shannon_original = alpha_diversity('shannon', otutable_phylum.T)
shannon_filtered = alpha_diversity('shannon', otutable_phylum_filtered.T)

print("Original Shannon Diversity:")
print(shannon_original.describe())
print("\nFiltered Shannon Diversity:")
print(shannon_filtered.describe())


# %% codecell
# Kolmogorov-Smirnov test checks whether the distributions before and after filtering differ
stats.kstest(shannon_original, shannon_filtered)

# %% md
# p-value > 0.05, meaning the distribution changed insignificantly

# %% codecell
# check that taxonomic composition is not affected after filtering
phylum_composition_original = (otutable_phylum
                               .join(tax_table_phylum['Phylum'])
                               .groupby(by='Phylum').sum())
phylum_composition_filtered = (otutable_phylum_filtered
                               .join(tax_table_phylum['Phylum'])
                               .groupby(by='Phylum').sum())

phylum_composition_original_rel = phylum_composition_original.div(
    phylum_composition_original.sum(axis=0), axis=1)
phylum_composition_filtered_rel = phylum_composition_filtered.div(
    phylum_composition_filtered.sum(axis=0), axis=1)


# %%
def prepare_plotting_df(composition_df, label, id_vars):
    composition_df = composition_df.reset_index()
    composition_df_melted = composition_df.melt(
        id_vars=id_vars, var_name='Sample', value_name='Abundance')
    composition_df_melted['Dataset'] = label
    return composition_df_melted


# %%
df_plot_original = prepare_plotting_df(phylum_composition_original_rel, 'Original', 'Phylum')
df_plot_filtered = prepare_plotting_df(phylum_composition_filtered_rel, 'Filtered', 'Phylum')

df_plot_combined = pd.concat([df_plot_original, df_plot_filtered])

# %% codecell
fig = px.bar(df_plot_combined, x='Sample', y='Abundance', color='Phylum', facet_col='Dataset', barmode='stack',
             title="Taxonomic Composition Before and After Filtering",
             labels={'Abundance': 'Relative Abundance', 'Sample': 'Sample'})

fig.update_layout(title_text='Taxonomic Composition Changes',
                  yaxis_range=[0, 1])
fig.show()
# %%
save_fig(fig, figures_dir / 'taxonomic_comp_phylum.html')

# %% markdown
# Taxonomic composition is not affected

# %% md
# ### Checking that filtering the dataset didn't affect alpha-diversity and composition on genus level

# %% codecell
# create genus tables from filtered data
tax_table_phylum_filtered = tax_table_phylum.reindex(otutable_phylum_filtered.index)
tax_table_genus_filtered = tax_table_phylum_filtered.dropna(subset='Genus')
otutable_genus_filtered = otutable_phylum_filtered.reindex(tax_table_genus_filtered.index)
assert (tax_table_genus.index == otutable_genus.index).all(
), "ASV in taxonomy and otutable don't match!"

# %% codecell
print('ASVs with known genus before filtering: ', len(otutable_genus))
print('ASVs with known genus after filtering: ', len(otutable_genus_filtered))


# %% codecell
# check that alpha diversity is not affected after filtering
shannon_original = alpha_diversity('shannon', otutable_genus.T)
shannon_filtered = alpha_diversity('shannon', otutable_genus_filtered.T)

print("Original Shannon Diversity:")
print(shannon_original.describe())
print("\nFiltered Shannon Diversity:")
print(shannon_filtered.describe())


# %% codecell
# Kolmogorov-Smirnov test checks whether the distributions before and after filtering differ
stats.kstest(shannon_original, shannon_filtered)


# %% md
# p-value > 0.05, meaning the distribution changed insignificantly

# %% codecell
# check that taxonomic composition is not affected after filtering
genus_composition_original = (otutable_genus
                              .join(tax_table_genus['Genus'])
                              .groupby(by='Genus').sum())
genus_composition_filtered = (otutable_genus_filtered
                              .join(tax_table_genus['Genus'])
                              .groupby(by='Genus').sum())

genus_composition_original_rel = genus_composition_original.div(
    genus_composition_original.sum(axis=0), axis=1)
genus_composition_filtered_rel = genus_composition_filtered.div(
    genus_composition_filtered.sum(axis=0), axis=1)

df_plot_original = prepare_plotting_df(genus_composition_original_rel, 'Original', 'Genus')
df_plot_filtered = prepare_plotting_df(genus_composition_filtered_rel, 'Filtered', 'Genus')

df_plot_combined = pd.concat([df_plot_original, df_plot_filtered])
# %% codecell
fig = px.bar(df_plot_combined, x='Sample', y='Abundance', color='Genus', facet_col='Dataset', barmode='stack',
             title="Taxonomic Composition Before and After Filtering",
             labels={'Abundance': 'Relative Abundance', 'Sample': 'Sample'})

fig.update_layout(title_text='Taxonomic Composition Changes',
                  yaxis_range=[0, 1])
fig.show()

# %%
save_fig(fig, figures_dir / 'taxonomic_comp_genus.html')

# %% markdown
# Taxonomic composition is not affected


# %% md
# # Save filtered data

# %% codecell
out_dir = Path('data/processed')
out_dir.mkdir(parents=True, exist_ok=True)

otutable_phylum_filtered.to_csv(out_dir / 'otutable.csv')
tax_table_phylum.to_csv(out_dir / 'taxtable.csv')
tax_table_genus_filtered.to_csv(out_dir / 'taxtable_genus.csv')
otutable_genus_filtered.to_csv(out_dir / 'otutable_genus.csv')

# %% markdown
# # Dimension reduction
# %% codecell
# perform PCA analysis
features_scaled = StandardScaler().fit_transform(otutable_phylum_filtered.T)

pca = PCA()
pca_model = pca.fit_transform(features_scaled)

cluster_values = pd.DataFrame(pca_model, columns=[f'PC{i + 1}' for i in range(pca_model.shape[1])])
cluster_values = cluster_values.join(metadata.reset_index())
# %% codecell
# PCA scree plot
explained_variance = pca.explained_variance_ratio_
cumulative_explained_variance = np.cumsum(explained_variance)

fig = px.line(x=range(1, len(cumulative_explained_variance) + 1),
              y=cumulative_explained_variance, markers=True)
fig.update_layout(
    title='PCA Elbow Plot',
    xaxis_title='Number of Components',
    yaxis_title='Cumulative Explained Variance',
    showlegend=False
)

fig.show()
save_fig(fig, figures_dir / 'pca_scree.html')

# %% codecell
by = 'rhizosphere'  # rhizosphere, soil_type, season
fig = px.scatter(cluster_values, x='PC1', y='PC2',
                 color=by, hover_data='index',
                 title='PCA of Soil Metagenomics Samples',
                 labels={'index': 'sample'})

fig.show()
save_fig(fig, figures_dir / 'pca.html')

# %% markdown
# Only 22% of variance is described by first 2 components (as seen from scree plot).
# Also the data doesn't cluster well.
#
# Below is the UMAP representation, which clusters data better,
# but it may struggle to preserve the balance between global and local structure.
# %% codecell
# perform UMAP analysis
features_scaled = StandardScaler().fit_transform(otutable_phylum_filtered.T)

umap = UMAP(n_components=2, min_dist=.5, n_neighbors=20)
umap_model = umap.fit_transform(features_scaled)

cluster_values = pd.DataFrame(
    umap_model, columns=[f'UMAP{i + 1}' for i in range(umap_model.shape[1])])
cluster_values = cluster_values.join(metadata.reset_index())

# %% codecell
fig = px.scatter(cluster_values, x='UMAP1', y='UMAP2',
                 color='rhizosphere', symbol='soil_type', hover_data='index',
                 title='UMAP of Soil Metagenomics Samples',
                 labels={'index': 'sample'})
fig.show()
save_fig(fig, figures_dir / 'umap.html')

# %% markdown
# ## PCA and UMAP of only Chern taiga (sample site N1)

# %% codecell

otutable_taiga = otutable_phylum_filtered.T[(metadata['sampling_site'] == 'N1')
                                            & (metadata['season'] != 'Summer')].T


# %% codecell
# perform PCA analysis
features_scaled = StandardScaler().fit_transform(otutable_taiga.T)

pca = PCA()
pca_model = pca.fit_transform(features_scaled)

cluster_values = pd.DataFrame(pca_model, columns=[f'PC{i + 1}' for i in range(pca_model.shape[1])])
cluster_values = cluster_values.join(metadata.reindex(otutable_taiga.columns).reset_index())
# %% codecell
# PCA scree plot
explained_variance = pca.explained_variance_ratio_
cumulative_explained_variance = np.cumsum(explained_variance)

fig = px.line(x=range(1, len(cumulative_explained_variance) + 1),
              y=cumulative_explained_variance, markers=True)
fig.update_layout(
    title='PCA Elbow Plot',
    xaxis_title='Number of Components',
    yaxis_title='Cumulative Explained Variance',
    showlegend=False
)

fig.show()
save_fig(fig, figures_dir / 'pca_scree_N1.html')
# %% codecell
fig = px.scatter(cluster_values, x='PC1', y='PC2',
                 color='rhizosphere', symbol='season', hover_data='index',
                 title='PCA of Soil Metagenomics Samples',
                 labels={'index': 'sample'})
fig.show()
save_fig(fig, figures_dir / 'pca_N1.html')

# %% codecell
features_scaled = StandardScaler().fit_transform(otutable_taiga.T)

umap = UMAP(n_components=2)
umap_model = umap.fit_transform(features_scaled)

cluster_values = pd.DataFrame(
    umap_model, columns=[f'UMAP{i + 1}' for i in range(umap_model.shape[1])])
cluster_values = cluster_values.join(metadata.reindex(otutable_taiga.columns).reset_index())
# %% codecell
fig = px.scatter(cluster_values, x='UMAP1', y='UMAP2',
                 color='rhizosphere', symbol='season', hover_data='index',
                 title='UMAP of Soil Metagenomics Samples',
                 labels={'index': 'sample'})
fig.show()
save_fig(fig, figures_dir / 'umap_N1.html')

# %% markdown
# ## PCA and UMAP of only control soil -- transitional ecotone forest (sample site N2)

# %% codecell
otutable_taiga = otutable_phylum_filtered.T[metadata['sampling_site'] == 'N2'].T


# %% codecell
# perform PCA analysis
features_scaled = StandardScaler().fit_transform(otutable_taiga.T)

pca = PCA()
pca_model = pca.fit_transform(features_scaled)

cluster_values = pd.DataFrame(pca_model, columns=[f'PC{i + 1}' for i in range(pca_model.shape[1])])
cluster_values = cluster_values.join(metadata.reindex(otutable_taiga.columns).reset_index())

# %% codecell
# PCA scree plot
explained_variance = pca.explained_variance_ratio_
cumulative_explained_variance = np.cumsum(explained_variance)

fig = px.line(x=range(1, len(cumulative_explained_variance) + 1),
              y=cumulative_explained_variance, markers=True)
fig.update_layout(
    title='PCA Elbow Plot',
    xaxis_title='Number of Components',
    yaxis_title='Cumulative Explained Variance',
    showlegend=False
)

fig.show()
save_fig(fig, figures_dir / 'pca_scree_N2.html')
# %% codecell
fig = px.scatter(cluster_values, x='PC1', y='PC2',
                 color='rhizosphere', symbol='season', hover_data='index',
                 title='PCA of Soil Metagenomics Samples',
                 labels={'index': 'sample'})
fig.show()
save_fig(fig, figures_dir / 'pca_N2.html')

# %% codecell
features_scaled = StandardScaler().fit_transform(otutable_taiga.T)

umap = UMAP(n_components=2)
umap_model = umap.fit_transform(features_scaled)

cluster_values = pd.DataFrame(
    umap_model, columns=[f'UMAP{i + 1}' for i in range(umap_model.shape[1])])
cluster_values = cluster_values.join(metadata.reindex(otutable_taiga.columns).reset_index())

# %% codecell
fig = px.scatter(cluster_values, x='UMAP1', y='UMAP2',
                 color='rhizosphere', symbol='season', hover_data='index',
                 title='UMAP of Soil Metagenomics Samples',
                 labels={'index': 'sample'})
fig.show()
save_fig(fig, figures_dir / 'umap_N2.html')

# %% md
# This sampling site shows good clusterization both by soil type (rhizosphere vs free soil) and
# season (fall vs spring)
