import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from scipy.stats import chi2_contingency, mannwhitneyu
from skbio.diversity import beta_diversity
from skbio.diversity.alpha import observed_otus, shannon
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from statsmodels.stats.multitest import multipletests

def save_fig(fig_, savepath):
    fig_.update_layout(dragmode='pan', margin=dict(l=30, r=30, t=30, b=30))
    fig_.write_html(savepath, config={'scrollZoom': True, 'displaylogo': False})


pio.templates.default = 'plotly_white'

# %%
otutable = pd.read_csv('data/processed/otutable.csv', index_col=0)
taxtable = pd.read_csv('data/processed/taxtable.csv', index_col=0)
taxtable_genus = pd.read_csv('data/processed/taxtable_genus.csv', index_col=0)
otutable_genus = pd.read_csv('data/processed/otutable_genus.csv', index_col=0)
metadata = pd.read_csv('data/phylo_sample-metadata.txt', sep='\t', index_col=0)


# %% md
# In the dataset, the samples differ by 'sampling_site', 'rhizosphere', and 'season'.
#
# At first, we'll check independency of the features.
#
# Later, we'll access alpha diversity, beta diversity, relative abundance of phylum and genera,
# differential abundance, depending on these factors and
# perform functional profiling of differing and abundant phyla and genera
# via literature, trying to explain seen differences.

# %% md
# # Checking features independency

# %%
relevant_metadata = metadata.reindex(otutable.columns)

# %%
# display sample groups
print(json.dumps({str(k): list(v) for k, v in
                  relevant_metadata.groupby(
    ['sampling_site', 'season', 'rhizosphere']
).groups.items()}, indent=2))


# %%
# checking independence of features in metadata
contingency_table_rhizosphere_season = pd.crosstab(
    relevant_metadata['rhizosphere'], relevant_metadata['season'])
contingency_table_rhizosphere_sampling_site = pd.crosstab(
    relevant_metadata['rhizosphere'], relevant_metadata['sampling_site'])
contingency_table_season_sampling_site = pd.crosstab(
    relevant_metadata['season'], relevant_metadata['sampling_site'])

# Chi-Square Tests
chi2_rhizosphere_season, p_rhizosphere_season, _, _ = chi2_contingency(
    contingency_table_rhizosphere_season)
chi2_rhizosphere_sampling_site, p_rhizosphere_sampling_site, _, _ = chi2_contingency(
    contingency_table_rhizosphere_sampling_site)
chi2_season_sampling_site, p_season_sampling_site, _, _ = chi2_contingency(
    contingency_table_season_sampling_site)

# Print p-values
print(f"p-value (rhizosphere vs season): {p_rhizosphere_season}")
print(f"p-value (rhizosphere vs sampling_site): {p_rhizosphere_sampling_site}")
print(f"p-value (season vs sampling_site): {p_season_sampling_site}")


# %% md
# - 'rhizosphere' and 'season' are dependent (lower p-value of 0.0001).
# - 'rhizosphere' and 'sampling_site', as well as 'season' and 'sampling_site',
# are not significantly dependent (higher p-values).
#
# Dependencies between 'rhizosphere' and 'season' imply that the microbial community structure
# varies significantly between different seasons within each type of rhizosphere environment.
#
# Therefore, we will perform analysis of rhizosphere separately in each season.
# And we will analyse microbiomial diversity by sampling_sites separately
# (only N1 and N2).

# %% md
# # Beta Diversity

# %%
otutable = otutable.loc[(otutable != 0).any(axis=1)]

hellinger_transformed_data = np.sqrt(otutable.T.apply(lambda x: x / x.sum(), axis=0))


distance_matrix = beta_diversity(
    'braycurtis', hellinger_transformed_data, ids=hellinger_transformed_data.index)

ordination_results = pcoa(distance_matrix)

explained_variances = ordination_results.proportion_explained

perc_pc1 = explained_variances.iloc[0] * 100
perc_pc2 = explained_variances.iloc[1] * 100

# %%
fig = px.scatter(ordination_results.samples.join(metadata),
                 x='PC1', y='PC2',
                 color='rhizosphere',
                 symbol='season',
                 text='sampling_site',
                 title='PCoA of Beta Diversity in All Samples (Bray-Curtis dissimilarity)')
fig.update_traces(textposition='bottom center')
fig.update_layout(
    xaxis_title=f'PC1 ({perc_pc1:.2f}%)',
    yaxis_title=f'PC2 ({perc_pc2:.2f}%)',
)
fig.show()

# %%
figures_dir = Path('results/figures/')
save_fig(fig, figures_dir / 'PCoA_beta_diversity.html')

# %% md
# Site T3 shows the most distinct separation between rhizosphere and soil samples,
# especially in the Fall.
# Sites N1 and N2 have more overlapping rhizosphere and soil samples, indicating less distinct
# microbial community differences between these environments.
# Seasonal Differences: There is some seasonal variation, with Fall samples showing more distinct
# clustering than Spring samples.


# %%
result = permanova(distance_matrix, metadata, column='rhizosphere', permutations=999)
print(result)


# %% md
# PERMANOVA shows statistically significant differences between rhizosphere and free soil samples.


# %%
dm_df = pd.DataFrame(distance_matrix.data, index=distance_matrix.ids, columns=distance_matrix.ids)

within_rhizosphere = []
within_bulk_soil = []

# Loop through the distance matrix to extract within-group distances
for i in range(len(dm_df)):
    for j in range(i + 1, len(dm_df)):
        sample_i = dm_df.index[i]
        sample_j = dm_df.columns[j]
        group_i = metadata.loc[sample_i, 'rhizosphere']
        group_j = metadata.loc[sample_j, 'rhizosphere']

        if group_i == group_j:
            if group_i == 'rhizosphere':
                within_rhizosphere.append(dm_df.iloc[i, j])
            elif group_i == 'soil':
                within_bulk_soil.append(dm_df.iloc[i, j])

# Combine the lists into a DataFrame
dissimilarity_df = pd.DataFrame({
    'Dissimilarity': within_rhizosphere + within_bulk_soil,
    'Compartment': ['Rhizosphere'] * len(within_rhizosphere) + ['Bulk Soil'] * len(within_bulk_soil)
})

# Create the violin plot using Plotly
fig = px.violin(dissimilarity_df, x='Compartment', y='Dissimilarity', box=True,
                color='Compartment', points='outliers',
                title='Within-Group Dissimilarity in Rhizosphere and Bulk Soil')
fig.show()
save_fig(fig, figures_dir / 'braycurtis_dissimilarity.html')

# %%
stat, p_value = mannwhitneyu(within_rhizosphere, within_bulk_soil, alternative='two-sided')

print(f"Wilcoxon rank-sum test statistic: {stat}")
print(f"p-value: {p_value}")


# %% md
# Dissimilarity within groups (i.e. differences between different communities within rhizosphere
# and free soil) showed no statistical differences.


# %% md
# # Alpha Diversity


# %% md
# ## Total

# %%
# Calculate observed species richness and Shannon's diversity index
alpha_diversity = pd.DataFrame(index=otutable.T.index)
alpha_diversity['observed_otus'] = otutable.T.apply(observed_otus, axis=1)
alpha_diversity['shannon'] = otutable.T.apply(shannon, axis=1)

# Merge alpha diversity metrics with metadata
alpha_diversity = alpha_diversity.join(metadata)

# %%
fig = px.violin(
    alpha_diversity,
    x='rhizosphere',
    y='observed_otus',
    color='rhizosphere',
    # facet_col='season',
    box=True,
    points='all',
    title="Observed Species Richness Total",
    labels={'observed_otus': "Observed Species Richness", 'sampling_site': 'Sampling Site'}
)
fig.show()
save_fig(fig, figures_dir / 'observed_species_richness_all.html')

# %%
stat, p_value = mannwhitneyu(
    alpha_diversity[alpha_diversity['rhizosphere'] == 'rhizosphere']['observed_otus'],
    alpha_diversity[alpha_diversity['rhizosphere'] == 'soil']['observed_otus'],
    alternative='two-sided')

print(f"Wilcoxon rank-sum test statistic: {stat}")
print(f"p-value: {p_value}")

# %%
fig = px.violin(
    alpha_diversity,
    x='rhizosphere',
    y='shannon',
    color='rhizosphere',
    # facet_col='season',
    box=True,
    points='all',
    title="Shannon's Diversity Index Total",
    labels={'shannon': "Shannon's Index", 'sampling_site': 'Sampling Site'},
)

fig.show()
save_fig(fig, figures_dir / 'shannon_index_all.html')


# %%
stat, p_value = mannwhitneyu(
    alpha_diversity[alpha_diversity['rhizosphere'] == 'rhizosphere']['shannon'],
    alpha_diversity[alpha_diversity['rhizosphere'] == 'soil']['shannon'],
    alternative='two-sided')

print(f"Wilcoxon rank-sum test statistic: {stat}")
print(f"p-value: {p_value}")

# %% md
# Alpha diversity within rhizosphere and free soil didn't show any statistical differences
# with all samples combined.


# %% md
# ## Per sampling site

# %%
sampling_sites = alpha_diversity['sampling_site'].unique()


# Number of comparisons (sampling sites)
num_comparisons = len(sampling_sites)

results = []

for site in sampling_sites:
    # Subset data for the current sampling site
    subset = alpha_diversity[alpha_diversity['sampling_site'] == site]

    # Perform Wilcoxon signed-rank test for Shannon diversity index
    stat, pval = mannwhitneyu(subset[subset['rhizosphere'] == 'rhizosphere']['shannon'],
                              subset[subset['rhizosphere'] == 'soil']['shannon'])

    # Store the results
    results.append({
        'Sampling Site': site,
        'Statistic': stat,
        'P-value': pval
    })

# Convert results to a DataFrame for easier interpretation
results_df = pd.DataFrame(results)

# Apply Bonferroni correction
reject, corrected_pvals, _, _ = multipletests(results_df['P-value'], method='bonferroni')

# Add corrected p-values to results DataFrame
results_df['Corrected P-value'] = corrected_pvals
results_df['Reject Null Hypothesis'] = reject

# Print or visualize the corrected results
results_df

# %%
fig = px.violin(
    alpha_diversity,
    x='sampling_site',
    y='shannon',
    color='rhizosphere',
    # facet_col='season',
    box=True,
    points='all',
    title="Shannon's Diversity Index by Sampling Site",
    labels={'shannon': "Shannon's Index", 'sampling_site': 'Sampling Site'}
)

fig.add_trace(go.Scatter(
    x=['N1', 'T1', 'T3', 'N2'],
    y=[7.25, 7.25, 7.25, 7.25],
    mode="text",
    text=["*", "", "*", "*"],
    textposition="top center",
    textfont=dict(size=16),
    showlegend=False
))

fig.show()
save_fig(fig, figures_dir / 'shannon_by_sample_site.html')


# %%
sampling_sites = alpha_diversity['sampling_site'].unique()


# Number of comparisons (sampling sites)
num_comparisons = len(sampling_sites)

results = []

for site in sampling_sites:
    # Subset data for the current sampling site
    subset = alpha_diversity[alpha_diversity['sampling_site'] == site]

    # Perform Wilcoxon signed-rank test for Shannon diversity index
    stat, pval = mannwhitneyu(subset[subset['rhizosphere'] == 'rhizosphere']['observed_otus'],
                              subset[subset['rhizosphere'] == 'soil']['observed_otus'])

    # Store the results
    results.append({
        'Sampling Site': site,
        'Statistic': stat,
        'P-value': pval
    })

# Convert results to a DataFrame for easier interpretation
results_df = pd.DataFrame(results)

# Apply Bonferroni correction
reject, corrected_pvals, _, _ = multipletests(results_df['P-value'], method='bonferroni')

# Add corrected p-values to results DataFrame
results_df['Corrected P-value'] = corrected_pvals
results_df['Reject Null Hypothesis'] = reject

# Print or visualize the corrected results
results_df


# %%
fig = px.violin(
    alpha_diversity,
    x='sampling_site',
    y='observed_otus',
    color='rhizosphere',
    # facet_col='season',
    box=True,
    points='all',
    title="Observed Species Richness by Sampling Site",
    labels={'observed_otus': "Observed Species Richness", 'sampling_site': 'Sampling Site'}
)
fig.add_trace(go.Scatter(
    x=['N1', 'T1', 'T3', 'N2'],
    y=[2100, 2100, 2100, 2100],
    mode="text",
    text=["*", "", "*", "*"],
    textposition="top center",
    textfont=dict(size=16),
    showlegend=False
))

fig.show()
save_fig(fig, figures_dir / 'observed_species_richness_by_sample_site.html')

# %% md
# Alpha diversity was higher in free soil compared to rhizosphere in samples **N1** and **N2**,
# which are taiga and transitional ecotonic forest, respectively.
#
# Alpha diversity between rhizosphere and soil haven't shown statistical
# differences in **T1** (taiga).
#
# And in **T3**, which is a comparatively oligotrophic sampling site, alpha diversity of
# rhizosphere was significantly higher, than in free soil.
