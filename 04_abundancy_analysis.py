from collections import Counter

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from faprotax_search import faprotax_search

pio.templates.default = 'plotly_white'


def save_fig(fig_, savepath):
    fig_.update_layout(dragmode='pan', margin=dict(l=30, r=30, t=30, b=30))
    fig_.write_html(savepath, config={'scrollZoom': True, 'displaylogo': False})

# %% md
# # Differential abundance analysis


# %% md
# This notebook/script plots the results of differential abundance analysis
# and performs functional annotation of abundant genera via FAPROTAX database.

# %% md
# ## Total
# %%
daa = pd.read_csv('results/DAA/DAA_rhizosphere.csv', index_col=0)
taxtable = pd.read_csv('data/processed/taxtable.csv', index_col=0)


# %%
merged_data = pd.merge(daa, taxtable[['Phylum', 'Genus']],
                       left_index=True, right_index=True, how='left')
merged_data = merged_data[merged_data['pvalue'] < 0.05]

# %%
# Filter out rows with NaN in Phylum or Genus columns
merged_data = merged_data.dropna(subset=['Phylum'])

# %%
top_phyla_positive = merged_data.groupby(
    'Phylum')['log2FoldChange'].mean().nlargest(5).reset_index()
top_phyla_negative = merged_data.groupby(
    'Phylum')['log2FoldChange'].mean().nsmallest(15).reset_index()

# %%
phyla_df = pd.concat([top_phyla_positive, top_phyla_negative]).sort_values('log2FoldChange')
phyla_df['color'] = phyla_df['log2FoldChange'] < 0
fig = px.bar(phyla_df, y='Phylum', x='log2FoldChange', color='color',
             title='Rhizosphere vs soil, by all samples')
fig.update_layout(showlegend=False)
fig.show()
save_fig(fig, 'results/figures/daa_phyla.html')


# %%
top_genera_positive = merged_data.groupby(
    'Genus')['log2FoldChange'].mean().nlargest(20).reset_index()
top_genera_negative = merged_data.groupby(
    'Genus')['log2FoldChange'].mean().nsmallest(20).reset_index()
# %%
genus_df = pd.concat([top_genera_positive, top_genera_negative]).sort_values('log2FoldChange')
genus_df['color'] = genus_df['log2FoldChange'] < 0
fig = px.bar(genus_df, y='Genus', x='log2FoldChange', color='color',
             title='Rhizosphere vs soil, by all samples')
fig.update_layout(showlegend=False)
fig.show()
save_fig(fig, 'results/figures/daa_genera.html')

# %% md
# ## By sampling site


# %%
# for each of the most differentially abundant phyla and genera in all sampling sites combined
# show, in which sampling sites is it also among top differentially abundant taxons and are they
# differentially abundant in the same direction (corr_phyla, corr_genera), or the opposite
# (opp_phyla, opp_genera)

corresponding_phyla = []
corresponding_phyla_samples = []
opposite_phyla = []
opposite_phyla_samples = []
corresponding_genera = []
corresponding_genera_samples = []
opposite_genera = []
opposite_genera_samples = []


for sampling_site in ['N1', 'T1', 'N2', 'T3']:
    daa = pd.read_csv(f'results/DAA/DAA_rhizosphere_{sampling_site}.csv', index_col=0)
    taxtable = pd.read_csv('data/processed/taxtable.csv', index_col=0)
    merged_data = pd.merge(daa, taxtable[['Phylum', 'Genus']],
                           left_index=True, right_index=True, how='left')
    merged_data = merged_data[merged_data['pvalue'] < 0.05]
    merged_data = merged_data.dropna(subset=['Phylum'])
    top_phyla_positive_ss = merged_data.groupby(
        'Phylum')['log2FoldChange'].mean().nlargest(10).reset_index()
    top_phyla_negative_ss = merged_data.groupby(
        'Phylum')['log2FoldChange'].mean().nsmallest(10).reset_index()

    top_genera_positive_ss = merged_data.groupby(
        'Genus')['log2FoldChange'].mean().nlargest(20).reset_index()
    top_genera_negative_ss = merged_data.groupby(
        'Genus')['log2FoldChange'].mean().nsmallest(20).reset_index()

    for phylum in top_phyla_positive['Phylum']:
        if phylum in top_phyla_positive_ss['Phylum'].values:
            corresponding_phyla.append(phylum)
            corresponding_phyla_samples.append(sampling_site)
        elif phylum in top_phyla_negative_ss['Phylum'].values:
            opposite_phyla.append(phylum)
            opposite_phyla_samples.append(sampling_site)

    for phylum in top_phyla_negative['Phylum']:
        if phylum in top_phyla_positive_ss['Phylum'].values:
            opposite_phyla.append(phylum)
            opposite_phyla_samples.append(sampling_site)
        elif phylum in top_phyla_negative_ss['Phylum'].values:
            corresponding_phyla.append(phylum)
            corresponding_phyla_samples.append(sampling_site)

    for genus in top_genera_positive['Genus']:
        if genus in top_genera_positive_ss['Genus'].values:
            corresponding_genera.append(genus)
            corresponding_genera_samples.append(sampling_site)
        elif genus in top_genera_negative_ss['Genus'].values:
            opposite_genera.append(genus)
            opposite_genera_samples.append(sampling_site)

    for genus in top_genera_negative['Genus']:
        if genus in top_genera_positive_ss['Genus'].values:
            opposite_genera.append(genus)
            opposite_genera_samples.append(sampling_site)
        elif genus in top_genera_negative_ss['Genus'].values:
            corresponding_genera.append(genus)
            corresponding_genera_samples.append(sampling_site)


# %%
# Organize results into a dataframe
corr_phyla = pd.DataFrame({
    'Corresponding_Phyla': corresponding_phyla,
    'Corresponding_Phyla_Samples': corresponding_phyla_samples, }
).groupby(
    'Corresponding_Phyla').apply(
        lambda x: ','.join(x['Corresponding_Phyla_Samples']),
        include_groups=False).reset_index()

opp_phyla = pd.DataFrame({
    'Opposite_Phyla': opposite_phyla,
    'Opposite_Phyla_Samples': opposite_phyla_samples, }).groupby(
        'Opposite_Phyla').apply(
            lambda x: ','.join(x['Opposite_Phyla_Samples']),
            include_groups=False).reset_index()

corr_genera = pd.DataFrame({
    'Corresponding_Genera': corresponding_genera,
    'Corresponding_Genera_Samples': corresponding_genera_samples, }).groupby(
        'Corresponding_Genera').apply(
            lambda x: ','.join(x['Corresponding_Genera_Samples']),
            include_groups=False).reset_index()

opp_genera = pd.DataFrame({
    'Opposite_Genera': opposite_genera,
    'Opposite_Genera_Samples': opposite_genera_samples
}).groupby(
    'Opposite_Genera').apply(
        lambda x: ','.join(x['Opposite_Genera_Samples']),
        include_groups=False).reset_index()


# %% md
# # Relative abundance analysis

# %%
otutable = pd.read_csv('data/processed/otutable.csv', index_col=0)
taxtable = pd.read_csv('data/processed/taxtable.csv', index_col=0)
metadata = pd.read_csv('data/phylo_sample-metadata.txt', index_col=0, sep='\t')
metadata = metadata.reindex(otutable.columns)


# %%
phylum_composition = (otutable
                      .join(taxtable['Phylum'])
                      .groupby(by='Phylum').sum())

# %%
phylum_composition_rel = phylum_composition.div(
    phylum_composition.sum(axis=0), axis=1)


# %%
plot_df_total = phylum_composition_rel.T.join(metadata[['rhizosphere']]
                                              ).groupby('rhizosphere').mean().T.sort_values(
    'rhizosphere', ascending=True).reset_index().rename(
        columns={'index': 'Phylum'}
)


# %%
df_melted = plot_df_total.melt(id_vars='Phylum', var_name='Sample',
                               value_name='Relative Abundance')

# %%
df_melted['Sample'] = pd.Categorical(df_melted['Sample'], categories=[
                                     'rhizosphere', 'soil'], ordered=True)

colors = [
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692"
]

phylum_color_map = {phylum: colors[i % len(colors)] for i, phylum in
                    enumerate(plot_df_total['Phylum'].unique())}

df_melted['color'] = df_melted['Phylum'].map(phylum_color_map)

fig = go.Figure()
for phylum in plot_df_total['Phylum'].unique():
    df_phylum = df_melted[df_melted['Phylum'] == phylum]
    fig.add_trace(go.Bar(
        x=df_phylum['Sample'],
        y=df_phylum['Relative Abundance'],
        text=phylum,
        textposition='inside',
        name=phylum,
        marker_color=phylum_color_map[phylum],
        hoverinfo='y+name',
    ))

fig.update_layout(
    title='Stacked Bar Plot of Phyla in Rhizosphere vs. Soil',
    xaxis_title='Sample Type',
    yaxis_title='Relative Abundance',
    barmode='stack'
)

fig.show()
save_fig(fig, 'results/figures/relative_abundancy_phyla.html')


# %%
genus_composition = (otutable
                     .join(taxtable['Genus']).dropna(subset='Genus')
                     .groupby(by='Genus').sum())

# %%
genus_composition_rel = genus_composition.div(
    genus_composition.sum(axis=0), axis=1)


# %%
plot_df_total = genus_composition_rel.T.join(metadata[['rhizosphere']]
                                             ).groupby('rhizosphere').mean().T.sort_values(
    'rhizosphere', ascending=True).reset_index().rename(
        columns={'index': 'Genus'}
)


# %%
df_melted = plot_df_total.melt(id_vars='Genus', var_name='Sample',
                               value_name='Relative Abundance')

# %%
df_melted['Sample'] = pd.Categorical(df_melted['Sample'], categories=[
                                     'rhizosphere', 'soil'], ordered=True)

colors = [
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692"
]

genus_color_map = {phylum: colors[i % len(colors)] for i, phylum in
                   enumerate(plot_df_total['Genus'].unique())}

df_melted['color'] = df_melted['Genus'].map(genus_color_map)

fig = go.Figure()

for genus in plot_df_total['Genus'].unique():
    df_genus = df_melted[df_melted['Genus'] == genus]
    fig.add_trace(go.Bar(
        x=df_genus['Sample'],
        y=df_genus['Relative Abundance'],
        text=genus,
        textposition='inside',
        name=genus,
        marker_color=genus_color_map[genus],
        hoverinfo='y+name',
    ))

fig.update_layout(
    title='Stacked Bar Plot of Genera in Rhizosphere vs. Soil',
    xaxis_title='Sample Type',
    yaxis_title='Relative Abundance',
    barmode='stack',
    showlegend=True,
)
fig.show()
save_fig(fig, 'results/figures/relative_abundancy_genera.html')


# %% md
# # Functional annotation

# %%
up_genus = genus_df[genus_df['log2FoldChange'] > 0]['Genus'].tolist()
down_genus = genus_df[genus_df['log2FoldChange'] < 0]['Genus'].tolist()

# %%
up_genus_functions = Counter([i for gen in up_genus for i in faprotax_search(gen)
                              if not '*' in i])
up_genus_functions

# %%
down_genus_functions = Counter([i for gen in down_genus for i in faprotax_search(gen)
                                if not '*' in i])
down_genus_functions


# %% md
# In the **rhizosphere**, there are more taxa involved in:
#
# - Nitrogen fixation
# - Plant pathogenesis
# - Ligninolysis, cellulolysis, xylanolysis
# - Methanol oxidation
#
# In **free soil**, there are more taxa involved in:
#
# - Ammonia and nitrite oxidation
# - Various sulfur respiration processes
# - Iron respiration and iron oxidation
