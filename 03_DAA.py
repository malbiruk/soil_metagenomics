from pathlib import Path

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


# %% md
# # Differential Abundance Analysis


# %% md
# This notebook/script calculates differential abundance of ASVs for all dataset
# and within each sample.

# %%
otutable = pd.read_csv('data/processed/otutable.csv', index_col=0)
metadata = pd.read_csv('data/phylo_sample-metadata.txt', sep='\t', index_col=0)

metadata = metadata.reindex(otutable.columns)

# Activate pandas2ri conversion
pandas2ri.activate()

# Import DESeq2 package
deseq2 = importr('DESeq2')

# Convert data to R dataframe
otu_table_r = pandas2ri.py2rpy(otutable)
metadata_r = pandas2ri.py2rpy(metadata)


# %%
# Create DESeq2 dataset
dds = deseq2.DESeqDataSetFromMatrix(
    countData=otu_table_r,
    colData=metadata_r,
    design=ro.Formula('~rhizosphere'))


# %%
# Run DESeq2
dds_processed = deseq2.DESeq(dds)


# %%
r_results = deseq2.results(dds_processed, contrast=ro.StrVector(
    ['rhizosphere', 'rhizosphere', 'soil']))
ro.globalenv['r_results'] = r_results
r_results_df = ro.r("as.data.frame(r_results)")

# %%
with localconverter(ro.default_converter + pandas2ri.converter):
    results_df = ro.conversion.rpy2py(r_results_df)

# %%
out_dir = Path('results/DAA')
out_dir.mkdir(parents=True, exist_ok=True)
results_df.to_csv(out_dir / 'DAA_rhizosphere.csv')


# %%
# by samples

def daa_analysis(outtable_, metadata_, outfile):
    pandas2ri.activate()
    deseq2 = importr('DESeq2')
    otu_table_r = pandas2ri.py2rpy(outtable_)
    metadata_r = pandas2ri.py2rpy(metadata_)
    dds = deseq2.DESeqDataSetFromMatrix(
        countData=otu_table_r,
        colData=metadata_r,
        design=ro.Formula('~rhizosphere'))
    dds_processed = deseq2.DESeq(dds)

    r_results = deseq2.results(dds_processed, contrast=ro.StrVector(
        ['rhizosphere', 'rhizosphere', 'soil']))
    ro.globalenv['r_results'] = r_results
    r_results_df = ro.r("as.data.frame(r_results)")

    with localconverter(ro.default_converter + pandas2ri.converter):
        results_df = ro.conversion.rpy2py(r_results_df)

    results_df.to_csv(outfile)


# %%
for sampling_site in metadata.sampling_site.unique():
    sub_otutable = otutable.T.loc[metadata[metadata.sampling_site == sampling_site].index].T
    sub_metadata = metadata.reindex(sub_otutable.columns)
    daa_analysis(sub_otutable, sub_metadata, out_dir / f'DAA_rhizosphere_{sampling_site}.csv')
