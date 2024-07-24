# Soil Metagenomics Analysis

This repository includes all the code created during the 2-day hackathon [Bioinformatics Bootcamp 2024](https://pish.itmo.ru/genomics-bootcamp) organized by ITMO University (Saint-Petersburg). The data provided was an OTU table of 16s rRNA-seq (`phylo_otutable.csv`), a table of ASV-taxonomy (`phylo_taxtable.csv`), and metadata (`phylo_sample-metadata.txt`) containing information about the sampling site, rhizosphere or free soil, and season of metagenomics from sampling sites: **N1**, **T1** — taiga, **N2** — transitional ecotonic forest, **T3** — relatively oligotrophic environment, described in [(Polyakov et al., 2022)](https://www.mdpi.com/2073-4395/12/8/1760).


## Installation

To install the required packages, run:

```bash
pip install -r requirements.txt
```

## Directory Structure
- `results/`
  - Contains results of the analysis, including figures and differential abundance analysis (DAA) results.
- `data/`
  - Contains raw and processed data.
- `FAPROTAX_1.2.10`
  - [FAPROTAX database](http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Home)

## Notebooks and Scripts
- `01_preprocessing_and_dimension_reduction.ipynb` / `01_preprocessing_and_dimension_reduction.py`:
  - **Purpose**: Preprocesses the soil metagenomics dataset and performs dimension reduction.
  - **Steps**:
    - Loads and merges data tables.
    - Filters the dataset based on taxonomic levels and abundance.
    - Checks the impact of filtering on alpha diversity and taxonomic composition.
    - Performs PCA and UMAP for dimension reduction and visualizes the results.
- `02_diversity.ipynb` / `02_diversity.py`:
  - **Purpose**: Analyzes alpha and beta diversity of the soil metagenomics samples.
  - **Steps**:
    - Checks the independence of metadata features.
    - Calculates and visualizes beta diversity using PCoA.
    - Performs PERMANOVA to test for significant differences.
    - Analyzes within-group dissimilarity.
    - Calculates and visualizes alpha diversity (observed species richness and Shannon's diversity index).
- `03_DAA.ipynb` / `03_DAA.py`:
  - **Purpose**: Performs differential abundance analysis (DAA) using DESeq2.
  - **Steps**:
    - Loads the processed OTU table and metadata.
    - Runs DESeq2 to identify differentially abundant ASVs between rhizosphere and soil samples.
    - Saves the DAA results for all samples and for each sampling site.
- `04_abundancy_analysis.ipynb` / `04_abundancy_analysis.py`:
  - **Purpose**: Analyzes the relative abundance of phyla and genera and performs functional annotation.
  - **Steps**:
    - Loads the DAA results and taxonomic table.
    - Identifies top differentially abundant phyla and genera.
    - Visualizes the relative abundance of phyla and genera in rhizosphere vs. soil.
    - Performs functional annotation of abundant genera using the FAPROTAX database.
- `faprotax_search.py`:
  - **Purpose**: Provides a function (`faprotax_search(search_term, file_path)`) to search the FAPROTAX database for functional annotations of taxa.

## Conclusions

The beta diversity analysis showed that the microbial composition of the rhizosphere differs significantly from that of free soil.

In oligotrophic biomes, the rhizosphere microbiome is more diverse compared to the free soil microbiome, whereas in nutrient-rich soils, the opposite trend is observed.

In the rhizosphere, taxa involved in methanol oxidation, nitrogen fixation, plant pathogenesis, and the decomposition of plant-specific substances are more common.
