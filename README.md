Overview

This repository contains analysis scripts used for integrative topology analysis of overlapping RING1B and H3K27ac HiChIP interaction networks.

The workflow includes:

* construction of HiChIP interaction graphs
* identification of overlapping chromatin interaction clusters
* graph topology analysis
* dimensionality reduction and clustering
* interaction directionality analysis
* integration with RNA-seq differential expression data

The analyses were performed in R.

Repository structure
├── RING1B–H3K27ac HiChIP topology analysis.R
├── RING1B–H3K27ac interaction directionality analysis.R
├── Integration of HiChIP clusters with gene expression changes.R
└── README.md

Analysis workflow

1. RING1B–H3K27ac HiChIP topology analysis.R

This script performs the core topology analysis of RING1B and H3K27ac HiChIP interaction networks.

Main steps include:

1. Build HiChIP interaction graphs
2. Identify connected interaction clusters
3. Compute graph topology features
4. Convert anchors to genomic regions
5. Identify overlapping hub networks
6. Compute structural similarity metrics
7. Perform PCA, UMAP, and K-means clustering

2. RING1B–H3K27ac interaction directionality analysis.R

This script evaluates directional consistency between overlapping RING1B and H3K27ac interaction clusters.

Main steps include:

1. Identify shared anchors
2. Compute anchor-relative genomic distances
3. Calculate anchor-level directionality
4. Assign cluster-level directionality categories
5. Visualize cluster direction distributions

3. Integration of HiChIP clusters with gene expression changes.R

This script integrates HiChIP interaction clusters with RNA-seq differential expression data.

Main steps include:

1. Map RING1B anchors to nearby TSS regions
2. Assign genes to overlapping interaction clusters
3. Integrate RNA-seq differential expression data
4. Annotate cluster directionality
5. Visualize transcriptional shifts across clusters

Required R packages
data.table
dplyr
igraph
GenomicRanges
IRanges
purrr
tibble
FactoMineR
umap
readr
ggplot2

Notes

* Intermediate CSV files generated from earlier steps are used as inputs for downstream analyses.
* UMAP cluster labels may vary slightly between runs and should be manually verified before downstream biological interpretation.
* Scripts were organized as independent analysis modules for clarity and reproducibility.

Citation

If you use this code or analysis framework, please cite the associated manuscript.
