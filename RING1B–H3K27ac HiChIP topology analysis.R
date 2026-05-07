# ============================================================
# RING1B–H3K27ac HiChIP topology analysis
# ============================================================
# Input files should contain the following columns:
# seqnames1
# start1
# seqnames2
# start2

# Example input file names used in the script:
# r1b_hichip.csv
# k27ac_hichip.csv
#
# Required packages:
#   data.table
#   dplyr
#   igraph
#   GenomicRanges
#   IRanges
#   purrr
#   tibble
#   FactoMineR
#   umap
#   ggplot2
#
# Required input columns:
#   seqnames1, start1
#   seqnames2, start2
#
# Analysis overview:
#   1. Build HiChIP interaction graphs
#   2. Identify connected hub clusters
#   3. Compute graph topology features
#   4. Convert anchors to genomic regions
#   5. Identify overlapping hub networks
#   6. Compute structural similarity
#   7. Perform PCA + UMAP + K-means clustering
# ============================================================


# ============================================================
# Load packages and input data
# ============================================================

library(data.table)
library(dplyr)
library(igraph)
library(GenomicRanges)
library(IRanges)
library(purrr)
library(tibble)
library(FactoMineR)
library(umap)
library(ggplot2)

r1b <- fread("r1b_hichip.csv")
k27ac <- fread("k27ac_hichip.csv")


# ============================================================
# Step 1: Build interaction graphs
# ============================================================

make_graph <- function(df) {
  
  edges <- data.frame(
    from = paste(df$seqnames1, df$start1, sep = ":"),
    to   = paste(df$seqnames2, df$start2, sep = ":")
  )
  
  graph_from_data_frame(edges, directed = FALSE)
}

r1b_graph <- make_graph(r1b)
k27ac_graph <- make_graph(k27ac)


# ============================================================
# Step 2: Identify hub clusters
# ============================================================

r1b_clusters <- components(r1b_graph)
k27ac_clusters <- components(k27ac_graph)


# ============================================================
# Step 3: Compute graph topology features
# ============================================================

compute_features <- function(graph, clusters) {
  
  cluster_list <- split(V(graph)$name, clusters$membership)
  
  features <- lapply(names(cluster_list), function(clid) {
    
    nodes <- cluster_list[[clid]]
    
    sg <- induced_subgraph(graph, nodes)
    
    if (vcount(sg) < 3) return(NULL)
    
    tibble(
      cluster_id = clid,
      n_nodes = vcount(sg),
      n_edges = ecount(sg),
      mean_degree = mean(degree(sg)),
      max_degree = max(degree(sg)),
      connectivity_index = mean(degree(sg)) / log1p(vcount(sg)),
      clustering = transitivity(sg, type = "average")
    )
  })
  
  bind_rows(features)
}


r1b_features <- compute_features(r1b_graph, r1b_clusters)
k27ac_features <- compute_features(k27ac_graph, k27ac_clusters)

r1b_features <- dplyr::rename(
  r1b_features,
  r1b_cluster = cluster_id
)

k27ac_features <- dplyr::rename(
  k27ac_features,
  k27ac_cluster = cluster_id
)

# ============================================================
# Step 4: Convert anchors to GRanges
# ============================================================

anchor_to_gr <- function(anchors) {
  
  chr <- sub(":.*", "", anchors)
  start <- as.integer(sub(".*:", "", anchors))
  end <- start + 5000
  
  GRanges(
    seqnames = chr,
    ranges = IRanges(start, end)
  )
}

r1b_cluster_ranges <- split(
  V(r1b_graph)$name,
  r1b_clusters$membership
) %>%
  lapply(anchor_to_gr)

k27ac_cluster_ranges <- split(
  V(k27ac_graph)$name,
  k27ac_clusters$membership
) %>%
  lapply(anchor_to_gr)


# ============================================================
# Step 5: Generate anchor GRanges
# ============================================================

r1b_anchor_gr <- unlist(
  GRangesList(r1b_cluster_ranges),
  use.names = FALSE
)

r1b_anchor_gr$cluster <- rep(
  seq_along(r1b_cluster_ranges),
  lengths(r1b_cluster_ranges)
)

k27ac_anchor_gr <- unlist(
  GRangesList(k27ac_cluster_ranges),
  use.names = FALSE
)

k27ac_anchor_gr$cluster <- rep(
  seq_along(k27ac_cluster_ranges),
  lengths(k27ac_cluster_ranges)
)


# ============================================================
# Step 6: Identify overlapping hub networks
# ============================================================

hits <- findOverlaps(
  r1b_anchor_gr,
  k27ac_anchor_gr,
  ignore.strand = TRUE
)

overlap_df <- tibble(
  
  r1b_cluster = r1b_anchor_gr$cluster[
    queryHits(hits)
  ],
  
  k27ac_cluster = k27ac_anchor_gr$cluster[
    subjectHits(hits)
  ]
  
) %>%
  distinct()

overlap_df <- overlap_df %>%
  mutate(
    r1b_cluster = as.character(r1b_cluster),
    k27ac_cluster = as.character(k27ac_cluster)
  )

r1b_features <- r1b_features %>%
  mutate(r1b_cluster = as.character(r1b_cluster))

k27ac_features <- k27ac_features %>%
  mutate(k27ac_cluster = as.character(k27ac_cluster))

overlap_df <- overlap_df %>%
  filter(
    r1b_cluster %in% r1b_features$r1b_cluster,
    k27ac_cluster %in% k27ac_features$k27ac_cluster
  )


# ============================================================
# Step 7: Merge features and compute structural distance
# ============================================================

overlap_features <- overlap_df %>%
  
  left_join(
    r1b_features,
    by = "r1b_cluster"
  ) %>%
  
  left_join(
    k27ac_features,
    by = "k27ac_cluster",
    suffix = c("_r1b", "_k27ac")
  )

overlap_features <- overlap_features %>%
  
  rowwise() %>%
  
  mutate(
    
    euclidean_dist = sqrt(
      
      sum(
        
        (
          c_across(
            c(
              n_nodes_r1b,
              mean_degree_r1b,
              max_degree_r1b,
              clustering_r1b
            )
          ) -
            
            c_across(
              c(
                n_nodes_k27ac,
                mean_degree_k27ac,
                max_degree_k27ac,
                clustering_k27ac
              )
            )
          
        )^2,
        
        na.rm = TRUE
      )
    )
  ) %>%
  
  ungroup()


# ============================================================
# Step 8: Log transformation and scaling
# ============================================================

log_features <- overlap_features %>%
  
  select(
    starts_with("n_nodes"),
    starts_with("mean_degree"),
    starts_with("max_degree"),
    starts_with("euclidean_dist")
  ) %>%
  
  mutate(
    across(
      everything(),
      log1p
    )
  )

scale_by_group <- function(df, prefix) {
  
  cols <- grep(
    paste0("^", prefix),
    colnames(df),
    value = TRUE
  )
  
  df[cols] <- scale(df[cols])
  
  return(df)
}

scaled_log_features <- log_features %>%
  
  scale_by_group("n_nodes_r1b") %>%
  scale_by_group("n_nodes_k27ac") %>%
  
  scale_by_group("mean_degree_r1b") %>%
  scale_by_group("mean_degree_k27ac") %>%
  
  scale_by_group("max_degree_r1b") %>%
  scale_by_group("max_degree_k27ac") %>%
  
  scale_by_group("euclidean_dist")

scaled_log_features <- scaled_log_features %>%
  
  dplyr::rename(
    n_nodes_r1b_log = n_nodes_r1b,
    n_nodes_k27ac_log = n_nodes_k27ac,
    mean_degree_r1b_log = mean_degree_r1b,
    mean_degree_k27ac_log = mean_degree_k27ac,
    max_degree_r1b_log = max_degree_r1b,
    max_degree_k27ac_log = max_degree_k27ac,
    euclidean_dist_log = euclidean_dist
  )

indexed_scaled_log_features <- bind_cols(
  overlap_features,
  scaled_log_features
)


# ============================================================
# Step 9: Feature selection
# ============================================================

required_cols <- c(
  "n_nodes_r1b_log",
  "n_nodes_k27ac_log",
  "mean_degree_r1b_log",
  "mean_degree_k27ac_log",
  "max_degree_r1b_log",
  "max_degree_k27ac_log",
  "euclidean_dist_log"
)

features <- indexed_scaled_log_features %>%
  select(all_of(required_cols))


# ============================================================
# Step 10: PCA and UMAP
# ============================================================

set.seed(42)

pca_res <- PCA(
  features,
  scale.unit = TRUE,
  ncp = 10,
  graph = FALSE
)

pc_scores <- pca_res$ind$coord[, 1:5]

umap_res <- umap(
  pc_scores,
  n_neighbors = 15,
  min_dist = 0.1
)

umap_df <- as.data.frame(umap_res$layout)

colnames(umap_df) <- c(
  "UMAP1",
  "UMAP2"
)


# ============================================================
# Step 11: K-means clustering
# ============================================================

set.seed(42)

kmeans_res <- kmeans(
  umap_df,
  centers = 6,
  nstart = 25
)

umap_df$cluster <- factor(
  kmeans_res$cluster
)


# ============================================================
# Step 12: Merge metadata
# ============================================================

umap_df <- bind_cols(
  
  umap_df,
  
  features,
  
  indexed_scaled_log_features %>%
    select(
      r1b_cluster,
      k27ac_cluster
    )
)


# ============================================================
# Step 13: Remove outlier clusters
# ============================================================

small_clusters <- names(
  which(
    table(umap_df$cluster) <= 44
  )
)

umap_clean <- umap_df %>%
  
  filter(
    !cluster %in% small_clusters
  ) %>%
  
  mutate(
    cluster = droplevels(cluster)
  )


# ============================================================
# Step 14: Visualization
# ============================================================

ggplot(
  umap_clean,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = cluster
  )
) +
  geom_point(size = 2) +
  theme_classic()

ggplot(
  umap_clean,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = n_nodes_k27ac_log
  )
) +
  geom_point(size = 2) +
  scale_color_gradientn(
    colors = c(
      "#6b5a94",
      "#dbdcdc",
      "#e8a776"
    )
  ) +
  theme_classic()

