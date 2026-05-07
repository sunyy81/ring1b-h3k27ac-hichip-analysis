# ============================================================
# RING1B–H3K27ac interaction directionality analysis
# ============================================================
# This script evaluates directional consistency between
# overlapping RING1B and H3K27ac interaction clusters.
#
# Input files:
#   overlap_df_filtered_r1b.csv
#   overlap_df_filtered_k27ac.csv
#
# These files were generated from:
#   RING1B–H3K27ac HiChIP topology analysis
#
# Required packages:
#   dplyr
#   ggplot2
#
# Required input columns:
#   r1b_cluster
#   k27ac_cluster
#   seqnames
#   start
#   end
#
# Analysis overview:
#   1. Identify shared anchors between cluster pairs
#   2. Compute anchor-relative genomic distances
#   3. Calculate anchor-level directional consistency
#   4. Assign cluster-level interaction directionality
#   5. Visualize cluster direction categories
# ============================================================


# ============================================================
# Load packages and input data
# ============================================================

library(dplyr)
library(ggplot2)

overlap_df_filtered_r1b <- read.csv(
  "overlap_df_filtered_r1b.csv"
)

overlap_df_filtered_k27ac <- read.csv(
  "overlap_df_filtered_k27ac.csv"
)


# ============================================================
# Step 1: Identify shared anchors
# ============================================================

shared_anchors <- inner_join(
  
  overlap_df_filtered_r1b,
  
  overlap_df_filtered_k27ac,
  
  by = c(
    "r1b_cluster",
    "k27ac_cluster",
    "seqnames",
    "start",
    "end"
  )
  
) %>%
  
  mutate(
    ref_pos = start,
    anchor_type = "exact"
  ) %>%
  
  select(
    r1b_cluster,
    k27ac_cluster,
    seqnames,
    ref_pos,
    anchor_type
  )


# ============================================================
# Step 2: Compute RING1B anchor distances
# ============================================================

shared_overlap_df_filtered_r1b <- overlap_df_filtered_r1b %>%
  
  semi_join(
    shared_anchors,
    by = c(
      "r1b_cluster",
      "k27ac_cluster",
      "seqnames"
    )
  ) %>%
  
  group_by(
    r1b_cluster,
    k27ac_cluster,
    seqnames
  ) %>%
  
  summarise(
    r1b_anchor_starts = list(start),
    .groups = "drop"
  )

r1b_direction <- shared_anchors %>%
  
  left_join(
    
    shared_overlap_df_filtered_r1b,
    
    by = c(
      "r1b_cluster",
      "k27ac_cluster",
      "seqnames"
    )
  ) %>%
  
  rowwise() %>%
  
  mutate(
    r1b_dists = list(
      r1b_anchor_starts - ref_pos
    )
  ) %>%
  
  ungroup()


# ============================================================
# Step 3: Compute H3K27ac anchor distances
# ============================================================

shared_overlap_df_filtered_k27ac <- overlap_df_filtered_k27ac %>%
  
  semi_join(
    shared_anchors,
    by = c(
      "r1b_cluster",
      "k27ac_cluster",
      "seqnames"
    )
  ) %>%
  
  group_by(
    r1b_cluster,
    k27ac_cluster,
    seqnames
  ) %>%
  
  summarise(
    k27ac_anchor_starts = list(start),
    .groups = "drop"
  )

k27ac_direction <- shared_anchors %>%
  
  left_join(
    
    shared_overlap_df_filtered_k27ac,
    
    by = c(
      "r1b_cluster",
      "k27ac_cluster",
      "seqnames"
    )
  ) %>%
  
  rowwise() %>%
  
  mutate(
    k27ac_dists = list(
      k27ac_anchor_starts - ref_pos
    )
  ) %>%
  
  ungroup()


# ============================================================
# Step 4: Merge anchor distance profiles
# ============================================================

shared_anchor_combined <- r1b_direction %>%
  
  select(
    r1b_cluster,
    k27ac_cluster,
    seqnames,
    ref_pos,
    r1b_dists
  ) %>%
  
  left_join(
    
    k27ac_direction %>%
      
      select(
        r1b_cluster,
        k27ac_cluster,
        seqnames,
        ref_pos,
        k27ac_dists
      ),
    
    by = c(
      "r1b_cluster",
      "k27ac_cluster",
      "seqnames",
      "ref_pos"
    )
  )


# ============================================================
# Step 5: Compute anchor-level directionality
# ============================================================

direction_per_anchor <- shared_anchor_combined %>%
  
  rowwise() %>%
  
  mutate(
    
    r1b_pos_frac = ifelse(
      length(r1b_dists) > 0,
      sum(r1b_dists > 0) / length(r1b_dists),
      NA_real_
    ),
    
    k27ac_pos_frac = ifelse(
      length(k27ac_dists) > 0,
      sum(k27ac_dists > 0) / length(k27ac_dists),
      NA_real_
    ),
    
    anchor_direction = case_when(
      
      r1b_pos_frac >= 2/3 &
        k27ac_pos_frac >= 2/3 ~ "Same",
      
      r1b_pos_frac <= 1/3 &
        k27ac_pos_frac <= 1/3 ~ "Same",
      
      r1b_pos_frac >= 2/3 &
        k27ac_pos_frac <= 1/3 ~ "Opposite",
      
      r1b_pos_frac <= 1/3 &
        k27ac_pos_frac >= 2/3 ~ "Opposite",
      
      TRUE ~ "Ambiguous"
    )
  ) %>%
  
  ungroup()


# ============================================================
# Step 6: Assign cluster-level directionality
# ============================================================

cluster_direction <- direction_per_anchor %>%
  
  group_by(
    r1b_cluster,
    k27ac_cluster
  ) %>%
  
  summarise(
    
    n_same = sum(
      anchor_direction == "Same"
    ),
    
    n_opposite = sum(
      anchor_direction == "Opposite"
    ),
    
    n_total = n(),
    
    cluster_direction = case_when(
      
      n_same / n_total >= 2/3 ~ "Same",
      
      n_opposite / n_total >= 2/3 ~ "Opposite",
      
      TRUE ~ "Ambiguous"
    ),
    
    .groups = "drop"
  )


# ============================================================
# Step 7: Export results
# ============================================================

write.csv(
  cluster_direction,
  "cluster_direction.csv",
  row.names = FALSE
)


# ============================================================
# Step 8: Summary statistics
# ============================================================

cluster_direction %>%
  
  count(cluster_direction) %>%
  
  mutate(
    prop = n / sum(n)
  )


# ============================================================
# Step 9: Visualization
# ============================================================

cluster_direction %>%
  
  count(cluster_direction) %>%
  
  ggplot(
    aes(
      x = cluster_direction,
      y = n,
      fill = cluster_direction
    )
  ) +
  
  geom_col(width = 0.6) +
  
  labs(
    x = NULL,
    y = "Number of cluster pairs"
  ) +
  
  theme_classic()
