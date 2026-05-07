# ============================================================
# Integration of HiChIP clusters with gene expression changes
# ============================================================
# This script integrates overlapping RING1B–H3K27ac
# interaction clusters with transcriptional changes from
# RING1A/B knockout RNA-seq data.
#
# Input files:
#   tss.csv
#   R1A_DEGs.csv
#   overlap_df_filtered_r1b.csv
#   umap_no_outliers.csv
#   cluster_direction.csv
#
# Required packages:
#   dplyr
#   readr
#   tibble
#   GenomicRanges
#   IRanges
#   ggplot2
#
# Required input columns:
#   tss.csv:
#       seqnames, tss, genes
#
#   R1A_DEGs.csv:
#       genes, log2FoldChange
#
# Analysis overview:
#   1. Map RING1B anchors to nearby TSS regions
#   2. Assign genes to overlapping HiChIP clusters
#   3. Integrate RNA-seq differential expression data
#   4. Annotate cluster directionality
#   5. Visualize transcriptional shifts across clusters
# ============================================================


# ============================================================
# Load packages and input data
# ============================================================

library(dplyr)
library(readr)
library(tibble)
library(GenomicRanges)
library(IRanges)
library(ggplot2)

tss <- read_csv(
  "tss.csv"
)

deg <- read_csv(
  "R1A_DEGs.csv"
)

overlap_df_filtered_r1b <- read_csv(
  "overlap_df_filtered_r1b.csv"
)

umap_no_outliers <- read_csv(
  "umap_no_outliers.csv"
)

cluster_direction <- read_csv(
  "cluster_direction.csv"
)


# ============================================================
# Step 1: Convert anchors and TSS regions to GRanges
# ============================================================

anchor_gr <- GRanges(
  
  seqnames = overlap_df_filtered_r1b$seqnames,
  
  ranges = IRanges(
    start = overlap_df_filtered_r1b$start,
    end = overlap_df_filtered_r1b$end
  ),
  
  r1b_cluster = overlap_df_filtered_r1b$r1b_cluster,
  
  k27ac_cluster = overlap_df_filtered_r1b$k27ac_cluster
)

tss_gr <- GRanges(
  
  seqnames = tss$seqnames,
  
  ranges = IRanges(
    start = tss$tss,
    end = tss$tss
  ),
  
  gene = tss$genes
)


# ============================================================
# Step 2: Identify anchor–TSS overlaps
# ============================================================

hits2 <- findOverlaps(
  anchor_gr,
  tss_gr
)

gene_anchor_df <- tibble(
  
  gene = mcols(tss_gr)$gene[
    subjectHits(hits2)
  ],
  
  r1b_cluster = mcols(anchor_gr)$r1b_cluster[
    queryHits(hits2)
  ],
  
  k27ac_cluster = mcols(anchor_gr)$k27ac_cluster[
    queryHits(hits2)
  ]
  
) %>%
  
  distinct()


# ============================================================
# Step 3: Integrate UMAP cluster annotations
# ============================================================

gene_group_df <- gene_anchor_df %>%
  
  left_join(
    
    umap_no_outliers %>%
      
      select(
        cluster,
        r1b_cluster,
        k27ac_cluster
      ),
    
    by = c(
      "r1b_cluster",
      "k27ac_cluster"
    )
  )


# ============================================================
# Step 4: Select target UMAP clusters
# ============================================================

# NOTE:
# UMAP cluster labels may vary between runs.
# Cluster identities should be manually verified
# before downstream analysis.

gene_group_df <- gene_group_df %>%
  
  filter(
    cluster %in% c(2, 3, 5, 6)
  )


# ============================================================
# Step 5: Integrate RNA-seq differential expression
# ============================================================

gene_group_df <- gene_group_df %>%
  
  left_join(
    
    deg,
    
    by = c(
      "gene" = "genes"
    )
  ) %>%
  
  filter(
    !is.na(log2FoldChange)
  )


# ============================================================
# Step 6: Remove duplicated entries
# ============================================================

gene_group_df <- gene_group_df %>%
  
  distinct(
    gene,
    r1b_cluster,
    k27ac_cluster,
    cluster,
    .keep_all = TRUE
  )


# ============================================================
# Step 7: Annotate cluster directionality
# ============================================================

gene_group_df <- gene_group_df %>%
  
  left_join(
    
    cluster_direction %>%
      
      select(
        r1b_cluster,
        k27ac_cluster,
        cluster_direction
      ),
    
    by = c(
      "r1b_cluster",
      "k27ac_cluster"
    )
  )


# ============================================================
# Step 8: Export integrated dataset
# ============================================================

write.csv(
  gene_group_df,
  "R1Banchor_R1Arna_gene_group_df.csv",
  row.names = FALSE
)


# ============================================================
# Step 9: Sanity check
# ============================================================

table(
  gene_group_df$cluster_direction,
  useNA = "ifany"
)


# ============================================================
# Step 10: Prepare visualization data
# ============================================================

plot_df <- gene_group_df %>%
  
  mutate(
    sig = abs(log2FoldChange) > 1
  )


# ============================================================
# Step 11: ECDF visualization
# ============================================================

ggplot(
  
  plot_df,
  
  aes(
    x = log2FoldChange
  )
  
) +
  
  stat_ecdf(
    
    aes(
      color = sig
    ),
    
    geom = "step",
    linewidth = 1
  ) +
  
  scale_color_manual(
    
    values = c(
      `TRUE` = "#e8a776",
      `FALSE` = "#8174a9"
    ),
    
    guide = "none"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  facet_wrap(
    
    ~ cluster,
    
    labeller = labeller(
      
      cluster = function(x) {
        
        paste0(
          
          "Cluster ",
          
          x,
          
          " (n = ",
          
          sapply(
            x,
            
            function(cl) {
              
              plot_df %>%
                
                filter(
                  cluster == cl,
                  sig
                ) %>%
                
                distinct(gene) %>%
                
                nrow()
            }
          ),
          
          ")"
        )
      }
    )
  ) +
  
  theme_classic() +
  
  labs(
    x = "log2FoldChange",
    y = "Cumulative fraction of genes"
  )
