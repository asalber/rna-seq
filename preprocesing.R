# Data preprocesing 
# Alfredo Sánchez Alberca (asalber@ceu.es)

# Libraries
library(tidyverse)
library(edgeR)

# Data loading
df_reads <- read_csv("data/counts.csv")

# Rename individuals
df_groups <- read_csv("data/groups.csv") %>%
    mutate(Group = as.factor(Group))
new_names <- df_groups %>%
    select("Id", "Name") %>%
    deframe()
df_reads <- rename_with(df_reads, ~ df_groups$Id, all_of(df_groups$Name))

# Create count matrix
matrix_reads <- as.data.frame(select(df_reads, -CLUSTER))
rownames(matrix_reads) <- df_reads$CLUSTER

# Add experimental group
df_reads_longer <- df_reads %>%
    pivot_longer(cols = C1:DG2T7, names_to = "Id", values_to = "Count") %>%
    left_join(df_groups, by = "Id")

# reads loading
reads <- DGEList(counts = matrix_reads, group = df_groups$Group)

## FILTERING
# Filtering genes with more than 1 reads per million in at least 7 individuals
keep <- rowSums(cpm(reads)> 1) >= 7

# Automatic filtering
#keep <- filterByExpr(reads)
reads_filtered <- reads[keep, , keep.lib.sizes=F]

# Filter also the longer data frame
df_reads_filtered_longer <- df_reads_longer %>%
    filter(CLUSTER %in% names(keep[keep==T]))

# Reads counts per million of reads
reads_filtered_cpm <- cpm(reads)
reads_filtered_cpm_longer <- as.data.frame(reads_filtered_cpm) %>%
    pivot_longer(cols = C1:DG2T7, names_to = "Id", values_to = "CPM") %>%
    left_join(df_groups, by = "Id")
