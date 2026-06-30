library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


OUTDIR <- "/PATH/"


# ============================================
# 1. Load Seurat Object & Extract Metadata with Individual
# ============================================
obj_path <- "/PATH/BALL_C1C2_BLAST_25PTs_IDXMRD.rds"

BALL.Obj <- readRDS(obj_path)
print(colnames(BALL.Obj@meta.data))
print(table(BALL.Obj$Stage))

# Extract sample-level metadata with Individual for pairing
sample_meta_df <- BALL.Obj@meta.data %>%
  select(Patient.Sample, Stage, Individual) %>%
  distinct() %>%
  mutate(Individual = as.character(Individual)) %>%
  as.data.frame()

cat("Sample metadata:")
print(sample_meta_df)

pair_check <- sample_meta_df %>%
  group_by(Individual) %>%
  summarise(n = n(), stages = paste(sort(Stage), collapse = "/"), .groups = "drop")
cat("Pairing check:")
print(pair_check)

n_cells_per_sample <- table(BALL.Obj$Patient.Sample)




# ============================================
# 2. Pseudo-bulk Peak Counts per Patient.Sample
# ============================================
DefaultAssay(BALL.Obj) <- "ATAC"

pseudo_bulk_peaks <- AggregateExpression(
  BALL.Obj,
  group.by = "Patient.Sample",
  assays = "ATAC",
  slot = "counts",
  verbose = TRUE
)$ATAC

orig_sample_names <- levels(factor(BALL.Obj$Patient.Sample))
sample_names <- colnames(pseudo_bulk_peaks)

rm(BALL.Obj)




# ============================================
# 3. CPM Normalization per Pseudo-bulk Sample
# ============================================
lib_sizes <- colSums(pseudo_bulk_peaks)

valid_samples <- lib_sizes > 0
sample_names_valid <- sample_names

if (any(!valid_samples)) {
  pseudo_bulk_peaks <- pseudo_bulk_peaks[, valid_samples]
  sample_names_valid <- sample_names[valid_samples]
  lib_sizes <- lib_sizes[valid_samples]
}

D <- Diagonal(x = 1e6 / lib_sizes)
cpm <- pseudo_bulk_peaks %*% D
colnames(cpm) <- sample_names_valid
rownames(cpm) <- rownames(pseudo_bulk_peaks)

rm(pseudo_bulk_peaks, D, lib_sizes, valid_samples)




# ============================================
# 4. Log2 Transform
# ============================================
log_cpm <- cpm
log_cpm@x <- log2(log_cpm@x + 1)
dimnames(log_cpm) <- list(rownames(cpm), sample_names_valid)

saveRDS(log_cpm, file.path(OUTDIR, "Normalized_ATACcount_peakLevel_log2CPM.rds"))

rm(cpm)




# ============================================
# 5. Compute Per-Sample Statistics 
# ============================================
sample_meta_df <- sample_meta_df[match(colnames(log_cpm), sample_meta_df$Patient.Sample), ]
sample_means <- colMeans(log_cpm)
log_cpm_dense <- as.matrix(log_cpm)
sample_medians <- apply(log_cpm_dense, 2, median)

sample_stats <- data.frame(
  Patient.Sample = colnames(log_cpm),
  Individual = sample_meta_df$Individual,
  Stage = sample_meta_df$Stage,
  N_Cells = as.vector(n_cells_per_sample[colnames(log_cpm)]),
  Mean_Accessibility = sample_means,
  Median_Accessibility = sample_medians,
  stringsAsFactors = FALSE
)

print(sample_stats)

write.csv(sample_stats, file.path(OUTDIR, "Normalized_ATACcount_PeakLevel_SampleStats_Paired.csv"), row.names = FALSE)




# ============================================
# 6. Paired Wilcoxon Signed-Rank Test
# ============================================
pair_df <- sample_stats %>%
  select(Individual, Stage, Mean_Accessibility, Median_Accessibility) %>%
  tidyr::pivot_wider(
    id_cols = Individual,
    names_from = Stage,
    values_from = c(Mean_Accessibility, Median_Accessibility)
  )

idx_means <- pair_df$Mean_Accessibility_IDX
mrd_means <- pair_df$Mean_Accessibility_MRD
diff_means <- mrd_means - idx_means


wilcox_mean <- wilcox.test(mrd_means, idx_means, paired = TRUE, exact = FALSE)
print(wilcox_mean)

# Save results
test_results <- data.frame(
  Metric = c("Mean"),
  N_pairs = nrow(pair_df),
  IDX_mean = c(mean(idx_means)),
  MRD_mean = c(mean(mrd_means)),
  Diff_mean = c(mean(diff_means)),
  V_statistic = c(wilcox_mean$statistic),
  P_value = c(wilcox_mean$p.value),
  stringsAsFactors = FALSE
)

write.csv(test_results, file.path(OUTDIR, "Normalized_ATACcount_PeakLevel_PairedWilcoxonResults.csv"), row.names = FALSE)




# ============================================
# 7. Plot
# ============================================
sample_stats_long <- sample_stats %>%
  tidyr::pivot_longer(
    cols = c(Mean_Accessibility),
    names_to = "Metric",
    values_to = "Accessibility"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Mean_Accessibility"),
                    labels = c("Mean")),
    Stage = factor(Stage, levels = c("IDX", "MRD"))
  )


p_mean_data <- sample_stats_long %>% filter(Metric == "Mean")
p_mean <- ggplot(p_mean_data, aes(x = Stage, y = Accessibility)) +
  geom_violin(aes(fill = Stage), trim = FALSE, color = NA, alpha = 0.9, show.legend = FALSE) +

  geom_boxplot(width = 0.15, fill = "white", color = "#222222", 
               outlier.shape = NA, linewidth = 0.6, show.legend = FALSE) +
  geom_line(aes(group = Individual), color = "gray30", alpha = 0.5, linewidth = 0.6) +
  geom_point(aes(fill = Stage), size = 3, shape = 21, color = "black", stroke = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = c("IDX" = "#B75151", "MRD" = "#5E758C")) +
  labs(
    title = "Peak-Level Accessibility: Paired Per-Sample Mean",
    subtitle = paste("Paired Wilcoxon p =", format(wilcox_mean$p.value, digits = 3, scientific = TRUE)),
    y = "Mean log2(CPM + 1)", x = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 13),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black")
  )

ggsave(file.path(OUTDIR, "Normalized_ATACcount_PeakLevel_PairedViolin_Mean.pdf"),p_mean, width = 4, height = 5.5, units = "in", dpi = 300)