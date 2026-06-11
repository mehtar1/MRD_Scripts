# =============================================================================
# Mutation Evidence in Myeloid-like and Non-B-lineage Blast Populations
# =============================================================================
#
# Description:
#   Visualises mutation and CNV evidence across leukemic blast compartments
#   to validate the blast identity of non-B-lineage (myeloid-like) cells.
#   Mutations are sourced from Mission Bio targeted scDNA-seq and then
#   confirmed in the multiome dataset using vartrix, so each cell can be
#   annotated with both its arrest-state projection label (from scRNA-seq)
#   and its somatic mutation status (from long-read vartrix output).
#
#   The script produces four outputs:
#     1. Tile plot  — mutations detected in myeloid-like blast cells, faceted
#                    by subtype / patient / gene, tile color encodes cell-type
#                    coverage, point shape encodes mutation evidence level
#     2. CNV plot (chr21, iAMP21 only) — median Numbat CNV probability for
#                    the chr21(amp) segment across cell types per sample
#     3. CNV plot (all events) — events shared across HSPC-like, lymphoid-like,
#                    and myeloid-like blast populations (high confidence)
#     4. Single-cell annotation panel — healthy reference UMAP, per-patient
#                    mutation status UMAPs, and chr21-amp probability UMAPs
#                    assembled with patchwork for a multi-panel figure
#
# Inputs (set in USER CONFIGURATION):
#   - gex_df.tsv                 : per-cell mutation evidence table
#   - metadata_df.tsv            : per-sample metadata (subtype, USI, stage)
#   - cnv_prob_df.tsv            : Numbat CNV probabilities per cell type
#   - sc_annotations_df.tsv      : per-cell UMAP coordinates and annotations
#   - healthy_cell_types_df.tsv  : healthy reference UMAP coordinates
#   - filtered_mutations_gt0_meta_df.tsv : full mutation annotation table
#   - Numbat joint_post_1.tsv files (one per sample, in numbat_dir)
#   - Per-sample mutation-by-barcode TSV files (in cell_barcode_mutation_dir)
#
# Outputs:
#   - main_fig.myeloid_events.pdf       : tile plot + chr21 CNV plot
#   - main_fig.single_cell_annotations.pdf : healthy ref + mutation + CNV UMAPs
#
# Dependencies:
#   tidyverse, Seurat, patchwork, ggh4x, scales, ggrepel
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(tidyverse)
library(Seurat)
library(patchwork)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Project root — all data paths are derived from here via here::here()
# Set this to the root of the VHR-ALL analysis directory
vhr_dir <- here::here()

# Sub-directory structure (relative to vhr_dir)
analysis_dir          <- here::here(vhr_dir, "data_analysis")
manuscript_dir        <- here::here(vhr_dir, "manuscript")
manuscript_data_dir   <- here::here(manuscript_dir, "data_used_in_plots")
main_fig1_data_dir    <- here::here(manuscript_data_dir, "main_fig1")
cell_barcode_mutation_dir <- here::here(main_fig1_data_dir, "mutations_by_cell_type")
numbat_dir            <- here::here(analysis_dir, "numbat_analysis", "outputs",
                                     "Step2_run_numbat")
plot_dir              <- here::here(manuscript_dir, "plot_pdfs")

# Sample-level focus for the single-cell annotation panel (Section 4).
# Set to the patient sample IDs and specific variant IDs to highlight.
# Example below shows PAVEPS IDX/MRD with the KRAS mutation.
SC_ANNO_SAMPLE_IDX   <- "PAVEPS_IDX"
SC_ANNO_SAMPLE_MRD   <- "PAVEPS_MRD"
SC_ANNO_VARIANT      <- "chr12:25245349:A:ACCT"   # KRAS variant GRCh38 ID

# Output figure dimensions
FIG_MYELOID_HEIGHT <- 4     # inches
FIG_MYELOID_WIDTH  <- 7.5   # inches
FIG_ANNO_HEIGHT    <- 8     # inches
FIG_ANNO_WIDTH     <- 15    # inches


# ── Shared aesthetics ─────────────────────────────────────────────────────────
# These encode cell type × coverage level combinations used in the tile plot.
# Colour encodes lineage; shade encodes number of cells with coverage at site.

shape_vector <- c("0" = 4, "1" = 1, "2+" = 19)

fill_vector <- c(
  "HSPC 0 cells"       = "#fee5d9", "HSPC <10 cells"      = "#fcae91",
  "HSPC 10+ cells"     = "#fb6a4a",
  "Lymphoid 0 cells"   = "#edf8e9", "Lymphoid <10 cells"  = "#bae4b3",
  "Lymphoid 10+ cells" = "#74c476",
  "Myeloid 0 cells"    = "#eff3ff", "Myeloid <10 cells"   = "#bdd7e7",
  "Myeloid 10+ cells"  = "#6baed6",
  "Normal 0 cells"     = "#f2f0f7", "Normal <10 cells"    = "#cbc9e2",
  "Normal 10+ cells"   = "#9e9ac8"
)

cell_type_with_coverage_levels <- names(fill_vector)

# Shared ggplot2 theme for tile plots
tile_theme <- theme_bw() %+replace% theme(
  strip.text.y.left  = element_text(angle = 0),
  strip.background   = element_blank(),
  strip.placement    = "outside",
  ggh4x.facet.nestline = element_line(),
  legend.position    = "right",
  legend.direction   = "vertical",
  legend.text        = element_text(size = 6),
  legend.title       = element_text(size = 8),
  panel.grid         = element_blank(),
  axis.text.y        = element_text(size = 6),
  axis.text.x        = element_text(angle = 45, hjust = 1)
)

# Shared ggplot2 theme for UMAP plots
umap_theme <- theme_bw() %+replace% theme(
  panel.grid  = element_blank(),
  axis.text   = element_blank(),
  axis.ticks  = element_blank()
)


# =============================================================================
# SECTION 1: Load data
# =============================================================================

cat("Loading data...\n")

gex_df                       <- readr::read_tsv(here::here(main_fig1_data_dir, "gex_df.tsv"))
metadata_df                  <- readr::read_tsv(here::here(main_fig1_data_dir, "metadata_df.tsv"))
cnv_prob_df                  <- readr::read_tsv(here::here(main_fig1_data_dir, "cnv_prob_df.tsv"))
sc_annotations_df            <- readr::read_tsv(here::here(main_fig1_data_dir, "sc_annotations_df.tsv"))
healthy_cell_types_df        <- readr::read_tsv(here::here(main_fig1_data_dir, "healthy_cell_types_df.tsv"))
filtered_mutations_gt0_meta_df <- readr::read_tsv(here::here(main_fig1_data_dir,
                                                               "filtered_mutations_gt0_meta_df.tsv"))

# Build sample IDs and Numbat file paths from metadata
samples <- metadata_df |>
  mutate(file_prefix = str_c(USI, stage, sep = "_")) |>
  pull(file_prefix)

joint_post_1_files_list <- samples |>
  map(\(x) here::here(numbat_dir, x, "joint_post_1.tsv")) |>
  setNames(samples)

# Split sc annotations by patient sample for per-sample plotting
sc_annotations_list <- sc_annotations_df |>
  split(sc_annotations_df$Patient.Sample)


# =============================================================================
# SECTION 2: Myeloid blast mutation tile plot
# =============================================================================
# Identifies mutations with evidence in myeloid-like blast cells, then
# plots a tile grid (cell type × stage) faceted by subtype / patient / gene.
# Tile fill encodes how many cells have coverage at the variant site;
# point shape encodes the number of cells with the mutation observed.
#
# myeloid_mutations_plot2 is the manuscript figure version. It uses cell-type
# fill colour with alpha-transparency for coverage level, which is cleaner
# in print than the explicit coverage colour scale in the exploratory version.

mutations_with_myeloid <- gex_df |>
  filter(updated_level_one == "Myeloid-like Blasts",
         mutation_evidence != "0") |>
  select(mutation, sample) |>
  unique()

myeloid_mutations_plot <- mutations_with_myeloid |>
  left_join(gex_df, by = c("mutation", "sample")) |>
  mutate(n_cells_color = factor(n_cells_color,
                                levels  = c("0 cells", "<10 cells", "10+ cells"),
                                ordered = TRUE)) |>
  ggplot(aes(x = updated_level_one, y = fct_rev(updated_stage))) +
  geom_tile(aes(fill = updated_level_one, alpha = n_cells_color),
            show.legend = TRUE) +
  geom_point(aes(shape = mutation_evidence)) +
  ggh4x::facet_nested(
    subtype + preferred_patient_name + VEP_SYMBOL.x ~ .,
    space     = "free",
    scales    = "free",
    drop      = TRUE,
    switch    = "y",
    solo_line = TRUE,
    strip     = ggh4x::strip_nested(size = "variable")
  ) +
  scale_shape_manual(values = shape_vector) +
  labs(x      = "Cell type",
       y      = NULL,
       fill   = "Cell type",
       alpha  = "Cells with\ncoverage at site",
       title  = "Myeloid-like blast mutations",
       shape  = "Cells with\nmutation observed") +
  guides(shape = guide_legend(override.aes = list(fill = NA))) +
  tile_theme %+replace%
  theme(plot.title.position = "plot",
        panel.border        = element_blank(),
        panel.background    = element_blank(),
        axis.text.x         = element_blank(),
        axis.ticks          = element_blank())


# =============================================================================
# SECTION 3: CNV probability plots (Numbat)
# =============================================================================
# Shared plot construction for CNV probability across cell types.
# Each point = one cell type in one sample; size encodes number of cells.

# Helper: build a CNV probability dot plot
make_cnv_plot <- function(plot_data, title) {
  ggplot(plot_data, aes(x = med, y = fct_rev(stage))) +
    geom_point(aes(size = n, color = level_one), shape = 21) +
    geom_point(aes(size = n, color = level_one), shape = ".", show.legend = TRUE) +
    ggh4x::facet_nested(
      subtype + preferred_patient_name ~ .,
      space     = "free",
      scales    = "free",
      drop      = TRUE,
      switch    = "y",
      solo_line = TRUE,
      strip     = ggh4x::strip_nested(size = "variable")
    ) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent,
                       expand = c(0, 0)) +
    labs(x      = "Median probability",
         y      = NULL,
         color  = "Cell type",
         size   = "Number of cells",
         title  = title) +
    coord_cartesian(clip = "off") +
    theme_bw() %+replace%
    theme(plot.title.position   = "plot",
          strip.text.y.left     = element_text(angle = 0),
          strip.background      = element_blank(),
          strip.placement       = "outside",
          strip.clip            = "off",
          ggh4x.facet.nestline  = element_line(),
          legend.direction      = "vertical",
          panel.grid.major.x    = element_blank(),
          panel.grid.minor.x    = element_blank(),
          panel.border          = element_blank(),
          panel.background      = element_blank(),
          axis.ticks.y          = element_blank())
}

# -- chr21 amplification (iAMP21 subtype only) --------------------------------
# For each iAMP21 sample, selects the chr21 amp segment with the highest
# median CNV probability in lymphoid-like blasts as the representative event.

CNV_chr21_amp_data <- metadata_df |>
  filter(subtype %in% c("iAMP21", "iAMP21 + Fusion")) |>
  mutate(multiome_sample_name = ifelse(str_starts(USI, "P"),
                                       sample_name,
                                       str_c(USI, stage, sep = "_"))) |>
  left_join(
    cnv_prob_df |>
      filter(level_one == "Lymphoid-like Blast", chr == "21", type == "(amp)") |>
      group_by(sample) |>
      summarize(seg_label = seg_label[which.max(med)]) |>
      left_join(cnv_prob_df, by = c("sample", "seg_label")),
    by = c("multiome_sample_name" = "sample")
  ) |>
  arrange(desc(n)) |>
  filter(!is.na(med)) |>
  mutate(level_one = ifelse(level_one == "Not-Blast", "Normal", level_one))

CNV_chr21_amp_plot <- make_cnv_plot(
  CNV_chr21_amp_data,
  "Probability of chr21(amp) across cell types"
)

# -- All high-confidence CNV events (shared across blast compartments) --------
# Keeps only CNV events where HSPC-like, lymphoid-like, AND myeloid-like
# blasts all show >= 90% median probability (but normal cells < 90%).
# This selects events that are clonal across blast lineages.

CNV_all_data <- metadata_df |>
  mutate(multiome_sample_name = ifelse(str_starts(USI, "P"),
                                       sample_name,
                                       str_c(USI, stage, sep = "_"))) |>
  left_join(
    cnv_prob_df |>
      pivot_wider(id_cols     = c(sample, seg_label, chr, type, relabel),
                  names_from  = level_one,
                  values_from = c(med, n)) |>
      filter(`med_HSPC-like Blast`      >= 0.9,
             `med_Lymphoid-like Blast`  >= 0.9,
             `med_Myeloid-like Blast`   >= 0.9,
             `med_Not-Blast`            <  0.9) |>
      pivot_longer(cols = starts_with(c("med_", "n_"))) |>
      separate(name, into = c("measure", "level_one"), sep = "_") |>
      filter(!is.na(value)) |>
      pivot_wider(id_cols    = c(sample, seg_label, chr, type, relabel, level_one),
                  names_from = measure,
                  values_from = value),
    by = c("multiome_sample_name" = "sample")
  ) |>
  arrange(desc(n)) |>
  filter(!is.na(med))

CNV_all_plot <- make_cnv_plot(
  CNV_all_data,
  "Probability of CNV events across cell types"
) %+replace%
  theme(strip.text.y.left  = element_text(hjust = 1, angle = 0, size = 8),
        legend.text         = element_text(size = 6),
        legend.title        = element_text(size = 8),
        axis.text.y         = element_text(size = 6))


# =============================================================================
# SECTION 4: Per-sample single-cell annotation functions
# =============================================================================
# These functions produce UMAP plots projecting patient cells onto the healthy
# donor reference space. They are used in the multi-panel figure assembly.

plot_healthy_reference <- function(healthy_df, show_legend = FALSE, seed = 1) {
  set.seed(seed)
  text_labels <- healthy_df |>
    group_by(cell_type) |>
    summarize(median_UMAP1 = median(UMAP_1),
              median_UMAP2 = median(UMAP_2),
              n_cells      = n(),
              .groups      = "drop")

  healthy_df |>
    ggplot(aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = cell_type), shape = 16,
               show.legend = ifelse(show_legend, TRUE, FALSE)) +
    ggrepel::geom_label_repel(
      data = text_labels,
      aes(x = median_UMAP1, y = median_UMAP2,
          label = cell_type, color = cell_type, size = n_cells),
      show.legend = FALSE
    ) +
    labs(title = "Healthy Reference Cell Types",
         x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    umap_theme
}

plot_mutation_status <- function(annot_df, variant, sample,
                                  healthy_df, display_mutation) {
  # Reads the per-barcode mutation TSV produced by vartrix, identifies cells
  # supporting the ALT allele, REF allele, or lacking coverage, then overlays
  # them on the healthy reference UMAP.
  mut_df <- readr::read_tsv(
    here::here(cell_barcode_mutation_dir,
               str_c(sample, ".gex.mutations_by_barcode.tsv")),
    show_col_types = FALSE
  ) |>
    pivot_longer(cols = starts_with("chr")) |>
    select(cell_barcode = barcode,
           gene.mutation = name,
           supports_ALT  = value) |>
    filter(!is.na(supports_ALT), gene.mutation == variant) |>
    mutate(supports_ALT = as.logical(supports_ALT))

  cells_support_ALT <- str_c(sample, "_",
                               mut_df |> filter(supports_ALT) |> pull(cell_barcode))
  cells_support_REF <- str_c(sample, "_",
                               mut_df |> filter(!supports_ALT) |> pull(cell_barcode))

  annot_df |>
    mutate(mutation_status = case_when(
      patient_bc %in% cells_support_ALT ~ "ALT",
      patient_bc %in% cells_support_REF ~ "REF",
      TRUE                              ~ "No coverage"
    )) |>
    mutate(mutation_status = factor(mutation_status,
                                    levels  = c("No coverage", "REF", "ALT"),
                                    ordered = TRUE)) |>
    arrange(mutation_status) |>
    ggplot(aes(x = refUMAP_1, y = refUMAP_2)) +
    geom_point(data = healthy_df,
               aes(x = UMAP_1, y = UMAP_2), shape = 16, color = "grey90") +
    geom_point(aes(color = mutation_status), shape = 16) +
    scale_color_manual(values = c("grey50", "black", "red"), drop = FALSE) +
    labs(x     = "UMAP 1",
         y     = "UMAP 2",
         color = "Allele",
         title = str_c(str_replace(sample, "_", " "), ": ", display_mutation)) +
    umap_theme
}

plot_specific_cnv_event <- function(nb, annot_df, sample, cnv_event,
                                     healthy_df, display_event) {
  # Overlays per-cell CNV probability for a specific Numbat segment
  # onto the healthy reference UMAP.
  plot_df <- nb |>
    filter(seg_label == cnv_event) |>
    mutate(cell = str_c(sample, "_", cell)) |>
    select(cell, p_cnv, seg_label) |>
    left_join(annot_df, by = c("cell" = "patient_bc")) |>
    arrange(desc(p_cnv))

  ggplot(plot_df, aes(x = refUMAP_1, y = refUMAP_2, color = p_cnv)) +
    geom_point(data = healthy_df,
               aes(x = UMAP_1, y = UMAP_2), shape = 16, color = "grey90") +
    geom_point(shape = 16) +
    scale_color_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x     = "UMAP 1",
         y     = "UMAP 2",
         color = "CNV probability",
         title = str_c(str_replace(sample, "_", " "), ": ", display_event)) +
    umap_theme
}


# =============================================================================
# SECTION 5: Build per-sample mutation and CNV UMAP plots
# =============================================================================
# Generates mutation-status UMAPs for every sample × mutation combination,
# and chr21-amp CNV UMAPs for iAMP21 samples, storing them in named lists
# for flexible figure assembly.

healthy_reference_plot <- plot_healthy_reference(healthy_cell_types_df)

# -- Per-sample mutation status UMAPs -----------------------------------------
compartment_mutations_df <- gex_df |>
  select(sample, mutation) |>
  unique()

compartment_mutations_list <- compartment_mutations_df |>
  split(f = compartment_mutations_df$sample)

compartment_mutation_plot_list <- names(compartment_mutations_list) |>
  map(\(samp) {
    compartment_mutations_list[[samp]]$mutation |>
      map(\(mutation) {
        display_mutation <- filtered_mutations_gt0_meta_df |>
          filter(sample == samp, scDNA_variant_GRCh38 == mutation) |>
          select(VEP_SYMBOL, VEP_HGVSp, VEP_Protein_position,
                 VEP_Amino_acids, VEP_Consequence) |>
          replace_na(list(VEP_Consequence = "ALT missing",
                          VEP_Amino_acids = "/")) |>
          mutate(VEP_Amino_acids = case_when(
            !str_detect(VEP_Amino_acids, "/") ~
              str_c(VEP_Amino_acids, "/", VEP_Amino_acids),
            .default = VEP_Amino_acids
          )) |>
          separate(VEP_Amino_acids, into = c("before", "after"), sep = "/") |>
          mutate(display_mutation = case_when(
            !is.na(VEP_HGVSp) ~
              str_c(VEP_SYMBOL, str_c("p.", before, VEP_Protein_position, after),
                    sep = " "),
            TRUE ~
              str_c(VEP_SYMBOL, " (", str_replace(VEP_Consequence, "_", " "), ")")
          )) |>
          pull(display_mutation)

        plot_mutation_status(
          annot_df         = sc_annotations_list[[samp]],
          variant          = mutation,
          sample           = samp,
          healthy_df       = healthy_cell_types_df,
          display_mutation = display_mutation
        )
      }) |>
      setNames(compartment_mutations_list[[samp]]$mutation)
  }) |>
  setNames(names(compartment_mutations_list))

# -- Chr21 amp CNV UMAPs (iAMP21 samples) -------------------------------------
chr21_amp_df <- metadata_df |>
  filter(subtype %in% c("iAMP21", "iAMP21 + Fusion")) |>
  mutate(multiome_sample_name = ifelse(str_starts(USI, "P"),
                                       sample_name,
                                       str_c(USI, stage, sep = "_"))) |>
  left_join(
    cnv_prob_df |>
      filter(level_one == "Lymphoid-like Blast", chr == "21", type == "(amp)") |>
      group_by(sample) |>
      summarize(seg_label = seg_label[which.max(med)]) |>
      left_join(cnv_prob_df, by = c("sample", "seg_label")),
    by = c("multiome_sample_name" = "sample")
  ) |>
  arrange(desc(n)) |>
  filter(!is.na(med)) |>
  select(multiome_sample_name, seg_label, relabel) |>
  unique()

chr21_amp_list <- chr21_amp_df |>
  split(f = chr21_amp_df$multiome_sample_name)

chr21_amp_plot_list <- names(chr21_amp_list) |>
  map(\(x) {
    samp      <- chr21_amp_list[[x]]$multiome_sample_name
    seg_label <- chr21_amp_list[[x]]$seg_label
    relabel   <- chr21_amp_list[[x]]$relabel

    if (!file.exists(joint_post_1_files_list[[samp]])) return(NA)

    this_joint_post <- readr::read_tsv(joint_post_1_files_list[[samp]],
                                        show_col_types = FALSE)

    plot_specific_cnv_event(nb          = this_joint_post,
                             annot_df   = sc_annotations_list[[samp]],
                             sample     = samp,
                             cnv_event  = seg_label,
                             healthy_df = healthy_cell_types_df,
                             display_event = relabel)
  }) |>
  setNames(names(chr21_amp_list))


# =============================================================================
# SECTION 6: Assemble and save manuscript figures
# =============================================================================
# Figure panel uses sample IDs and variant set in USER CONFIGURATION.

cat("Saving figures...\n")

# -- Main figure: myeloid mutations tile + chr21 CNV dot plot -----------------
ggsave(
  filename = here::here(plot_dir, "main_fig.myeloid_events.pdf"),
  plot     = myeloid_mutations_plot + CNV_chr21_amp_plot +
             plot_annotation(tag_levels = "a"),
  height   = FIG_MYELOID_HEIGHT * 2,
  width    = FIG_MYELOID_WIDTH  * 2
)

# -- Single-cell annotation panel ---------------------------------------------
# Retrieves specific plots for the focus patient set in USER CONFIGURATION
PAVEPS_IDX_KRAS    <- compartment_mutation_plot_list[[SC_ANNO_SAMPLE_IDX]][[SC_ANNO_VARIANT]]
PAVEPS_MRD_KRAS    <- compartment_mutation_plot_list[[SC_ANNO_SAMPLE_MRD]][[SC_ANNO_VARIANT]]
PAVEPS_IDX_chr21amp <- chr21_amp_plot_list[[SC_ANNO_SAMPLE_IDX]]
PAVEPS_MRD_chr21amp <- chr21_amp_plot_list[[SC_ANNO_SAMPLE_MRD]]

pdf(file = here::here(plot_dir, "main_fig.single_cell_annotations.pdf"),
    width  = FIG_ANNO_WIDTH,
    height = FIG_ANNO_HEIGHT)

healthy_reference_plot + PAVEPS_IDX_KRAS + PAVEPS_MRD_KRAS +
  guide_area() + PAVEPS_IDX_chr21amp + PAVEPS_MRD_chr21amp +
  plot_layout(guides = "collect", ncol = 3) +
  plot_annotation(tag_levels = list(c("c", "d", "e", "f", "g")))

dev.off()

cat("Done.\n")
