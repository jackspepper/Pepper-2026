# This script assumes the setup chunk in the primary script has been run to define all folders.

# ---- Packages ----
sessionpackages <- c(
  "dplyr", "ggplot2", "tidyr", "cowplot", "grid",
  "RColorBrewer", "MetBrewer", "scales", "tibble",
  "emmeans"
)


invisible(lapply(sessionpackages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package not installed: ", pkg, ". Please install it first."))
  }
  library(pkg, character.only = TRUE)
}))

# ---- Input ----
biod.dat <- read.csv("Figures/Boxplots_Biodistribution/Hh biodistribution.csv") %>%
  as.data.frame() %>%
  dplyr::filter(Group != "Hh+IAV")

# ---- Reshape to long format ----
biod.long <- biod.dat %>%
  tidyr::pivot_longer(
    cols = c("Nasal.Hh.density", "Lung.Hh.Density", "Blood.Hh.Density"),
    names_to = "Tissue",
    values_to = "Hh_Density"
  ) %>%
  dplyr::mutate(
    Tissue_moa = dplyr::case_when(
      Tissue == "Nasal.Hh.density" ~ "Nasal Tissue",
      Tissue == "Lung.Hh.Density"  ~ "Lungs",
      Tissue == "Blood.Hh.Density" ~ "Blood",
      TRUE ~ Tissue
    )
  ) %>%
  mutate(Timepoint = stringr::str_remove(Timepoint, "h"))

# ---- Limits of detection (LOD) ----
# Edit these values if your assay's LODs change.
LOD <- tibble::tibble(
  Tissue_moa = c("Nasal Tissue", "Lungs", "Blood"),
  lod        = c(0.1, 0.02, 0.012),
  lod_half   = lod / 2
)

# ---- Apply half-LOD substitution for values below LOD ----
# IMPORTANT: These adjusted densities are used for BOTH plotting AND statistical testing.
biod.long <- biod.long %>%
  dplyr::left_join(LOD, by = "Tissue_moa") %>%
  dplyr::mutate(
    Density_adj = dplyr::case_when(
      is.na(Hh_Density) ~ NA_real_,
      !is.na(lod) & Hh_Density < lod ~ lod_half,
      TRUE ~ Hh_Density
    )
  )

# ---- Factor levels for consistent ordering ----
my_levels <- c("1", "2", "6", "24", "48", "96")

# Force Timepoint ordering even when some levels are missing:
biod.long <- biod.long %>%
  dplyr::mutate(Timepoint = factor(Timepoint, levels = my_levels))

# ---- Colors ----
my_cols <- c(
  "Naive" = "grey70",
  "Hh"    = "#4FAFA7"
)

# ---- Stats: Linear model + emmeans contrasts ----

library(emmeans)

stat.test <- biod.long %>%
  dplyr::filter(!is.na(Density_adj), !is.na(Group), !is.na(Timepoint)) %>%
  dplyr::mutate(logDensity = log10(Density_adj)) %>%
  dplyr::group_split(Tissue_moa) %>%
  lapply(function(df) {

    tissue_name <- unique(df$Tissue_moa)

    # Fit linear model
    fit <- lm(logDensity ~ Group * Timepoint, data = df)

    # Get estimated marginal means
    em <- emmeans::emmeans(fit, ~ Group | Timepoint)

    # Pairwise comparison (Hh vs Naive per timepoint)
    contr <- emmeans::contrast(em, method = "revpairwise") %>%
      as.data.frame()

    # Clean output
    out <- contr %>%
      dplyr::mutate(
        Tissue_moa = tissue_name,
        Timepoint = as.character(Timepoint),
        p = p.value
      ) %>%
      dplyr::select(Tissue_moa, Timepoint, p)

    # Add y-position for plotting
    ymax <- df %>%
      dplyr::group_by(Timepoint) %>%
      dplyr::summarise(y.max = max(Density_adj, na.rm = TRUE), .groups = "drop")

    out <- dplyr::left_join(out, ymax, by = "Timepoint")

    return(out)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(Tissue_moa) %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p, method = "BH"),
    y.position = y.max * 1.6,
    p.signif = dplyr::case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::ungroup()

# ---- Plot function ----
make_plot <- function(tissue_name, show_y_axis = TRUE) {

  df <- biod.long %>% dplyr::filter(Tissue_moa == tissue_name)
  stat.df <- stat.test %>% dplyr::filter(Tissue_moa == tissue_name)

  # Defensive: compute LOD for tissue
  this_lod <- unique(df$lod)
  this_lod <- this_lod[!is.na(this_lod)]
  if (length(this_lod) == 0) this_lod <- NA_real_

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Timepoint,
      y = Density_adj,
      colour = Group
    )
  ) +
    ggplot2::geom_hline(
      yintercept = this_lod,
      linetype = "dashed"
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = interaction(Timepoint, Group)),
      position = ggplot2::position_dodge(0.75),
      colour = "black",
      outlier.shape = NA
    ) +
    ggplot2::geom_dotplot(
      aes(
        group = interaction(Timepoint, Group),
        fill = Group
      ),
      binaxis = "y",        # stack along y (your density axis)
      stackdir = "center",  # center the stack (nice symmetric look)
      position = ggplot2::position_dodge(width = 0.75),
      dotsize = 0.6,
      colour = "black",
      width = 0.1,
      stackratio = 0.5
    ) +
    ggplot2::scale_fill_manual(values = my_cols) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_log10(
      limits = c(0.005, 300),
      breaks = c(0.01, 0.1, 1, 10, 100),
      labels = c(0.01, 0.1, 1, 10, 100)
    ) +
    ggplot2::labs(
      title = tissue_name,
      x = NULL,
      y = "Hh Density (µg/mL)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "none",
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    ) +
    ggplot2::geom_text(
      data = stat.df %>% dplyr::filter(!is.na(p.signif)),
      ggplot2::aes(
        x = Timepoint,
        y = y.position,
        label = p.signif
      ),
      inherit.aes = FALSE,
      size = 6
    )

  if (!show_y_axis) {
    p <- p + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank()
    )
  }

  return(p)
}

# ---- Build the 3-panel figure ----
p1 <- make_plot("Nasal Tissue", show_y_axis = TRUE)
p2 <- make_plot("Lungs", show_y_axis = FALSE)
p3 <- make_plot("Blood", show_y_axis = FALSE)

plots_row <- cowplot::plot_grid(
  p1, p2, p3,
  nrow = 1,
  rel_widths = c(1.30, 1, 1)
)

xlab <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Timepoint (hrs)",
    fontface = "plain",
    size = 12
  )

final_plot <- cowplot::plot_grid(
  plots_row,
  xlab,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

print(final_plot)

# ---- Save ----
ggplot2::ggsave(
  "Figures/Boxplots_Biodistribution/Fig1C_Biodistribution.png",
  plot = final_plot,
  width = 14.5,
  height = 8,
  units = "cm"
)

# ---- Optional: save stats table ----
# write.csv(stat.test, "Biodistribution_stats_Wilcoxon_BH_FDR.csv", row.names = FALSE)
