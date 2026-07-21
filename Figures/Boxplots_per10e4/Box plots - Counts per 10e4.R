
# ===============================
# REFRACTORED MoA BOXPLOT PIPELINE
# ===============================

library(dplyr)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(ggsignif)
library(stringr)

# ===============================
# CONFIG
# ===============================

scenario_config <- list(
  Hh_only = list(
    include_patterns = c("Hh", "Saline", "Pam2CSK4"),
    exclude_patterns = c("IAV", "NTHi"),
    timepoints = c(2, 6, 24, 48, 72, 96, 144)
  ),
  IAV = list(
    include_patterns = c("Saline", "Hh", "Pam2CSK4"),
    exclude_patterns = c("NTHi"),
    timepoints = c(48, 72, 96, 144)
  ),
  NTHi = list(
    include_patterns = c("NTHi"),
    exclude_patterns = c(),
    timepoints = c(48, 72, 96, 144)
  ),
  Hh_NTHi_only = list(
    include_patterns = c("Hh", "Pam2CSK4", "Saline", "NTHi"),
    exclude_patterns = c("IAV"),timepoints = c(48, 72, 96, 144)
  )
)

# ===============================
# DATA
# ===============================

moa.dat <- read.csv("output/Proportion_CellIDPerSample.csv") %>%
  mutate(
    per10e4 = (Total_idPerSample * 10000) / Total_evtsPerSample,
    Treatment = str_replace_all(Treatment, "Pam1ug", "Pam2CSK4")
  ) %>%
  dplyr::arrange(Timepoint) %>%
  filter(!Cell_ID %in% c("Drop", "Unknown Myeloid"))

annotation_map <- function(x) {
  case_when(
    x == "B-cells" ~ "B-cells",
    x == "moDCs" ~ "Alv Mac",
    x == "CD4+ T-cells" ~ "CD4+ T",
    x == "Ly6C-CD8+ T-cells" ~ "Ly6C-CD8+ T",
    x == "Ly6C+CD8+ T-cells" ~ "Ly6C+CD8+ T",
    x == "Ly6C+ NK cells" ~ "Ly6C+NK",
    x == "Ly6C- NK cells" ~ "Ly6C-NK",
    x == "CD103+CD4+ T-cells" ~ "CD103+CD4+ T",
    x == "Ly6C(hi)CD64(hi) Monocytes" ~ "iMono",
    x == "cDC1" ~ "cDC1",
    x == "Monocytes" ~ "Mono",
    x == "Neutrophils" ~ "Neut",
    x == "CD206+ Macrophages" ~ "CD206+ Mac",
    x == "cDC2" ~ "cDC2",
    x == "pDCs" ~ "pDC",
    x == "CD64+MHCII- Macrophages" ~ "CD64+ Mac"
  )
}

moa.dat <- moa.dat %>% mutate(Annotation = annotation_map(Cell_ID))

df.sig <- read.csv("output/AllPer10e4ComparisonsIncludingSignificance.csv") %>%
  mutate(
    Annotation = annotation_map(Cell_ID),
    Treatment1 = str_replace_all(Treatment1, "Pam1ug", "Pam2CSK4"),
    Treatment2 = str_replace_all(Treatment2, "Pam1ug", "Pam2CSK4")
  )
# ===============================
# COLOURS
# ===============================


my_cols <- c(
  "Saline" = "grey",
  "Saline+IAV" = "grey",
  "Saline+IAV+NTHi" = "grey",
  "Hh" = "#4FAFA7",
  "Hh+IAV" = "#4FAFA7",
  "Hh+IAV+NTHi" = "#4FAFA7",
  "Pam2CSK4" = "#4A4280",
  "Pam2CSK4+IAV" = "#4A4280",
  "Pam2CSK4+IAV+NTHi" = "#4A4280"
)


# ===============================
# FILTER
# ===============================

filter_scenario <- function(df, config, treatment_levels) {


  df <- df %>% filter(Timepoint %in% config$timepoints)

  if (length(config$include_patterns) > 0) {
    df <- df %>% filter(grepl(paste(config$include_patterns, collapse = "|"), Treatment))
  }

  if (length(config$exclude_patterns) > 0) {
    df <- df %>% filter(!grepl(paste(config$exclude_patterns, collapse = "|"), Treatment))
  }

  df <- df %>%
    mutate(Treatment = factor(Treatment, levels = names(treatment_levels)),
           Tissue = factor(Tissue, levels = c("Nasal Tissue", "Lungs", "Spleen"))) %>%
    droplevels()

  df
}

# ===============================
# SIGNIFICANCE
# ===============================

compute_significance <- function(data, sig_df, name) {

  if (nrow(data) == 0) return(NULL)

  rng <- range(data$per10e4, na.rm = TRUE)
  value.range <- diff(rng)

  if (is.na(value.range) || value.range == 0) {
    value.range <- max(rng, na.rm = TRUE) * 0.1 + 1e-6
  }

  df.max <- data %>%
    group_by(Tissue, Timepoint, Treatment) %>%
    summarise(max_val = max(per10e4, na.rm = TRUE), .groups = "drop")

  sig <- sig_df %>%
    filter(
      Annotation == name,
      significance != ""
    ) %>%
    inner_join(df.max,
               by = c("Tissue", "Timepoint", "Treatment1" = "Treatment")) %>%
    rename(g1 = max_val) %>%
    inner_join(df.max,
               by = c("Tissue", "Timepoint", "Treatment2" = "Treatment")) %>%
    rename(g2 = max_val)

  sig <- sig %>%
    filter(!is.na(g1), !is.na(g2))

  if (nrow(sig) == 0) return(NULL)

  sig <- sig %>%
    mutate(max_val = pmax(g1, g2)) %>%
    group_by(Tissue, Timepoint) %>%
    dplyr::arrange(max_val, .by_group = TRUE) %>%
    mutate(
      row = row_number(),
      base = max(max_val, na.rm = TRUE),
      step = 0.15 * value.range,
      y_pos = base + row * step
    ) %>%
    ungroup()

  sig <- sig %>%
    filter(
      is.finite(y_pos),
      y_pos > 0,
      Treatment1 != Treatment2
    )

  if (nrow(sig) == 0) return(NULL)

  sig
}

# ===============================
# PLOT
# ===============================

plot_moa <- function(data, sig_df, name, scenario_name) {

  if (nrow(data) == 0) return(NULL)

  sig <- compute_significance(data, sig_df, name)

  data <- data %>%
    mutate(
      Tissue = trimws(Tissue),
      Tissue = factor(Tissue, c("Nasal Tissue", "Lungs", "Spleen")),
      Timepoint = as.character(Timepoint),
      Timepoint = case_when(
        Timepoint %in% c("2", "2h") ~ "2hrs",
        Timepoint %in% c("6", "6h") ~ "6hrs",
        Timepoint %in% c("24", "24h") ~ "24hrs",
        Timepoint %in% c("48", "48h") ~ "48hrs",
        Timepoint %in% c("72", "72h") ~ "72hrs",
        Timepoint %in% c("96", "96h") ~ "96hrs",
        Timepoint %in% c("144", "144h") ~ "144hrs"
      ),
      Treatment = factor(Treatment, names(my_cols))
    )

  if (!is.null(sig)) {
    sig <- sig %>%
      mutate(
        Tissue = factor(Tissue, c("Nasal Tissue", "Lungs", "Spleen")),
        Timepoint = case_when(
          Timepoint %in% c("2", "2h") ~ "2hrs",
          Timepoint %in% c("6", "6h") ~ "6hrs",
          Timepoint %in% c("24", "24h") ~ "24hrs",
          Timepoint %in% c("48", "48h") ~ "48hrs",
          Timepoint %in% c("72", "72h") ~ "72hrs",
          Timepoint %in% c("96", "96h") ~ "96hrs",
          Timepoint %in% c("144", "144h") ~ "144hrs"
        )
      )
  }

  # -----------------------------
  # STORE EXPORT DATA
  # -----------------------------
  data_export[[paste0(name, "_", scenario_name)]] <<- list(
    cells = data,
    sig   = sig
  )

  # -----------------------------
  # Y LIMITS INCLUDING SIGNIFICANCE
  # -----------------------------
  y_limits <- if (!is.null(sig) && nrow(sig) > 0) {
    range(c(data$per10e4, sig$y_pos), na.rm = TRUE)
  } else {
    range(data$per10e4, na.rm = TRUE)
  }

  tissues <- levels(data$Tissue)

  plots <- lapply(seq_along(tissues), function(i) {

    t <- tissues[i]
    df_t <- data %>% filter(Tissue == t)

    if (nrow(df_t) == 0) return(NULL)

    sig_t <- if (!is.null(sig)) sig %>% filter(Tissue == t) else NULL

    p <- ggplot(df_t, aes(Treatment, per10e4)) +
      geom_point(aes(col = Treatment)) +
      geom_boxplot(outliers = FALSE, fill = "transparent") +
      facet_wrap(~ Timepoint, nrow = 1, drop = FALSE) +
      scale_color_manual(values = my_cols) +
      scale_y_continuous(
        limits = y_limits,
        expand = expansion(mult = c(0.1, 0.1))
      ) +
      theme_classic() +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = -90),
        legend.position = "none"
      ) +
      labs(
        title = t,
        y = if (i == 1) expression("Cells per 10"^4) else NULL
      )

    if (!is.null(sig_t) && nrow(sig_t) > 0) {
      p <- p + geom_signif(
        data = sig_t,
        aes(annotations = significance),
        xmin = sig_t$Treatment1,
        xmax = sig_t$Treatment2,
        y_position = sig_t$y_pos,
        manual = TRUE,
        tip_length = 0.01,
        textsize = 4
      )
    }

    p
  })

  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0) return(NULL)

  combined <- wrap_plots(plots, nrow = 1)

  ggsave(
    file.path(getwd(), "Figures", "Boxplots_per10e4",
              paste0(name, "_SuppFig.png")),
    combined,
    width = 7.5,
    height = 3.7,
    dpi = 600
  )

  print(combined)
}

# ===============================
# DRIVER
# ===============================

data_export <- list()

run_all <- function() {
  for (scenario_name in names(scenario_config)) {

    message("Running: ", scenario_name)


    cfg <- scenario_config[[scenario_name]]
    data_sub <- filter_scenario(moa.dat, cfg, my_cols)

    if (nrow(data_sub) == 0) next


    for (ct in unique(data_sub$Annotation)) {
      plot_moa(
        data_sub %>% filter(Annotation == ct),
        df.sig,
        ct,
        scenario_name
      )
    }
  }
}

library(openxlsx)

bind_with_cell <- function(lst, slot_name) {
  dfs <- lapply(names(lst), function(nm) {
    df <- lst[[nm]][[slot_name]]
    if (is.null(df)) return(NULL)
    if (is.data.frame(df) && nrow(df) > 0) {
      df$Cell <- nm
      df
    } else {
      NULL
    }
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) data.frame() else do.call(rbind, dfs)
}

# ---------------- CELLS ----------------
wb_cells <- createWorkbook()

cells_combined <- bind_with_cell(data_export, "cells")
addWorksheet(wb_cells, "Combined")
writeData(wb_cells, "Combined", cells_combined)

for (nm in names(data_export)) {
  df <- data_export[[nm]]$cells
  if (is.null(df) || nrow(df) == 0) next

  safe_nm <- substr(gsub("[\\/:*?\\[\\]]", "_", nm), 1, 31)
  addWorksheet(wb_cells, safe_nm)
  writeData(wb_cells, safe_nm, df)
}

saveWorkbook(wb_cells, "output/Boxplots_per10e4/cells.xlsx", overwrite = TRUE)

# ---------------- SIGNIFICANCE ----------------
wb_sig <- createWorkbook()

sig_combined <- bind_with_cell(data_export, "sig")
addWorksheet(wb_sig, "Combined")
writeData(wb_sig, "Combined", sig_combined)

for (nm in names(data_export)) {
  df <- data_export[[nm]]$sig
  if (is.null(df) || nrow(df) == 0) next

  safe_nm <- substr(gsub("[\\/:*?\\[\\]]", "_", nm), 1, 31)
  addWorksheet(wb_sig, safe_nm)
  writeData(wb_sig, safe_nm, df)
}

saveWorkbook(wb_sig, "output/Boxplots_per10e4/significance.xlsx", overwrite = TRUE)


run_all()