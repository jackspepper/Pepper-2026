# This script assumes the setup chunk in the primary script has been run to define all folders.

sessionpackages <- (c("dplyr", "ggplot2", "MetBrewer", "RColorBrewer","tidyr"))
lapply(sessionpackages, library, character.only = TRUE)

met.dat <- read.csv("Figures/Boxplots_Biodistribution/MoA MET counts.csv")
met.dat <- met.dat %>% filter(!Treatment == "Pam")
met.dat <- met.dat %>% mutate(Pre_Treatment = case_when(Treatment == "Pam" ~ "Pam",
                                                        Treatment == "Hh" ~ "Hh",
                                                        Treatment == "PBS" ~ "Saline"))

my_cols = c("Saline" = "grey70",
             "Hh" = "#4FAFA7",
             "Pam" = "#4A4280")
my_facet = c("Nasal.Hh.density" = "Nasal Tissue",
             "Lung.Hh.Density" = "Lung",
             "Blood.Hh.Density" = "Blood")

# ---- Prepare data ----
met.dat <- met.dat %>%
  dplyr::mutate(
    Pre_Treatment = factor(Pre_Treatment, levels = c("Hh", "Saline", "Pam"))
  )

pval <- wilcox.test(
  MET.Counts..CFU.mL. ~ Pre_Treatment,
  data = met.dat %>% dplyr::filter(Pre_Treatment %in% c("Hh", "Saline")),
  exact = FALSE
)$p.value

# ---- Plot ----
p <- ggplot(met.dat,
            aes(x = Pre_Treatment,
                y = MET.Counts..CFU.mL.)) +
  
  # LOD line (match style)
  geom_hline(yintercept = 1.92, linetype = "dashed") +
  
  # Boxplot (no outliers, black outline)
  geom_boxplot(
    aes(group = Pre_Treatment),
    colour = "black",
    outlier.shape = NA,
    width = 0.5
  ) +
  
  # Dotplot (match biodistribution style)
  geom_dotplot(
    aes(fill = Pre_Treatment),
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.6,
    stackratio = 0.5,
    colour = "black",
    stroke = 0.4
  ) +
  
  # Colors
  scale_fill_manual(values = my_cols) +
  
  # Labels
  labs(
    y = "Middle Ear Tissue NTHi (CFU/mL)",
    x = NULL
  ) +
  
  # Classic theme + matching tweaks
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  ) +
  
  # Significance as TEXT (instead of brackets)
  geom_text(
    data = data.frame(
      x = 1.5,  # midpoint between Saline (1) and Hh (2)
      y = max(met.dat$MET.Counts..CFU.mL., na.rm = TRUE) * 1.4,
      label = sprintf("p = %.2g", pval)
    ),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  )


# ---- Save ----
ggsave(
  "Figures/Boxplots_Biodistribution/Fig1D_CountsMET.png",
  plot = p,
  height = 6.69,
  width = 4.42,
  units = "cm"
)
