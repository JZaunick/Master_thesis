#basic DE stats visualisatioin


suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

load("data/A2TEA_finished_potatofix2.RData")
OUT_PDF <- "figures/DE_gene_overview.pdf"
OUT_PNG <- NULL          # set to a path if you also want a PNG

conserved_HOGS <- HOG_level_list$all_species_overview %>% filter(if_all(ends_with("_total"), ~ . > 0)) %>% pull(HOG)


#1. Numerators — DE gene counts per species and category

de_counts <- HOG_DE.a2tea %>%
  filter(significant == "yes") %>%
  group_by(species) %>%
  summarise(
    `All DE genes`                = n(),
    `DE genes in HOGs`            = sum(HOG != "singleton"),
    `DE genes in conserved HOGs`  = sum(HOG %in% conserved_HOGS),
    .groups = "drop"
  )

# Denominators — how many genes exist in each category, per species

gene_universe <- HOG_DE.a2tea %>%
  group_by(species) %>%
  summarise(
    `All DE genes`                = n(),
    `DE genes in HOGs`            = sum(HOG != "singleton"),
    `DE genes in conserved HOGs`  = sum(HOG %in% conserved_HOGS),
    .groups = "drop"
  )

## Long format + percentages, and a printout to sanity-check

CATS <- c("All DE genes", "DE genes in HOGs", "DE genes in conserved HOGs")

plot_df <- de_counts %>%
  pivot_longer(all_of(CATS), names_to = "category", values_to = "n_de") %>%
  left_join(
    gene_universe %>%
      pivot_longer(all_of(CATS), names_to = "category", values_to = "n_total"),
    by = c("species", "category")
  ) %>%
  mutate(
    pct      = 100 * n_de / n_total,
    category = factor(category, levels = CATS),
    # "Solanum_pennellii" -> "S. pennellii"
    sp_label = sub("^([A-Z])[a-z]+_", "\\1. ", species),
    sp_label = factor(sp_label, levels = sort(unique(sp_label)))
  )

cat("\n--- check these denominators before using the % panel ---\n")
print(as.data.frame(plot_df[, c("species", "category", "n_de", "n_total", "pct")]),
      row.names = FALSE, digits = 4)
cat("\nFor reference, the old figure showed C. annuum at 23.3 / 24.0 / 34.5 %",
    "\n(pepper was unaffected by the potato fix, so those should still match).\n\n")



## ---- single faceted figure: rows = measure, columns = species ----------------
cat_short <- c("All DE genes"               = "all",
               "DE genes in HOGs"           = "in HOGs",
               "DE genes in conserved HOGs" = "in conserved HOGs")
DIR_PAL   <- c(up = "#D9822B", down = "#3F7CAC")   # softer orange / blue

plot_df <- plot_df %>%
  mutate(sp_label  = factor(sp_label, levels = sort(unique(sp_label))),
         cat_short = factor(cat_short[as.character(category)], levels = cat_short),
         direction = factor(direction, levels = c("up", "down")))

long <- bind_rows(
  plot_df %>% transmute(sp_label, cat_short, direction,
                        measure = "Number of DE genes", value = n_de,
                        seg_label = format(n_de, big.mark = ",", trim = TRUE)),
  plot_df %>% transmute(sp_label, cat_short, direction,
                        measure = "% DE genes", value = pct,
                        seg_label = NA_character_)
) %>% mutate(measure = factor(measure, levels = c("Number of DE genes", "% DE genes")))

tot_long <- long %>%
  group_by(sp_label, cat_short, measure) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(label = if_else(measure == "% DE genes",
                         sprintf("%.1f%%", value),
                         format(round(value), big.mark = ",", trim = TRUE)))

DE_gene_overview <- ggplot(long, aes(cat_short, value, fill = direction)) +
  geom_col(width = 0.8, position = position_stack()) +
  geom_text(aes(label = seg_label, group = direction), na.rm = TRUE,
            position = position_stack(vjust = 0.5), size = 2.5, colour = "white") +
  geom_text(data = tot_long, aes(cat_short, value, label = label),
            inherit.aes = FALSE, vjust = -0.45, size = 2.7, colour = "grey25") +
  facet_grid(measure ~ sp_label, scales = "free_y", switch = "y") +
  scale_fill_manual(values = DIR_PAL, labels = c("up-regulated", "down-regulated")) +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.14))) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.spacing.x    = unit(0.6, "lines"),
    panel.spacing.y    = unit(1.2, "lines"),
    axis.line.x        = element_line(colour = "grey25", linewidth = 0.4),
    axis.title         = element_blank(),
    axis.text.x        = element_text(angle = 35, hjust = 1, colour = "grey15", size = 8.5),
    strip.placement    = "outside",
    strip.text.x       = element_text(face = "italic", size = 10.5),
    strip.text.y.left  = element_text(angle = 90, size = 10.5),
    legend.position    = "top",
    legend.title       = element_blank()
  )

ggsave("figures/DE_gene_overview.pdf", DE_gene_overview,
       width = 8.4, height = 6.8, device = cairo_pdf, bg = "white")
