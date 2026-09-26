#!/usr/bin/env Rscript
# Draw the manuscript figures from manuscript/figures/data/*.tsv (written by make_figure_data.py).
#
#   Rscript manuscript/figures/plot_figures.R
#
# Writes fig-<name>.pdf (vector, for the journal) and fig-<name>.png (for the HTML render) next to
# this script. Every number drawn or printed comes from the data files; only labels are written here.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
dat <- function(f) read_tsv(file.path(here, "data", f), show_col_types = FALSE)

# One hue, two steps: the light bar is the whole (calls before, genomes carrying), the dark bar the
# part that survives (calls kept, genomes exporting). Text stays in ink, never in the series colour.
C_WHOLE <- "#c5daf3"
C_PART <- "#2a78d6"
INK <- "#1a1a19"
INK_2 <- "#52514e"
GRID <- "#e6e5e0"

theme_fig <- function() {
  theme_minimal(base_size = 8.5, base_family = "Helvetica") +
    theme(
      text = element_text(colour = INK),
      axis.text = element_text(colour = INK_2),
      axis.title = element_text(colour = INK_2, size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 8.5, hjust = 0),
      plot.title.position = "plot",
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.key.size = unit(3, "mm"),
      plot.margin = margin(4, 8, 4, 4)
    )
}

save_fig <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(here, paste0("fig-", name, ".pdf")), p, width = w_mm, height = h_mm,
         units = "mm", device = cairo_pdf)
  ggsave(file.path(here, paste0("fig-", name, ".png")), p, width = w_mm, height = h_mm,
         units = "mm", dpi = 300, device = ragg::agg_png, bg = "white")
  # Vector for the preprint PDF and the HTML; the PNG stays for Word, which renders SVG poorly.
  ggsave(file.path(here, paste0("fig-", name, ".svg")), p, width = w_mm, height = h_mm,
         units = "mm", device = svglite::svglite, bg = "white")
  message("  fig-", name, ".{pdf,png,svg}")
}

# Two overlaid horizontal bars per row with a "part of whole" label, the shared form of two figures.
paired_bars <- function(d, whole_lab, part_lab, xlab) {
  d_long <- bind_rows(
    mutate(d, n = whole, what = whole_lab),
    mutate(d, n = part, what = part_lab)
  ) |> mutate(what = factor(what, levels = c(whole_lab, part_lab)))
  ggplot(d_long, aes(x = n, y = label)) +
    geom_col(data = filter(d_long, what == whole_lab), aes(fill = what), width = 0.72) +
    geom_col(data = filter(d_long, what == part_lab), aes(fill = what), width = 0.72) +
    geom_text(data = d, aes(x = whole, y = label, label = text), inherit.aes = FALSE,
              hjust = -0.12, size = 2.6, colour = INK) +
    scale_fill_manual(values = setNames(c(C_WHOLE, C_PART), c(whole_lab, part_lab))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = xlab, y = NULL) +
    theme_fig()
}

# --- Route calls by lineage, before and after lineage confirmation -------------------------------
routes <- dat("route_calls_by_lineage.tsv")

route_panel <- function(route_name, title, fold_at, other_unit) {
  d <- filter(routes, route == route_name)
  # Lineages with only a couple of calls, all removed, fold into one row so the long tail of
  # unnamed genera does not dominate the panel. Kept rows are never folded.
  small <- d$calls_family_check <= fold_at & d$calls_lineage_policy == 0
  if (sum(small) > 1) {
    d <- bind_rows(
      d[!small, ],
      tibble(route = route_name,
             lineage = sprintf("%d other %s", sum(small), other_unit),
             calls_family_check = sum(d$calls_family_check[small]),
             calls_lineage_policy = 0L)
    )
  }
  d <- d |>
    mutate(whole = calls_family_check, part = calls_lineage_policy,
           text = sprintf("%d of %d", part, whole),
           is_other = grepl("^[0-9]+ other ", lineage),
           italic = !is_other & !grepl("oxidisers|[0-9]|^unnamed", lineage)) |>
    arrange(is_other, whole) |>
    mutate(label = factor(lineage, levels = unique(c(lineage[is_other], lineage[!is_other]))))
  total <- sprintf("%s (%d → %d calls)", title, sum(d$whole), sum(d$part))
  # Genus and order names are italic; placeholder names (DQIP01) and the grouped rows are not.
  it <- setNames(d$italic, d$lineage)
  lab <- function(x) lapply(x, function(v) {
    fam <- sub("^unnamed (\\S+) genera$", "\\1", v)
    if (fam != v) bquote(unnamed ~ italic(.(fam)) ~ genera)
    else if (it[[v]]) bquote(italic(.(v))) else v
  })
  paired_bars(d, "Family-level check", "Lineage policy", NULL) +
    labs(title = total) +
    scale_y_discrete(labels = lab)
}

p_aceto <- route_panel("acetoclastic", "Acetoclastic, by genus", 2, "genera")
p_hydro <- route_panel("hydrogenotrophic", "Hydrogenotrophic, by order", 0, "orders")
p_methyl <- route_panel("methylotrophic", "Methylotrophic, by order", 0, "orders")
n_rows <- function(p) nlevels(p$data$label)
fig_routes <- (p_aceto / p_hydro / p_methyl) +
  plot_layout(heights = c(n_rows(p_aceto), n_rows(p_hydro), n_rows(p_methyl)), guides = "collect") +
  plot_annotation(theme = theme(legend.position = "top", legend.justification = "left")) &
  labs(x = "Genomes with a route call")
save_fig(fig_routes, "route-lineage", 174, 175)

# --- What --contig-level costs, against assembly fragmentation -----------------------------------
contig <- dat("contig_level_loss.tsv") |>
  mutate(lost = 100 * (routed_genome_level - routed_contig_level) / routed_genome_level,
         contigs = factor(contigs, levels = contigs),
         top = sprintf("%.0f %%", lost),
         xlab = sprintf("%s\n%d → %d", contigs, routed_genome_level, routed_contig_level))
fig_contig <- ggplot(contig, aes(x = contigs, y = lost)) +
  geom_col(fill = C_PART, width = 0.62) +
  geom_text(aes(label = top), vjust = -0.5, size = 2.8, colour = INK) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, " %"), expand = expansion(mult = c(0, 0.08))) +
  scale_x_discrete(labels = setNames(contig$xlab, contig$contigs)) +
  labs(x = "Contigs in the genome",
       y = "Methanogens losing every route") +
  theme_fig() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
        axis.title.x = element_text(margin = margin(t = 4)))
save_fig(fig_contig, "contig-loss", 85, 70)

# --- Carries the family versus exports it --------------------------------------------------------
NAMES <- c(
  "HYDROL-PROT" = "Protein hydrolysis", "HYDROL-CARB" = "Carbohydrate hydrolysis",
  "HYDROL-LIP" = "Lipid hydrolysis",
  "ACID-ACE" = "Acetate formation", "ACID-BUT" = "Butyrate formation",
  "ACID-PROP" = "Propionate formation", "ACID-ETA" = "Ethanolamine use",
  "ACID-H2" = "Hydrogen formation",
  "ACET-WL" = "Wood–Ljungdahl", "ACET-SYN-PROP" = "Syntrophic propionate oxidation",
  "MEG-ACETO" = "Acetoclastic", "MEG-HYDRO" = "Hydrogenotrophic", "MEG-METHYL" = "Methylotrophic",
  "MEG-CORE" = "Methanogenesis core"
)
STAGES <- c("hydrolysis", "acidogenesis", "acetogenesis", "methanogenesis")
sec <- dat("secretion_by_module.tsv")
stopifnot(all(sec$module %in% names(NAMES)))
sec <- sec |>
  mutate(whole = genomes_carrying, part = genomes_exporting,
         text = ifelse(part > 0, sprintf("%d of %d (%.1f %%)", part, whole, 100 * part / whole),
                       sprintf("0 of %d", whole)),
         stage = factor(stage, levels = STAGES, labels = tools::toTitleCase(STAGES))) |>
  arrange(desc(stage), whole) |>
  mutate(label = factor(sprintf("%s  %s", NAMES[module], module),
                        levels = unique(sprintf("%s  %s", NAMES[module], module))))
fig_export <- paired_bars(sec, "Carries the family", "Exports it", "Genomes, of 1,401 digester MAGs") +
  facet_grid(stage ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", colour = INK),
        panel.spacing.y = unit(2.5, "mm"))
save_fig(fig_export, "export", 174, 95)
