library(forcats)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(here)

# Core wrapping function
wrap.it <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

# Call this function with a list or vector
wrap.labels <- function(x, len) {
  if (is.list(x)) {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

# Core truncation function
truncate.it <- function(x, len) {
  sapply(x, function(y) {
    if (nchar(y) > len) {
      paste0(substr(y, 1, len), "...")
    } else {
      y
    }
  }, USE.NAMES = FALSE)
}

# Call this function with a list or vector to truncate descriptions
truncate.labels <- function(x, len) {
  if (is.list(x)) {
    lapply(x, truncate.it, len)
  } else {
    truncate.it(x, len)
  }
}

# Collapse the count range (split[1], split[2]] down to a gap of `gap` units so one
# large bar doesn't squash the rest
.split_counts <- function(count, split, gap) {
  ifelse(count > split[2], count - (split[2] - split[1]) + gap, count)
}


# Function to create bar plots for GO enrichment results
go.barplot <- function(
  sumfile, 
  colors, 
  figsize, 
  rev, 
  output_filename, 
  title = NULL,
  title_size = 24, 
  axis_text_size = 18,
  axis_title_size = 20,
  strip_text_size = 20,
  desc_text_size = 16,
  wrap_desc = NULL, 
  truncate_desc = NULL, 
  transparent = FALSE,
  exclude = NULL,          # GO term descriptions to drop from the plot
  split_axis = NULL,       # c(from, to): cut the gene-count axis between these values
  split_gap = 7,           # width of the break, in gene-count units
  split_tick_interval = 10
  ) {
  # Load the data
  df <- read.table(
    file = sumfile,
    sep = "\t", quote = "", row.names = NULL, stringsAsFactors = FALSE, header = TRUE
  )
  omitted <- filter(df, Description %in% exclude)
  df <- filter(df, !Description %in% exclude)
  df <- transform(
    df,
    "Enrichment score" = -log10(p)
  )

  if (is.numeric(wrap_desc)) {
    df["Description"] <- lapply(df["Description"], wrap.labels, len = wrap_desc)
  } else if (is.numeric(truncate_desc)) {
    df["Description"] <- lapply(df["Description"], truncate.labels, len = truncate_desc)
  }

  df$Description <- as.character(df$Description)

  df <- mutate(df, Description = fct_reorder(Description, df$`Enrichment score`))

  if (rev) {
    sw <- NULL
    hj <- 1
  } else {
    sw <- "y"
    hj <- 0
  }

  if (rev && !is.null(split_axis)) {
    stop("split_axis is not supported with rev = TRUE")
  }

  df$.display_count <- if (is.null(split_axis)) df$Count else .split_counts(df$Count, split_axis, split_gap)

  theme_elements <- list(
    text = element_text(family = "Verdana"),
    axis.text.x = element_text(color = "black", size = axis_text_size),
    axis.text.y = element_text(color = "black", size = desc_text_size),
    legend.position = "none",
    axis.title.x = element_text(size = axis_title_size),
    axis.title.y = element_blank(),
    # With a split axis the border is drawn by panel.background so that the
    # white patch behind the "//" (a layer) sits on top of it; panel.border
    # would be drawn over the patch.
    panel.background = element_rect(
      fill = "white", linewidth = 1,
      color = if (is.null(split_axis)) NA else "black"
    ),
    panel.grid = element_blank(),
    panel.border = if (is.null(split_axis)) {
      element_rect(color = "black", fill = NA, linewidth = 1)
    } else {
      element_blank()
    },
    strip.background = element_rect(color = "black", linewidth = 1, fill = colors[1]),
    strip.text = element_text(size = strip_text_size, face = "bold"),
    plot.title = element_text(size = title_size, hjust = hj)
  )

  if (transparent) {
    theme_elements$rect <- element_rect(fill = "transparent")
  }

  p <- ggbarplot(df,
    x = "Description", y = ".display_count", fill = "p",
    orientation = "horiz", width = 0.8, ylab = "Gene count"
  ) +
    facet_grid(
      rows = vars(Ontology), scales = "free_y",
      space = "free_y", drop = TRUE, switch = sw
    ) +
    scale_fill_gradient(
      high = colors[1],
      low = colors[2],
      name = "p-value",
      limits = c(0, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04)
    ) +
    do.call(theme, theme_elements) +
    scale_x_discrete(position = "top")

  if (!is.null(split_axis)) {
    # Ticks below the break are labelled as-is; ticks above it are shifted down to
    # their display position but keep their true count as the label
    upper_ticks <- seq(
      0, ceiling(max(df$Count) / split_tick_interval) * split_tick_interval,
      by = split_tick_interval
    )
    upper_ticks <- upper_ticks[upper_ticks > split_axis[2]]
    tick_labels <- c(seq(0, split_axis[1], by = split_tick_interval), upper_ticks)

    gap_mid <- split_axis[1] + split_gap / 2

    # "//" across the count axis at the bottom edge of the last facet panel
    slashes <- data.frame(Ontology = tail(sort(unique(df$Ontology)), 1))

    p <- p +
      # White block through each bar that crosses the break, wider than the bar
      # so it also covers the outline
      geom_tile(
        data = filter(df, Count > split_axis[2]),
        aes(x = Description, y = gap_mid),
        inherit.aes = FALSE, fill = "white", width = 0.9, height = split_gap
      ) +
      geom_hline(
        yintercept = gap_mid, color = "grey35",
        linewidth = 0.6, linetype = "dotted"
      ) +
      geom_label(
        data = slashes,
        aes(x = -Inf, y = gap_mid, label = "//"),
        inherit.aes = FALSE, size = 1.25 * axis_text_size / .pt, fontface = "bold",
        fill = "white", border.colour = NA, label.padding = unit(0.05, "lines")
      ) +
      scale_y_continuous(
        breaks = .split_counts(tick_labels, split_axis, split_gap),
        labels = tick_labels
      ) +
      expand_limits(y = .split_counts(max(tick_labels), split_axis, split_gap)) +
      coord_flip(clip = "off")
  }

  if (rev) {
    p <- p + scale_y_reverse() + scale_x_discrete(position = "bottom")
  }

  # Note the excluded terms in the top-right corner of their own ontology's panel
  if (nrow(omitted) > 0) {
    omitted_note <- omitted %>%
      group_by(Ontology) %>%
      summarise(label = paste0(
        "Omitted: ", Description, " (", Count, " genes)", collapse = "\n"
      ))
    p <- p + geom_label(
      data = omitted_note,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = 1.02, vjust = 1.3,
      size = 0.8 * desc_text_size / .pt, family = "Verdana", fontface = "italic",
      fill = "white"
    )
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  bg_setting <- if (transparent) "transparent" else "white"

  ggsave(output_filename,
    width = figsize[1],
    height = figsize[2],
    units = "cm",
    bg = bg_setting,
    dpi = 600
  )

}

# Hyphy analysis
go.barplot(
  sumfile = here("results/go_enrichment/relaxed/summary_relaxed.tsv"),
  title = expression("Enriched GO terms, " * bold("relaxed") * " genes"),
  colors = c("#F5EFF9", "#926AB0"),
  figsize = c(30, 26),
  wrap_desc = 35,
  desc_text_size = 22,
  title_size = 28,
  rev = FALSE,
  exclude = "membrane",
  transparent = TRUE,
  output_filename = here("figures/figure_2/relax_relaxed_go_barplot.png")
)

go.barplot(
  sumfile = here("results/go_enrichment/intensified/summary_intensified.tsv"),
  title = expression("Enriched GO terms, " * bold("intensified") * " genes"),
  colors = c("#F5EFF9", "#926AB0"),
  figsize = c(30, 26),
  wrap_desc = 35,
  desc_text_size = 22,
  title_size = 28,
  rev = FALSE,
  exclude = "nucleus",
  split_axis = c(55, 265),
  transparent = TRUE,
  output_filename = here("figures/figure_2/relax_intensified_go_barplot.png")
)

# go.barplot(
#   sumfile = here("results/go_enrichment/busted_ph_orb/summary_busted_ph_orb.tsv"),
#   colors = c("#F5EFF9", "#926AB0"),
#   figsize = c(20, 16),
#   rev = FALSE,
#   output_filename = here("figures/figure_4/busted_ph_orb_go_barplot.png"),
#   wrap_desc = 33,
#   desc_text_size = 20,
#   axis_text_size = 26,
#   axis_title_size = 30,
#   strip_text_size = 26
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/busted_ph_non_orb/summary_busted_ph_non_orb.tsv"),
#   colors = c("#F5EFF9", "#926AB0"),
#   figsize = c(20, 16),
#   desc_text_size = 20,
#   wrap_desc = 30,
#   rev = FALSE,
#   output_filename = here("figures/figure_4/busted_ph_non_orb_go_barplot.png"),
#   axis_text_size = 26,
#   axis_title_size = 30,
#   strip_text_size = 26
# )

# # Log odds ratio analysis
# go.barplot(
#   sumfile = here("results/go_enrichment/dup_bg/summary_dup_bg.tsv"),
#   colors = c("#FEE9E7", "#FA8072"),
#   figsize = c(28,28),
#   title = "Genes more likely to be duplicated in non-orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/supplementary/duplication_nonorb_go_barplot.png"),
#   wrap_desc = 35
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/dup_fg/summary_dup_fg.tsv"),
#   colors = c("#E0F2DB", "#639B51"),
#   figsize = c(30, 25),
#   title = "Genes more likely to be duplicated in orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/supplementary/duplication_orb_go_barplot.png"),
#   wrap_desc = 35
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/loss_bg/summary_loss_bg.tsv"),
#   colors = c("#DFE9F2", "#4682B4"),
#   figsize = c(25, 15),
#   title = "Genes more likely to be missing in non-orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/supplementary/loss_nonorb_go_barplot.png"),
#   wrap_desc = 35
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/loss_fg/summary_loss_fg.tsv"),
#   colors = c("#F9F0D9", "#DAA520"),
#   figsize = c(28, 30),
#   title = "Genes more likely to be missing in orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/supplementary/loss_orb_go_barplot.png"),
#   wrap_desc = 35
# )

# # PhyloGLM

# go.barplot(
#   sumfile = here("results/go_enrichment/pglm_non_orb/summary_pglm_non_orb.tsv"),
#   colors = c("#F5EFF9", "#926AB0"),
#   rev = FALSE,
#   output_filename = here("figures/supplementary/phyloglm_bg_ptep_go_barplot.png"),
#   wrap_desc = 35,
#   figsize = c(25, 20),
#   title = "PhyloGLM genes, copy number ~ non-orb-weaving"
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/pglm_orb/summary_pglm_orb.tsv"),
#   colors = c("#F5EFF9", "#926AB0"),
#   rev = FALSE,
#   output_filename = here("figures/supplementary/phyloglm_fg_udiv_go_barplot.png"),
#   wrap_desc = 35,
#   figsize = c(25, 20),
#   title = "PhyloGLM genes, copy number ~ orb-weaving"
# )
