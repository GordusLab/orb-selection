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
  wrap_desc = NULL, 
  truncate_desc = NULL, 
  transparent = FALSE
  ) {
  # Load the data
  df <- read.table(
    file = sumfile,
    sep = "\t", quote = "", row.names = NULL, stringsAsFactors = FALSE, header = TRUE
  )
  print(rev)
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

  theme_elements <- list(
    text = element_text(family = "Verdana"),
    axis.text.x = element_text(color = "black", size = axis_text_size),
    axis.text.y = element_text(color = "black", size = axis_text_size),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(color = "black", size = 1, fill = colors[1]),
    strip.text = element_text(size = strip_text_size, face = "bold"),
    plot.title = element_text(size = title_size, hjust = hj)
  )

  if (transparent) {
    theme_elements$rect <- element_rect(fill = "transparent")
  }

  p <- ggbarplot(df,
    x = "Description", y = "Count", fill = "p",
    orientation = "horiz", width = 0.8, ylab = "Gene count"
  ) +
    facet_grid(
      rows = vars(Ontology), scales = "free_y",
      space = "free_y", drop = T, switch = sw
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

  max_count <- suppressWarnings(max(as.numeric(df$Count), na.rm = TRUE))

  if (rev) {
    if (is.finite(max_count) && max_count <= 10) {
      p <- p + scale_y_reverse(breaks = c(10, 5, 0), limits = c(10, 0)) + scale_x_discrete(position = "bottom")
    } else {
      p <- p + scale_y_reverse() + scale_x_discrete(position = "bottom")
    }
  } else if (is.finite(max_count) && max_count <= 10) {
    p <- p + scale_y_continuous(breaks = c(0, 5, 10), limits = c(0, 10))
  }

  if (!is.null(title)) {
    p <- p + theme(plot.title = element_blank())
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
  sumfile = here("results/go_enrichment/relax_relaxed/summary_relax_relaxed_hits.txt"),
  colors = c("#F5EFF9", "#926AB0"),
  figsize = c(25, 24),
  wrap_desc = 45,
  rev = FALSE,
  output_filename = here("figures/figure_3/relax_relaxed_go_barplot.png")
)

go.barplot(
  sumfile = here("results/go_enrichment/relax_intensified/summary_relax_intensified_hits.txt"),
  colors = c("#F5EFF9", "#926AB0"),
  figsize = c(25, 21),
  rev = FALSE,
  output_filename = here("figures/figure_3/relax_intensified_go_barplot.png")
)

go.barplot(
  sumfile = here("results/go_enrichment/busted_ph/summary_busted_ph_hits.txt"),
  colors = c("#F5EFF9", "#926AB0"),
  figsize = c(26, 20),
  rev = FALSE,
  output_filename = here("figures/figure_4/busted_ph_go_barplot.png"),
  axis_text_size = 26,
  axis_title_size = 30,
  strip_text_size = 26
)

go.barplot(
  sumfile = here("results/go_enrichment/busted_ph_rev/summary_busted_ph_rev_hits.txt"),
  colors = c("#F5EFF9", "#926AB0"),
  figsize = c(26, 20),
  wrap_desc = 30,
  rev = FALSE,
  output_filename = here("figures/figure_4/busted_ph_rev_go_barplot.png"),
  axis_text_size = 26,
  axis_title_size = 30,
  strip_text_size = 26
)

# # Log odds ratio analysis
# go.barplot(
#   sumfile = here("results/go_enrichment/odds_ratio_test/duplication_nonorb/summary_duplication_nonorb.txt"),
#   colors = c("#FEE9E7", "#FA8072"),
#   figsize = c(25, 20),
#   title = "Genes more likely to be duplicated in non-orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/figure_5/duplication_nonorb_go_barplot.png"),
#   truncate_desc = 50
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/odds_ratio_test/duplication_orb/summary_duplication_orb.txt"),
#   colors = c("#E0F2DB", "#639B51"),
#   figsize = c(30, 25),
#   title = "Genes more likely to be duplicated in orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/figure_5/duplication_orb_go_barplot.png"),
#   truncate_desc = 50
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/odds_ratio_test/loss_nonorb/summary_loss_nonorb.txt"),
#   colors = c("#DFE9F2", "#4682B4"),
#   figsize = c(25, 15),
#   title = "Genes more likely to be missing in non-orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/figure_5/loss_nonorb_go_barplot.png"),
#   truncate_desc = 50
# )

# go.barplot(
#   sumfile = here("results/go_enrichment/odds_ratio_test/loss_orb/summary_loss_orb.txt"),
#   colors = c("#F9F0D9", "#DAA520"),
#   figsize = c(25, 20),
#   title = "Genes more likely to be missing in orbweavers",
#   rev = FALSE,
#   output_filename = here("figures/figure_5/loss_orb_go_barplot.png"),
#   truncate_desc = 50
# )
