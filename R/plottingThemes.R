#' @import ggplot2
#' @import viridis
#' @export
# Theme for single plots
promor_theme <- ggplot2::theme_classic() +
  ggplot2::theme(
    panel.border = element_rect(
      fill = NA,
      colour = "grey40",
      size = 0.5
    ),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey40"),
    axis.line = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
#' @export
# Theme for faceted plots
promor_facet_theme <- ggplot2::theme_classic() +
  ggplot2::theme(
    panel.border = element_rect(
      fill = NA,
      colour = "grey40",
      size = 0.5
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "grey40"),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      colour = "grey20",
      hjust = 0.01,
      face = "bold",
      vjust = 0
    )
  )

#' @export
# Pick viridis colors
set_col <- function(palette,
                    n,
                    direction = 1) {
  # Assign palette
  if (palette == "inferno") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.6,
      end = 1,
      direction = direction,
      option = "B"
    )(n)
  } else if (palette == "magma") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.4,
      end = 1,
      direction = direction,
      option = "A"
    )(n)
  } else if (palette == "plasma") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.55,
      end = 1,
      direction = direction,
      option = "C"
    )(n)
  } else if (palette == "cividis") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.2,
      end = 1,
      direction = direction,
      option = "E"
    )(n)
  } else if (palette == "rocket") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.3,
      end = 1,
      direction = direction,
      option = "F"
    )(n)
  } else if (palette == "mako") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.5,
      end = 1,
      direction = direction,
      option = "G"
    )(n)
  } else if (palette == "turbo") {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.3,
      end = 1,
      direction = direction,
      option = "H"
    )(n)
  } else {
    pal_col <- viridis::viridis_pal(
      alpha = 1,
      begin = 0.5,
      end = 1,
      direction = direction,
      option = "D"
    )(n)
  }
  return(pal_col)
}
