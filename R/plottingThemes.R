#' @import ggplot2
promor_theme <- ggplot2::theme_classic(base_size = 12,
                                       base_family = "Helvetica")+
  ggplot2::theme(panel.border = element_rect(fill = NA,
                                             colour = "grey",
                                             size = 0.5),
                 legend.title = element_blank(),
                 axis.ticks = element_line(colour = "grey"),
                 axis.line = element_line(colour = "grey",
                                          size = 0.5))
promor_facet_theme <- ggplot2::theme_classic(base_size = 12,
                                       base_family = "Helvetica")+
  ggplot2::theme(panel.border = element_rect(fill = NA,
                                             colour = "grey",
                                             size = 0.5),
                 legend.title = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y= element_blank(),
                 axis.ticks = element_line(colour = "grey"),
                 axis.line = element_line(colour = "grey",
                                          size = 0.5),
                 strip.background = element_blank(),
                 strip.text = element_text(colour = "grey20",
                                           hjust = 0.01,
                                           face = "bold",
                                           vjust = 0 ))
