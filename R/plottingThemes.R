#' @import ggplot2
#' @export
promor_theme <- ggplot2::theme_classic(base_size = 12,
                                       base_family = "sans")+
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


arctic <- c("#0F1D2B", "#17456B","#419FB7", "#F5F6EE", "#FFB000")

sunflower <- c("#004358", "#1F8A70","#BEDB39", "#FFE11A", "#FD7400")

moss <- c("#D9C5C9", "#0F261E","#1A4032", "#688C7B", "#94A681")

islander <- c("#5E0042", "#2C2233","#005869", "#00856A", "#8DB500")

chesapeake <- c("#FFBC67", "#DA727E","#AC6C82", "#685C79", "#455C7B")

