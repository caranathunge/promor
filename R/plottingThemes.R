#' @import ggplot2
#' @import viridis
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

set_col <- function(palette){
  #Assign palette
  if (palette == "inferno"){
    pal.col = viridis::inferno(15)
    } else if (palette == "magma"){
      pal.col = viridis::magma(18)
      } else if (palette == "plasma"){
        pal.col = viridis::plasma(134)
        } else if (palette == "cividis"){
          pal.col = viridis::cividis(11)
          } else if (palette == "rocket"){
            pal.col = viridis::rocket(14)
            } else if (palette == "mako"){
              pal.col = viridis::mako(67)
              } else if (palette == "turbo"){
                pal.col = viridis::turbo(80)
                } else {
                  pal.col = viridis::viridis(79)
                  }
return(pal.col)
}
