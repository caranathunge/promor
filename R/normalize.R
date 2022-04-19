# Normalize intensities ----------------------------------------------------
#' Normalize intensity data
#' @author Chathurani Ranathunge
#' @description This function normalizes data using a user-specified
#' normalization method.
#' @import limma
#' @export

normalize.data <- function(x,
                      method = "quantile"){
  norm_df <- limma::normalizeBetweenArrays(x,
                                           method = method)
  return(norm_df)
}

# Visualize normalization effects -----------------------------------------
#' Visualize the effect of normaliztion on the data
#' @author Chathurani Ranathunge
#' @description This function helps visualize the impact of normalization with
#' box plots or density plots.
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
norm.plot <- function(normalized,
                      original,
                      type = "box",
                      alpha = 0.25,
                      lwd = 0.1,
                      xlabel = "",
                      ylabel="",
                      xlab.size = 7,
                      ylab.size = 7,
                      palette = "YlGnBu",
                      strip.txt.size = 5,
                      strip.ln.size = 0.5,
                      lab.text.size = 5,
                      save = FALSE,
                      filename = "Norm_plot",
                      filetype = "pdf",
                      dpi = 80,
                      width = 7,
                      height = 7){

  #Pre-prossesing data for plotting
  norm1 <- reshape2::melt(normalized, na.rm = FALSE)
  norm1$normstage <- "After normalization"
  orig1 <- reshape2::melt(original, na.rm =FALSE)
  orig1$normstage <- "Before normalization"
  plot_data<- rbind(orig1, norm1)
  colnames(plot_data) <- c("prot", "sample", "intensity", "normstage")
  plot_data$group<- sapply(strsplit(as.character(plot_data[,"sample"]),'_'),
                           getElement,1)
  plot_data$normstage <- factor(plot_data$normstage,
                                levels = c("Before normalization",
                                           "After normalization"))

  if(type == "density"){
    norm_plot<- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = intensity,
                                         color = sample)) +
      ggplot2::geom_density(lwd = lwd)+
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel)+
      #ggplot2::scale_fill_brewer(palette = palette)+
      ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "none",
            text= element_text(size = lab.text.size),
            axis.line.x = element_line(size = 0.1),
            axis.line.y = element_line(size = 0.1),
            axis.ticks.x = element_line(size = 0.1),
            axis.ticks.y = element_line(size = 0.1),
            axis.title.x = element_text(size = xlab.size),
            axis.title.y= element_text(size = ylab.size),
            strip.text = element_text(size= strip.txt.size),
            strip.background = element_rect(colour = "grey",
                                            size = strip.ln.size))+
      ggplot2::facet_wrap( ~normstage)


#Default: Boxplots

  }else{
    norm_plot<- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = sample,
                                         y = intensity,
                                         fill = group)) +
      geom_boxplot(color = "grey30",
                   alpha = 0.9,
                   outlier.shape = 1,
                   outlier.stroke = 0.1,
                   outlier.size = 0.2,
                   outlier.color = "grey30",
                   lwd = lwd)+
      ggplot2::coord_flip()+
      ggplot2::facet_wrap( ~normstage)+
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel)+
      ggplot2::scale_fill_brewer(palette = palette)+
      ggplot2::theme_classic()+
      ggplot2::theme(text= element_text(size = lab.text.size),
            legend.title = element_blank(),
            axis.line.x = element_line(size = 0.1),
            axis.line.y = element_line(size = 0.1),
            axis.ticks.x = element_line(size = 0.1),
            axis.ticks.y = element_line(size = 0.1),
            axis.title.x = element_text(size = xlab.size),
            axis.title.y= element_text(size = ylab.size),
            strip.text = element_text(size= strip.txt.size),
            strip.background = element_rect(colour = "grey",
                                            size = strip.ln.size))

  }
  if(save == TRUE){
    ggsave(paste0(filename,".", filetype),
           norm_plot,
           dpi = dpi,
           width = width,
           height = height)
    }else{
      return(norm_plot)
    }
  }
