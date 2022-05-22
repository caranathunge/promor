# Predictor plots ----------------------------------------------------------------
#' Visualize predictor (protein) variation among conditions
#' @description This function visualizes protein intensity differences among
#' conditions with plots.
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import RColorBrewer
#'
#'
#' @param df A \code{model.df} object from performing \code{pre_process}.
#' @param type Type of plot to generate. Choices are "box" or "density." Default
#' is \code{"box."}
#' @param text.size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Color palette for box plots. Default is \code{"YlGnBu."}
#' @param nrow The number of rows to print the plots.
#' @param ncol The number of columns to print the plots.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file.name file.name File name to save the plot.
#' Default is \code{"Predictor_plot."}
#' @param file.type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot.width Width of the plot. Default is \code{7}.
#' @param plot.height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details This function visualizes condition-wise differences in protein
#' intensity using boxplots and/or density plots.
#'
#' @return A \code{ggplot2} object
#' @seealso
#' \itemize{
#' \item \code{pre_process}, \code{rem_predictors}}
#' @examples
#' \dontrun{
#'
#'  ## Box plots
#'  predictor_plot(model_df, type = "box", nrow = 2, ncol = 4)
#'
#'  ## Density plots
#'  predictor_plot(model_df, type = "density", nrow = 2, ncol = 4)
#'
#'  }
#' @export
predictor_plot <- function(df,
                      type = "box",
                      text.size = 10,
                      palette = "YlGnBu",
                      nrow,
                      ncol,
                      save = FALSE,
                      file.name = "Norm_plot",
                      file.type = "pdf",
                      dpi = 80,
                      plot.width = 7,
                      plot.height = 7){

  #Create plot data
  modeldf_melted <- reshape2::melt(df)
  colnames(modeldf_melted) <- c("Condition", "Protein", "Intensity")

  if (type == "density"){
    pred_plot <- ggplot2::ggplot(modeldf_melted, aes(x=Intensity,
                                                     color=Condition)) +
      ggplot2::geom_density(alpha = 0.25,
                            lwd = text.size * 0.04)+
      ggplot2::scale_color_brewer(palette=palette)+
      ggplot2::facet_wrap(~Protein,
                                   scale = "free",
                                   nrow = nrow,
                                   ncol = ncol)+
      ggplot2::xlab("") +
      ggplot2::ylab("")+
      ggplot2::theme_classic()+
      ggplot2::theme(panel.border = element_rect(fill = NA,
                                                 colour = "grey",
                                                 size = 0.5),
                     legend.title = element_blank(),
                     axis.ticks = element_line(colour = "grey"),
                     axis.title.x = element_text(size = text.size),
                     axis.title.y = element_text(size = text.size),
                     axis.text = element_text(size = text.size/2),
                     axis.line = element_line(colour = "grey",
                                              size = 0.5),
                     legend.text = element_text(size = text.size),
                     strip.background = element_blank(),
                     strip.text = element_text(size = text.size,
                                               hjust = 0.01,
                                               face = "bold",
                                               vjust = 0 ))

    return(pred_plot)
  }else{
    pred_plot <- ggplot2::ggplot(modeldf_melted,
                               ggplot2::aes(x = Condition,
                                            y = Intensity,
                                            fill = Condition)) +
      geom_boxplot(color = "grey30",
                   alpha = 0.7,
                   outlier.shape = 1,
                   outlier.stroke = 0.1,
                   outlier.size = text.size * 0.04,
                   outlier.color = "grey30",
                   lwd = text.size * 0.02)+
      ggplot2::facet_wrap( ~Protein, scales = "free",
                           nrow = nrow,
                           ncol = ncol)+
      ggplot2::xlab("") +
      ggplot2::ylab("")+
      ggplot2::scale_color_brewer(palette = palette)+
      ggplot2::theme_classic()+
      ggplot2::theme(panel.border = element_rect(fill = NA,
                                               colour = "grey",
                                               size = 0.5),
                   text = element_text(size = text.size),
                   legend.title = element_blank(),
                   legend.position = "none",
                   axis.line.x = element_line(size = 0.1),
                   axis.line.y = element_line(size = 0.1),
                   axis.ticks.x = element_line(size = 0.1),
                   axis.ticks.y = element_line(size = 0.1),
                   axis.title.x = element_text(size = text.size),
                   axis.title.y= element_text(size = text.size),
                   axis.line = element_line(colour = "grey",
                                            size = 0.5),
                   strip.background = element_blank(),
                   strip.text = element_text(size = text.size,
                                             hjust = 0.01,
                                             face = "bold",
                                             vjust = 0 ))
    return(pred_plot)
  }
  if(save == TRUE){
    ggsave(paste0(file.name,".", file.type),
           pred_plot,
           dpi = dpi,
           width = plot.width,
           height = plot.height)
  }


}
#-------------------------------------------------------------------------------
#' Variable importance plot
#'
#plot_data <- varImp(rfFit)
#plot_data <- plot_data$importance
#plot_data$protein <- rownames(plot_data)
#rownames(plot_data) <- NULL
#names(plot_data) <- c('importance', 'protein')
#plot_data
#ggplot( data = plot_data, aes(x=reorder(protein, importance), y=importance)) +
 # geom_segment( aes(xend=protein, yend=0)) +
  #geom_point( size=4, color="orange") +
  #coord_flip() +
  #theme_bw() +
  #xlab("")
