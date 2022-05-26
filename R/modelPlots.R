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
                      file.name = "Predictor_plot",
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
#' @description This function visualizes variable importance in models
#'
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import RColorBrewer
#'
#'
#' @param model.list A \code{model.list} object from performing
#' \code{train_models}.
#' @param ... Additional arguments to be passed on to
#' \code{\link[caret:varImp]{varImp}}.
#' @param type Type of plot to generate. Choices are "bar" or "lollipop."
#' Default is \code{"lollipop."}
#' @param text.size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Color palette for plots.
#' @param nrow The number of rows to print the plots.
#' @param ncol The number of columns to print the plots.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file.name file.name File name to save the plot.
#' Default is \code{"VarImp_plot."}
#' @param file.type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot.width Width of the plot. Default is \code{7}.
#' @param plot.height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details \itemize{\item \code{varimp_plot} produces a list of plots showing
#' variable importance measures calculated from models generated with different
#' machine-learning algorithms.
#' \item Note: Variables are ordered by variable importance in
#' descending order and by default, importance values are scaled to 0 and 100.
#' This can be changed by providing \code{scale = FALSE}. See
#' \code{\link[caret:varImp]{varImp}} for more information.}
#'
#' @return A list of \code{ggplot2} objects.
#' @seealso
#' \itemize{
#' \item \code{train_models}
#' \item \code{\link[caret:varImp]{caret: varImp}}}
#' @examples
#' \dontrun{
#'
#' ##Train models using default settings
#' model_list <- train_models(split_df)
#'
#' ## Lollipop plots
#' varimp_plot(model_list)
#'
#' ## Bar plots
#' varimp_plot(model_list, type = "bar")
#'
#' ## Do not scale variable importance values
#' varimp_plot(model_list, scale  = FALSE)
#'
#'
#'  }
#' @export
varimp_plot <- function(model.list,
                        ...,
                        type = "lollipop",
                        text.size = 10,
                        palette = "YlGnBu",
                        nrow,
                        ncol,
                        save = FALSE,
                        file.name = "VarImp_plot",
                        file.type = "pdf",
                        dpi = 80,
                        plot.width = 7,
                        plot.height = 7){

  #Calculate variable importance with VarImp for each ML algorithm-based
  #model in the list
  vimp <- lapply(model.list,
                 function(x) tryCatch(caret::varImp(x,
                                                    ...),
                                      error = function(e) NULL))


  #Drop Null items from the list
  vimp <- Filter(Negate(is.null), vimp)


  #Extract importance estimates from the list
  #imp <- sapply(vimp, function(x) x['importance'], USE.NAMES = TRUE)

  #Make necessary changes to 'importance data frames' in the list before
  #plotting.
  plot_imp_data <- lapply(vimp,
                          function(x) {x <- x$importance[1];
                          x$protein <- rownames(x);
                          rownames(x) <- NULL ;
                          colnames(x) <- c("Importance", "Protein");
                          x$Importance[is.nan(x$Importance)] <- NA;
                          x = x[rowSums(is.na(x)) == 0,];
                          x})

  #Drop empty data frames before proceeding
  plot_imp_data <- Filter(function(x) dim(x)[1] > 0, plot_imp_data)

  #Make plots based on type

  if (type == "bar"){
    #Make bar plots
    vi_plots <- lapply(plot_imp_data, function (t)
      ggplot2::ggplot(data = t,
                      aes(x = reorder(Protein, Importance),
                          y = Importance,
                          fill = Importance)) +
        ggplot2::geom_bar(stat = "identity")+
        ggplot2::coord_flip()+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "right",
                       legend.direction="vertical",
                       legend.key.height = grid::unit(0.8, "cm"),
                       legend.key.width = grid::unit(0.2, "cm"),
                       legend.title = element_blank(),
                       legend.text = element_text(size = text.size*0.7,
                                                  face = "bold"),
                       axis.text.y = element_text(size = text.size,
                                                  face = "bold"),
                       axis.line.x = element_line(size = 0.1),
                       axis.line.y = element_line(size = 0.1),
                       axis.ticks.x = element_line(size = 0.1),
                       axis.ticks.y = element_line(size = 0.1),
                       axis.title.y= element_text(size = text.size,
                                                  face = "bold"),
                       axis.line = element_line(colour = "grey",
                                                size = 0.5),
                       strip.background = element_blank(),
                       strip.text = element_text(size = text.size,
                                                 hjust = 0.01,
                                                 face = "bold",
                                                 vjust = 0 ),
                       panel.border = element_rect(fill = NA,
                                                   colour = "grey",
                                                   size = 0.5))+
        ggplot2::scale_fill_viridis_c(direction = -1)+
        xlab("")
    )
    #Add plot titles
    vi_plots <- lapply(seq_along(vi_plots), function(i) {
      vi_plots[[i]] + ggtitle(gsub(".importance","", names(vi_plots)[i]))
    })

  }else{
    #Make lollipop plots
    vi_plots <- lapply(plot_imp_data, function (t)
        ggplot2::ggplot(data = t,
                      aes(x = reorder(Protein, Importance),
                          y = Importance,
                          color = Importance)) +
        ggplot2::geom_segment(aes( xend=Protein,
                                   yend = 0),
                              lwd = 2) +
        ggplot2::geom_point( aes(color = Importance),
                             size=text.size,
                             show.legend = FALSE)+
        ggplot2::coord_flip()+
        ggplot2::geom_label(label = round(t$Importance,
                                          digits = 1),
                            fill = NA,
                            label.size = NA,
                            colour = "white")+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "right",
                       legend.direction="vertical",
                       legend.key.height = grid::unit(0.8, "cm"),
                       legend.key.width = grid::unit(0.2, "cm"),
                       legend.title = element_blank(),
                       legend.text = element_text(size = text.size*0.7,
                                                  face = "bold"),
                       axis.text.y = element_text(size = text.size,
                                                  face = "bold"),
                       axis.line.x = element_line(size = 0.1),
                       axis.line.y = element_line(size = 0.1),
                       axis.ticks.x = element_line(size = 0.1),
                       axis.ticks.y = element_line(size = 0.1),
                       axis.title.y= element_text(size = text.size,
                                                  face = "bold"),
                       axis.line = element_line(colour = "grey",
                                                size = 0.5),
                       strip.background = element_blank(),
                       strip.text = element_text(size = text.size,
                                                 hjust = 0.01,
                                                 face = "bold",
                                                 vjust = 0 ),
                       panel.border = element_rect(fill = NA,
                                                   colour = "grey",
                                                   size = 0.5))+
        ggplot2::scale_color_viridis_c(direction  = -1)+
        xlab(""))
    #Add plot titles
    vi_plots <- lapply(seq_along(vi_plots), function(i) {
      vi_plots[[i]] + ggtitle(gsub(".importance","",names(vi_plots)[i]))
    })
  }

  if(save == TRUE){
    ggplot2::ggsave(paste0(file.name,".", file.type),
                    marrangeGrob(grobs = vi_plots,
                                 nrow = nrow,
                                 ncol = ncol,
                                 top=""),
                    dpi = dpi)
    grid.arrange(grobs = vi_plots, newpage = TRUE)

  }else{

    grid.arrange(grobs = vi_plots, newpage = TRUE)

  }


}
