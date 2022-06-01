# Predictor plots ----------------------------------------------------------------
#' Visualize predictor (protein) variation among conditions
#' @description This function visualizes protein intensity differences among
#' conditions with plots.
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import viridis
#'
#'
#' @param model.df A \code{model.df} object from performing \code{pre_process}.
#' @param type Type of plot to generate. Choices are "box" or "density." Default
#' is \code{"box."}
#' @param text.size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis: scale_color_viridis]{scale_color_viridis}}
#' for available options.
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
predictor_plot <- function(model.df,
                      type = "box",
                      text.size = 10,
                      palette = "viridis",
                      nrow,
                      ncol,
                      save = FALSE,
                      file.name = "Predictor_plot",
                      file.type = "pdf",
                      dpi = 80,
                      plot.width = 7,
                      plot.height = 7){

  #Create plot data
  modeldf_melted <- reshape2::melt(model.df)
  colnames(modeldf_melted) <- c("Condition", "Protein", "Intensity")



  #Make density plots
  if (type == "density"){
    pred_plot <- ggplot2::ggplot(modeldf_melted, aes(x = Intensity,
                                                     color = Condition)) +
      ggplot2::geom_density(alpha = 0.25,
                            lwd = text.size * 0.07)+
      viridis::scale_color_viridis(discrete = TRUE,
                                   option = palette,
                                   begin = 0.3,
                                   end = 0.7)+
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


  }else{
    #Make box plots
    pred_plot <- ggplot2::ggplot(modeldf_melted,
                               ggplot2::aes(x = Condition,
                                            y = Intensity,
                                            col = Condition)) +
      geom_boxplot(alpha = 0.7,
                   outlier.shape = 1,
                   outlier.stroke = 0.1,
                   outlier.size = text.size * 0.1,
                   outlier.color = "grey50",
                   lwd = text.size * 0.07)+
      viridis::scale_color_viridis(discrete = TRUE,
                                   option = palette,
                                   begin = 0.3,
                                   end = 0.7)+
      ggplot2::facet_wrap( ~Protein,
                           scales = "free",
                           nrow = nrow,
                           ncol = ncol)+
      ggplot2::xlab("") +
      ggplot2::ylab("")+
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

  }
  if(save == TRUE){
    ggplot2::ggsave(paste0(file.name,".", file.type),
           pred_plot,
           dpi = dpi,
           width = plot.width,
           height = plot.height)
  }

  return(pred_plot)
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
#' @importFrom viridis scale_color_viridis
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
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis: scale_color_viridis]{scale_color_viridis}}
#' for available options.
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
                        palette = "viridis",
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
                       legend.text = element_text(size = text.size * 0.7,
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
        viridis::scale_fill_viridis(option = palette,
                                     direction = -1)+
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
                              lwd = text.size * 0.2) +
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
                       legend.text = element_text(size = text.size * 0.7,
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
          viridis::scale_color_viridis(option = palette,
                                       direction = -1)+
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
#-------------------------------------------------------------------------------
#' Model performance plot
#' @description This function visualizes model performance
#'
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import caret
#' @import viridis
#'
#'
#' @param model.list A \code{model.list} object from performing
#' \code{train_models}.
#' @param type Type of plot to generate. Choices are "box" or "dot."
#' Default is \code{"box."} for boxplots.
#' @param text.size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis: scale_color_viridis]{scale_color_viridis}}
#' for available options.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file.name file.name File name to save the plot.
#' Default is \code{"Performance_plot."}
#' @param file.type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot.width Width of the plot. Default is \code{7}.
#' @param plot.height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details \itemize{\item \code{performance_plot} uses resampling results from
#' models included in the \code{model.list} to generate plots of model
#' performance.
#' \item The default metrics used for classification based models are "Accuracy"
#' and "Kappa."
#' \item These metric types can be changed by providing additional arguments to
#' the \code{train_models} function. See \code{\link[caret: train]{train}} and
#' \code{\link[caret: trnControl]{trnControl}} for more information.}
#'
#' @return A \code{ggplot2} object.
#' @seealso
#' \itemize{
#' \item \code{train_models}
#' \item\code{\link[caret: resample]{caret: resamples}}
#' \item \code{\link[caret:train]{caret: train}}
#' \item \code{\link[caret:trnControl]{caret: trnControl}}
#' }
#' @examples
#' \dontrun{
#'
#' ##Train models using default settings
#' model_list <- train_models(split_df)
#'
#' ## Generate box plots to visualize performance of different methods
#' performance_plot(model_list)
#'
#' ## Generate dot plots
#' performance_plot(model_list, type = "dot")
#'  }
#' @export
performance_plot <- function(model.list,
                             type = "box",
                             text.size = 10,
                             palette = "viridis",
                             save = FALSE,
                             file.name = "Performance_plot",
                             file.type = "pdf",
                             plot.width = 7,
                             plot.height = 7,
                             dpi = 80){


  #Extract resample data from the model.list object
  resample_data <- resamples(model.list)

  #Convert the performance data into long-form format for plotting.
  plot_data <- reshape2::melt(resample_data$values)

  #Extract method and metric information from the variable column and add them
  #to new columns.
  plot_data$method <- sapply(strsplit(as.character(plot_data$variable),
                                      '~'),
                             "[", 1)
  plot_data$metric <- sapply(strsplit(as.character(plot_data$variable),
                                      '~'),
                             "[", 2)

  #Remove the variable column as it is no longer needed.
  plot_data$variable <- NULL

  #Assign palette color
  pal_col <- set_col(palette, 1)

  #Make plots
  if(type == "dot"){
    #Make dot plots
    perform_plot <- ggplot2::ggplot(plot_data, aes(x = method,
                                                   y = value
                                                   )) +
      ggplot2::stat_summary(fun.data = "mean_se",
                   lwd = text.size * 0.1,
                   color = pal_col)+
      ggplot2::stat_summary(fun = mean,
                   geom = "pointrange",
                   size = text.size * 0.1,
                   pch = 16,
                   position = "identity",
                   color = pal_col)+

      ggplot2::coord_flip()+
      ggplot2::facet_wrap( ~metric, scales = "free")+
      ggplot2::xlab("") +
      ggplot2::ylab("")+
      promor::promor_theme+
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
                     axis.title.y = element_text(size = text.size),
                     axis.line = element_line(colour = "grey",
                                              size = 0.5),
                     strip.background = element_blank(),
                     strip.text = element_text(size = text.size,
                                               hjust = 0.01,
                                               face = "bold",
                                               vjust = 0 ))


    }else{
    #Make box plots
    perform_plot <- ggplot2::ggplot(plot_data,
                                    ggplot2::aes(x = method,
                                                 y = value
                                                 )) +
      geom_boxplot(color = pal_col,
                   width = text.size * 0.03,
                   alpha = 0.7,
                   outlier.shape = 1,
                   outlier.stroke = 0.1,
                   outlier.size = text.size * 0.1,
                   outlier.color = "grey50",
                   lwd = text.size * 0.1)+
      ggplot2::facet_wrap( ~metric,
                           scales = "free")+
      ggplot2::coord_flip()+
      ggplot2::xlab("") +
      ggplot2::ylab("")+
      promor::promor_theme+
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
                     axis.title.y = element_text(size = text.size),
                     axis.line = element_line(colour = "grey",
                                              size = 0.5),
                     strip.background = element_blank(),
                     strip.text = element_text(size = text.size,
                                               hjust = 0.01,
                                               face = "bold",
                                               vjust = 0 ))

    }



  if(save == TRUE){
    ggplot2::ggsave(paste0(file.name,".", file.type),
           perform_plot,
           dpi = dpi,
           width = plot.width,
           height = plot.height)
  }
  return(perform_plot)
}


#-------------------------------------------------------------------------------
#' ROC plot
#' @description This function generates Receiver Operator Characteristic (ROC)
#' curves to evaluate models
#'
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import caret
#' @importFrom pROC roc
#' @importFrom viridis scale_color_viridis
#'
#'
#' @param probability.list A \code{probability.list} object from performing
#' \code{test_models} with \code{type = "prob"}.
#' @param split.df A \code{split.df} object from performing \code{split_df}
#' @param ... Additional arguments to be passed on to
#' \code{\link[pROC: roc]{roc}}.
#' @param multiple.plots Logical. If \code{FALSE} plots all ROC curves
#' representing algorithms included in the \code{probability.list} in a single
#' plot.
#' @param text.size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis: scale_color_viridis]{scale_color_viridis}}
#' for available options.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file.name file.name File name to save the plot.
#' Default is \code{"ROC_plot."}
#' @param file.type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot.width Width of the plot. Default is \code{7}.
#' @param plot.height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details \itemize{\item \code{roc_plot} first uses probabilities generated
#' during \code{test_models} to build a ROC object.
#' \item Next, relevant information is extracted from the ROC object to
#' plot the ROC curves.}
#'
#' @return A \code{ggplot2} object.
#' @seealso
#' \itemize{
#' \item \code{test_models}
#' \item\code{\link[pROC: roc]{pROC: roc}}
#' }
#' @examples
#' \dontrun{
#'
#' ##Test models
#' model_probs <- test_models(model.list = model_fit, split.df = split_df,
#'                           type = "prob")
#'
#' ## Plot all ROC curves in one plot
#' roc_plot(model_probs, multiple.plots = FALSE)
#'
#' ## Plot ROC curves separately for each algorithm
#' roc_plot(model_probs)
#'  }
#' @export
roc_plot <- function(probability.list,
                     split.df,
                     ... ,
                     multiple.plots = TRUE,
                     text.size = 10,
                     palette = "viridis",
                     save = FALSE,
                     file.name = "ROC_plot",
                     file.type = "pdf",
                     plot.width = 7,
                     plot.height = 7,
                     dpi = 80){

  #Generate ROC objects and add them to a list
  roc_list <- lapply(probability.list,
                     function(x) pROC::roc(x[[1]],
                                           response = split.df$test$Condition,
                                           percent = TRUE,
                                           ...))


  #Extract auc values and create a dataframe for plot labels
  auc_list <- lapply(roc_list, function(x) round(x$auc, 1))
  auc_labels <- as.data.frame(do.call("rbind", auc_list))
  auc_labels$method <- rownames(auc_labels)
  rownames(auc_labels) <- NULL
  names(auc_labels) <- c("auc", "method")

  #Extract sensitivities and specificities from ROC object and make a data frame
  # for plotting
  sens_list <- lapply(roc_list, function(x) x$sensitivities)
  sens <- do.call("rbind",sens_list)
  sens_melted <- reshape2::melt(sens)
  names(sens_melted) <- c("method", "number", "sensitivity")

  spec_list <- lapply(roc_list, function(x) x$specificities)
  spec <- do.call("rbind", spec_list)
  spec_melted <- reshape2::melt(spec)
  names(spec_melted) <- c("method", "number", "specificity")

  #Merge and create the final data frame for plotting
  roc_plotdata <- merge(sens_melted, spec_melted)

  #Order data frame by sensitivity. Important for geom_step
  roc_plotdata <- roc_plotdata[order(roc_plotdata$method,
                                     roc_plotdata$sensitivity),]

  #Make ROC curves

  rocplots <- ggplot2::ggplot(roc_plotdata, aes(x = specificity,
                                        y = sensitivity,
                                        colour = method)) +
    ggplot2::geom_step(direction = "hv", lwd = text.size * 0.1)+
    ggplot2::scale_x_reverse(lim=c(100, 0),)+
    ggplot2::xlab("Specificity (%)") +
    ggplot2::ylab("Sensitivity (%)")+
    viridis::scale_color_viridis(discrete = TRUE,
                                 option = palette,
                                 begin = 0.3,
                                 end = 0.7)+
    ggplot2::theme_bw()+
    ggplot2::geom_abline(intercept = 100,
                slope = 1,
                color='grey60',
                linetype = 2,
                show.legend = FALSE)


  if(multiple.plots == FALSE){
  rocplots1 <- rocplots + ggplot2::theme(legend.title = element_blank())

  }else{
  rocplots1 <- rocplots +
    ggplot2::theme(legend.position =  "none",
            strip.background = element_blank(),
            strip.text = element_text(size = text.size,
                                      hjust = 0.01,
                                      face = "bold",
                                      vjust = 0 ))+
    ggplot2::facet_wrap(~method)+
    ggplot2::geom_label(data = auc_labels,
                 aes(label = paste0("AUC: ", auc, "%")),
                 x = Inf,
                 y = -Inf,
                 hjust = 1.2,
                 vjust = -1,
                 inherit.aes = FALSE,
                 label.size = NA)
  }
  if(save == TRUE){
    ggplot2::ggsave(paste0(file.name,".", file.type),
                    rocplots1,
                    dpi = dpi,
                    width = plot.width,
                    height = plot.height)
  }

  return(rocplots1)

}
