# Feature plots ----------------------------------------------------------------
#' Visualize feature (protein) variation among conditions
#' @description This function visualizes protein intensity differences among
#' conditions (classes) using box plots or density distribution plots.
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import viridis
#'
#'
#' @param model_df A \code{model_df} object from performing \code{pre_process}.
#' @param type Type of plot to generate. Choices are "box" or "density." Default
#' is \code{"box."}
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param n_row Number of rows to print the plots.
#' @param n_col Number of columns to print the plots.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_name file_name File name to save the plot.
#' Default is \code{"Feature_plot."}
#' @param file_type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{7}.
#' @param plot_height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details This function visualizes condition-wise differences in protein
#' intensity using boxplots and/or density plots.
#'
#' @return A \code{ggplot2} object
#' @seealso
#' \itemize{
#' \item \code{pre_process}, \code{rem_feature}}
#' @examples
#' \dontrun{
#'
#' ## Create a model_df object with default settings.
#' covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#'
#' ## Feature variation - box plots
#' feature_plot(covid_model_df, type = "box", n_row = 4, n_col = 2)
#'
#' ## Density plots
#' feature_plot(covid_model_df, type = "density")
#'
#' ## Change color palette
#' feature_plot(covid_model_df, type = "density", n_row = 4, n_col = 2, palette = "rocket")
#' }
#' @export
feature_plot <- function(model_df,
                         type = "box",
                         text_size = 10,
                         palette = "viridis",
                         n_row,
                         n_col,
                         save = FALSE,
                         file_name = "Feature_plot",
                         file_type = "pdf",
                         dpi = 80,
                         plot_width = 7,
                         plot_height = 7) {

  #binding for global variable
  Intensity <- Condition <- NULL

  # Create plot data
  modeldf_melted <- reshape2::melt(model_df)
  colnames(modeldf_melted) <- c("Condition", "Protein", "Intensity")

  # n_row and n_col if not defined
  if (missing(n_row) & missing(n_col)) {
    n_row <- length(unique(modeldf_melted$Protein))
    n_col <- 1
  }


  # Make density plots
  if (type == "density") {
    pred_plot <- ggplot2::ggplot(modeldf_melted, aes(
      x = Intensity,
      color = Condition
    )) +
      ggplot2::geom_density(
        alpha = 0.25,
        lwd = text_size * 0.07
      ) +
      viridis::scale_color_viridis(
        discrete = TRUE,
        option = palette,
        begin = 0.3,
        end = 0.7
      ) +
      ggplot2::facet_wrap(~Protein,
        scale = "free",
        nrow = n_row,
        ncol = n_col
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      promor_facet_theme() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_text(size = text_size),
        axis.title.y = element_text(size = text_size),
        axis.text = element_text(size = text_size * 0.5),
        legend.text = element_text(size = text_size)
      )
  } else {
    # Make box plots
    pred_plot <- ggplot2::ggplot(
      modeldf_melted,
      ggplot2::aes(
        x = Condition,
        y = Intensity,
        col = Condition
      )
    ) +
      geom_boxplot(
        alpha = 0.7,
        outlier.shape = 1,
        outlier.stroke = 0.1,
        outlier.size = text_size * 0.1,
        outlier.color = "grey50",
        lwd = text_size * 0.07
      ) +
      viridis::scale_color_viridis(
        discrete = TRUE,
        option = palette,
        begin = 0.3,
        end = 0.7
      ) +
      ggplot2::facet_wrap(~Protein,
        scales = "free",
        nrow = n_row,
        ncol = n_col
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      promor_facet_theme() +
      ggplot2::theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = text_size),
        axis.title.y = element_text(size = text_size)
      )
  }
  if (save == TRUE) {
    ggplot2::ggsave(paste0(file_name, ".", file_type),
      pred_plot,
      dpi = dpi,
      width = plot_width,
      height = plot_height
    )
  }

  return(pred_plot)
}
# Variable importance plot-----------------------------------------------------
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
#' @param model_list A \code{model_list} object from performing
#' \code{train_models}.
#' @param ... Additional arguments to be passed on to
#' \code{\link[caret:varImp]{varImp}}.
#' @param type Type of plot to generate. Choices are "bar" or "lollipop."
#' Default is \code{"lollipop."}
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param n_row Number of rows to print the plots.
#' @param n_col Number of columns to print the plots.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_name file_name File name to save the plot.
#' Default is \code{"VarImp_plot."}
#' @param file_type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{7}.
#' @param plot_height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details \itemize{\item \code{varimp_plot} produces a list of plots showing
#' variable importance measures calculated from models generated with different
#' machine-learning algorithms.
#' \item Note: Variables are ordered by variable importance in
#' descending order, and by default, importance values are scaled to 0 and 100.
#' This can be changed by specifying \code{scale = FALSE}. See
#' \code{\link[caret:varImp]{varImp}} for more information.}
#'
#' @return A list of \code{ggplot2} objects.
#' @seealso
#' \itemize{
#' \item \code{train_models}, \code{rem_feature}
#' \item \code{\link[caret:varImp]{caret: varImp}}
#' }
#'
#' @examples
#' \dontrun{
#'
#' ## Create a model_df object
#' covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#'
#' ## Split the data frame into training and test data sets
#' covid_split_df <- split_data(covid_model_df)
#'
#' ## Fit models based on the default list of machine learning (ML) algorithms
#' covid_model_list <- train_models(covid_split_df)
#'
#' ## Variable importance - lollipop plots
#' varimp_plot(covid_model_list)
#'
#' ## Bar plots
#' varimp_plot(covid_model_list, type = "bar")
#'
#' ## Do not scale variable importance values
#' varimp_plot(covid_model_list, scale = FALSE)
#'
#' ## Change color palette
#' varimp_plot(covid_model_list, palette = "magma")
#' }
#' @export
varimp_plot <- function(model_list,
                        ...,
                        type = "lollipop",
                        text_size = 10,
                        palette = "viridis",
                        n_row,
                        n_col,
                        save = FALSE,
                        file_name = "VarImp_plot",
                        file_type = "pdf",
                        dpi = 80,
                        plot_width = 7,
                        plot_height = 7) {

  #binding for global variable
  Protein <- Importance <- NULL

  # Assign n_row and n_col
  if (missing(n_row) & missing(n_col)) {
    n_row <- length(model_list)
    n_col <- 1
  }

  # Calculate variable importance with VarImp for each ML algorithm-based
  # model in the list
  vimp <- lapply(
    model_list,
    function(x) {
      tryCatch(caret::varImp(
        x,
        ...
      ),
      error = function(e) NULL
      )
    }
  )


  # Drop Null items from the list
  vimp <- Filter(Negate(is.null), vimp)

  # Make necessary changes to 'importance data frames' in the list before
  # plotting.
  plot_imp_data <- lapply(
    vimp,
    function(x) {
      x <- x$importance[1]
      x$protein <- rownames(x)
      rownames(x) <- NULL
      colnames(x) <- c("Importance", "Protein")
      x$Protein <- gsub("`|\\\\", "", x$Protein)
      x$Importance[is.nan(x$Importance)] <- NA
      x <- x[rowSums(is.na(x)) == 0, ]
      x
    }
  )

  # Drop empty data frames before proceeding
  plot_imp_data <- Filter(function(x) dim(x)[1] > 0, plot_imp_data)

  # Make plots based on type

  if (type == "bar") {
    # Make bar plots
    vi_plots <- lapply(plot_imp_data, function(t) {
      ggplot2::ggplot(
        data = t,
        aes(
          x = reorder(Protein, Importance),
          y = Importance,
          fill = Importance
        )
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_flip() +
        viridis::scale_fill_viridis(
          option = palette,
          direction = -1
        ) +
        xlab("") +
        promor_facet_theme() +
        ggplot2::theme(
          plot.title = element_text(
            size = text_size,
            face = "bold"
          ),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.margin = margin(0, 0, 0, 0, unit = "cm"),
          legend.key.width = grid::unit(0.8, "cm"),
          legend.key.height = grid::unit(0.2, "cm"),
          legend.title = element_text(size = text_size * 0.7),
          legend.text = element_text(size = text_size * 0.5),
          axis.text.y = element_text(
            size = text_size * 0.8,
            face = "bold"
          ),
          panel.grid.major.x = element_line(
            size = 0.1,
            color = "grey80"
          )
        ) +
        ggplot2::guides(fill = guide_colorbar(
          title.position = "top",
        ))
    })

    # Add plot titles
    vi_plots <- lapply(seq_along(vi_plots), function(i) {
      vi_plots[[i]] + ggtitle(gsub(".importance", "", names(vi_plots)[i]))
    })
  } else {
    # Make lollipop plots
    vi_plots <- lapply(plot_imp_data, function(t) {
      ggplot2::ggplot(
        data = t,
        aes(
          x = reorder(Protein, Importance),
          y = Importance,
          color = Importance
        )
      ) +
        ggplot2::geom_segment(aes(
          xend = Protein,
          y = 0,
          yend = Importance
        ),
        lwd = text_size * 0.2
        ) +
        ggplot2::geom_point(aes(color = Importance),
          size = text_size,
          show.legend = FALSE
        ) +
        ggplot2::coord_flip() +
        ggplot2::geom_label(
          label = round(t$Importance,
            digits = 1
          ),
          fill = NA,
          label.size = NA,
          colour = "white",
          size = text_size * 0.3
        ) +
        ggplot2::xlab("") +
        viridis::scale_color_viridis(
          option = palette,
          direction = -1
        ) +
        promor_facet_theme() +
        ggplot2::theme(
          plot.title = element_text(
            size = text_size,
            face = "bold"
          ),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.margin = margin(0, 0, 0, 0, unit = "cm"),
          legend.key.width = grid::unit(0.8, "cm"),
          legend.key.height = grid::unit(0.2, "cm"),
          legend.title = element_text(size = text_size * 0.7),
          legend.text = element_text(size = text_size * 0.5),
          axis.text.y = element_text(
            size = text_size * 0.8,
            face = "bold"
          ),
          axis.text.x = element_text(
            size = text_size * 0.9
          ),
          panel.grid.major.x = element_line(
            size = 0.1,
            color = "grey80"
          )
        ) +
        ggplot2::guides(colour = guide_colourbar(
          title = "Importance",
          title.position = "top",
          title.hjust = 0.5
        ))
    })

    # Add plot titles
    vi_plots <- lapply(seq_along(vi_plots), function(i) {
      vi_plots[[i]] + ggtitle(gsub(".importance", "", names(vi_plots)[i]))
    })
  }

  if (save == TRUE) {
    ggplot2::ggsave(paste0(file_name, ".", file_type),
      marrangeGrob(
        grobs = vi_plots,
        nrow = n_row,
        ncol = n_col,
        top = ""
      ),
      dpi = dpi,
      width = plot_width,
      height = plot_height,
    )
    grid.arrange(grobs = vi_plots, newpage = TRUE)
  } else {
    grid.arrange(grobs = vi_plots, newpage = TRUE)
  }
}
# Model performance plot-------------------------------------------------------
#' Model performance plot
#' @description This function generates plots to visualize model performance
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
#' @param model_list A \code{model_list} object from performing
#' \code{train_models}.
#' @param type Type of plot to generate. Choices are "box" or "dot."
#' Default is \code{"box."} for boxplots.
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_name file_name File name to save the plot.
#' Default is \code{"Performance_plot."}
#' @param file_type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{7}.
#' @param plot_height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details \itemize{\item \code{performance_plot} uses resampling results from
#' models included in the \code{model_list} to generate plots showing model
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
#' \item\code{\link[caret: resamples]{caret: resamples}}
#' \item \code{\link[caret:train]{caret: train}}
#' \item \code{\link[caret:trainControl]{caret: trainControl}}
#' }
#' @examples
#' \dontrun{
#'
#' ## Create a model_df object
#' covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#'
#' ## Split the data frame into training and test data sets
#' covid_split_df <- split_data(covid_model_df)
#'
#' ## Fit models based on the default list of machine learning (ML) algorithms
#' covid_model_list <- train_models(covid_split_df)
#'
#' ## Generate box plots to visualize performance of different ML algorithms
#' performance_plot(covid_model_list)
#'
#' ## Generate dot plots
#' performance_plot(covid_model_list, type = "dot")
#'
#' ## Change color palette
#' performance_plot(covid_model_list, type = "dot", palette = "inferno")
#' }
#' @export
performance_plot <- function(model_list,
                             type = "box",
                             text_size = 10,
                             palette = "viridis",
                             save = FALSE,
                             file_name = "Performance_plot",
                             file_type = "pdf",
                             plot_width = 7,
                             plot_height = 7,
                             dpi = 80) {

  #binding for global variable
  method <- value <- NULL

  # Extract resample data from the model_list object
  resample_data <- resamples(model_list)

  # Convert the performance data into long-form format for plotting.
  plot_data <- reshape2::melt(resample_data$values)

  # Extract method and metric information from the variable column and add them
  # to new columns.
  plot_data$method <- sapply(
    strsplit(
      as.character(plot_data$variable),
      "~"
    ),
    "[", 1
  )
  plot_data$metric <- sapply(
    strsplit(
      as.character(plot_data$variable),
      "~"
    ),
    "[", 2
  )

  # Remove the variable column as it is no longer needed.
  plot_data$variable <- NULL

  # Assign palette color
  pal_col <- set_col(palette, 1)

  # Make plots
  if (type == "dot") {
    # Make dot plots
    perform_plot <- ggplot2::ggplot(
      plot_data,
      aes(
        x = method,
        y = value
      )
    ) +
      ggplot2::stat_summary(
        fun.data = "mean_se",
        lwd = text_size * 0.1,
        color = pal_col
      ) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "pointrange",
        size = text_size * 0.1,
        pch = 16,
        position = "identity",
        color = pal_col
      ) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~metric,
        scales = "free"
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      promor_facet_theme() +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y = element_text(
          size = text_size * 0.8,
          face = "bold"
        )
      )
  } else {
    # Make box plots
    perform_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = method,
        y = value
      )
    ) +
      geom_boxplot(
        color = pal_col,
        width = text_size * 0.03,
        alpha = 0.7,
        outlier.shape = 1,
        outlier.stroke = 0.1,
        outlier.size = text_size * 0.1,
        outlier.color = "grey30",
        lwd = text_size * 0.1
      ) +
      ggplot2::facet_wrap(~metric,
        scales = "free"
      ) +
      ggplot2::coord_flip() +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      promor_facet_theme() +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y = element_text(
          size = text_size * 0.8,
          face = "bold"
        )
      )
  }

  if (save == TRUE) {
    ggplot2::ggsave(paste0(file_name, ".", file_type),
      perform_plot,
      dpi = dpi,
      width = plot_width,
      height = plot_height
    )
  }
  return(perform_plot)
}


# ROC plot--------------------------------------------------------------------
#' ROC plot
#' @description This function generates Receiver Operating Characteristic (ROC)
#' curves to evaluate models
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
#' @param probability_list A \code{probability_list} object from performing
#' \code{test_models} with \code{type = "prob"}.
#' @param split_df A \code{split_df} object from performing \code{split_data}
#' @param ... Additional arguments to be passed on to
#' \code{\link[pROC: roc]{roc}}.
#' @param multiple_plots Logical. If \code{FALSE} plots all ROC curves
#' representing algorithms included in the \code{probability_list} in a single
#' plot.
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_name file_name File name to save the plot.
#' Default is \code{"ROC_plot."}
#' @param file_type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{7}.
#' @param plot_height Height of the plot. Default is \code{7}.
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
#' ## Create a model_df object
#' covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#'
#' ## Split the data frame into training and test data sets
#' covid_split_df <- split_data(covid_model_df)
#'
#' ## Fit models using the default list of machine learning (ML) algorithms
#' covid_model_list <- train_models(covid_split_df)
#'
#' # Test a list of models on a test data set and output class probabilities,
#' covid_prob_list <- test_models(covid_model_list, covid_split_df, type = "prob")
#'
#' ## Plot ROC curves separately for each ML algorithm
#' roc_plot(covid_prob_list, covid_split_df)
#'
#' ## Plot all ROC curves in one plot
#' roc_plot(covid_prob_list, covid_split_df, multiple_plots = FALSE)
#'
#' ## Change color palette
#' roc_plot(covid_prob_list, covid_split_df, palette = "plasma")
#' }
#' @export
roc_plot <- function(probability_list,
                     split_df,
                     ...,
                     multiple_plots = TRUE,
                     text_size = 10,
                     palette = "viridis",
                     save = FALSE,
                     file_name = "ROC_plot",
                     file_type = "pdf",
                     plot_width = 7,
                     plot_height = 7,
                     dpi = 80) {

  #binding for global variable
  method <- auc <- NULL

  # Generate ROC objects and add them to a list
  roc_list <- lapply(
    probability_list,
    function(x) {
      pROC::roc(x[[1]],
        response = split_df$test$condition,
        percent = TRUE,
        ...
      )
    }
  )


  # Extract auc values and create a dataframe for plot labels
  auc_list <- lapply(
    roc_list,
    function(x) {
      round(
        x$auc,
        1
      )
    }
  )
  auc_labels <- as.data.frame(do.call(
    "rbind",
    auc_list
  ))
  auc_labels$method <- rownames(auc_labels)
  rownames(auc_labels) <- NULL
  names(auc_labels) <- c(
    "auc",
    "method"
  )

  # Extract sensitivities and specificities from ROC object and make a data
  # frame for plotting
  sens_list <- lapply(
    roc_list,
    function(x) x$sensitivities
  )
  sens <- do.call(
    "rbind",
    sens_list
  )
  sens_melted <- reshape2::melt(sens)
  names(sens_melted) <- c(
    "method",
    "number",
    "sensitivity"
  )

  spec_list <- lapply(
    roc_list,
    function(x) x$specificities
  )
  spec <- do.call(
    "rbind",
    spec_list
  )
  spec_melted <- reshape2::melt(spec)
  names(spec_melted) <- c(
    "method",
    "number",
    "specificity"
  )

  # Merge to create the final data frame for plotting
  roc_plotdata <- merge(
    sens_melted,
    spec_melted
  )

  # Order data frame by sensitivity. Important for geom_step
  roc_plotdata <- roc_plotdata[order(
    roc_plotdata$method,
    roc_plotdata$sensitivity
  ), ]

  # Make ROC curves
  rocplots <- ggplot2::ggplot(roc_plotdata, aes(
    x = specificity,
    y = sensitivity,
    colour = method
  )) +
    ggplot2::geom_step(
      direction = "hv",
      lwd = text_size * 0.1
    ) +
    ggplot2::scale_x_reverse(lim = c(100, 0), ) +
    viridis::scale_color_viridis(
      discrete = TRUE,
      option = palette,
      begin = 0.3,
      end = 0.7
    ) +
    ggplot2::xlab("Specificity (%)") +
    ggplot2::ylab("Sensitivity (%)") +
    ggplot2::geom_abline(
      intercept = 100,
      slope = 1,
      color = "grey60",
      linetype = 2,
      show.legend = FALSE
    )


  if (multiple_plots == FALSE) {
    rocplots1 <- rocplots +
      promor_theme() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = element_text(
          size = text_size * 0.9,
          face = "bold"
        ),
        axis.title = element_text(size = text_size)
      )
  } else {
    rocplots1 <- rocplots +
      ggplot2::facet_wrap(~method) +
      promor_facet_theme() +
      ggplot2::theme(
        text = element_text(size = text_size),
        legend.position = "none",
        axis.title.x = element_text(size = text_size),
        axis.title.y = element_text(
          size = text_size,
          angle = 90
        )
      ) +
      ggplot2::geom_label(
        data = auc_labels,
        aes(label = paste0(
          "AUC: ",
          auc,
          "%"
        )),
        x = Inf,
        y = -Inf,
        hjust = 1,
        vjust = -0.5,
        inherit.aes = FALSE,
        label.size = NA
      )
  }
  if (save == TRUE) {
    ggplot2::ggsave(paste0(
      file_name,
      ".",
      file_type
    ),
    rocplots1,
    dpi = dpi,
    width = plot_width,
    height = plot_height
    )
  }

  return(rocplots1)
}
