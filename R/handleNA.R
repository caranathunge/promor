# Visualize missing data -------------------------------------
#' Visualize missing data
#' @author Chathurani Ranathunge
#'
#' @description This function visualizes the patterns of missing value
#' occurrence using a heatmap.
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom stats reorder
#' @import viridis
#'
#' @param df A \code{raw_df} object (output from \code{\link{create_df}}).
#' @param protein_range The range of proteins to plot. Default is \code{ALL},
#' meaning all the proteins in the data frame.
#' @param sample_range The range of samples to plot. Default is \code{ALL},
#' meaning all the samples in the data frame.
#' @param reorder_x Logical. If \code{TRUE} samples on the x axis are reordered
#' using the function given in \code{x_fun}. Default is \code{FALSE}.
#' @param reorder_y Logical. If \code{TRUE} proteins in the y axis are reordered
#' using the function given in \code{y_fun}. Default is \code{FALSE}.
#' @param x_fun Function to reorder samples along the x axis. Possible options
#' include \code{mean} and \code{sum}. Default is \code{mean}.
#' @param y_fun Function to reorder proteins along the y axis. Possible options
#' include \code{mean} and \code{sum}. Default is \code{mean}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis: scale_color_viridis]{scale_color_viridis}}
#' for available options.
#' @param text_size Text size for axis labels. Default is \code{10}.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_type File type to save the heatmap. Default is \code{"pdf"}.
#' @param file_name File name to save the heatmap.
#' Default is \code{"MissingData_heatmap"}.
#' @param plot_width Width of the plot. Default is \code{15}.
#' @param plot_height Height of the plot. Default is \code{15}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details
#' This function visualizes patterns of missing value occurrence using a
#' heatmap. User can choose to reorder the axes using the available functions
#' (\code{x_fun}, \code{y_fun}) to better understand the underlying cause of
#' missing data.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @seealso \code{\link{create_df}}
#'
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(prot.groups = "./proteinGroups.txt")
#'
#' ## Missing data heatmap with default settings.
#' heatmap_na(raw)
#'
#' ## Missing data heatmap with x and y axes reordered by sum of intensity.
#' heatmap_na(raw,
#'   reorder_x = TRUE, reorder_y = TRUE, x_fun = sum,
#'   y_fun = sum
#' )
#' }
#'
#' @export
heatmap_na <- function(df,
                       protein_range = "ALL",
                       sample_range = "ALL",
                       reorder_x = FALSE,
                       reorder_y = FALSE,
                       x_fun = mean,
                       y_fun = mean,
                       palette = "viridis",
                       text_size = 10,
                       save = FALSE,
                       file_type = "pdf",
                       file_name = "MissingData_heatmap",
                       plot_width = 15,
                       plot_height = 15,
                       dpi = 80) {

  # Binding global variables to the local function
  value <- protgroup <- NULL

  # select the range of data to plot
  if (protein_range == "ALL" && sample_range != "ALL") {
    df <- as.matrix(df[, sample_range])
  } else if (protein_range != "ALL" && sample_range == "ALL") {
    df <- as.matrix(df[protein_range, ])
  } else if (protein_range != "ALL" && sample_range != "ALL") {
    df <- as.matrix(df[protein_range, sample_range])
  } else {
    df <- as.matrix(df)
  }

  # Convert the data into long format for plotting and make necessary changes
  # i.e. adding column headers etc.
  hmap_data <- reshape2::melt(df, na.rm = FALSE)
  hmap_data$mjprot <- sapply(
    strsplit(as.character(hmap_data[, 1]), ";"),
    getElement, 1
  )
  hmap_data[1] <- NULL
  hmap_df <- as.data.frame(hmap_data)
  colnames(hmap_df) <- c(
    "sample",
    "value",
    "protgroup"
  )
  # Set colors
  col.na <- set_col(palette, n = 1, direction = -1)
  col.val <- set_col(palette, n = 1, direction = 1)

  # Options for arranging the rows and the columns of the heat map
  if (reorder_x == TRUE & reorder_y == TRUE) {
    hmap <- ggplot2::ggplot(
      hmap_df,
      ggplot2::aes(
        x = stats::reorder(sample,
          value,
          na.rm = TRUE,
          FUN = x_fun
        ),
        y = stats::reorder(protgroup,
          value,
          na.rm = TRUE,
          FUN = y_fun
        ),
        fill = value
      )
    )
  } else if (reorder_x == FALSE & reorder_y == TRUE) {
    hmap <- ggplot2::ggplot(
      hmap_df,
      ggplot2::aes(
        x = sample,
        y = stats::reorder(protgroup,
          value,
          na.rm = TRUE,
          FUN = y_fun
        ),
        fill = value
      )
    )
  } else if (reorder_x == TRUE & reorder_y == FALSE) {
    hmap <- ggplot2::ggplot(
      hmap_df,
      ggplot2::aes(
        x = stats::reorder(sample,
          value,
          na.rm = TRUE,
          FUN = x_fun
        ),
        y = protgroup,
        fill = value
      )
    )
  } else {
    hmap <- ggplot2::ggplot(
      hmap_df,
      ggplot2::aes(
        x = sample,
        y = protgroup,
        fill = value
      )
    )
  }

  # Create heat map
  hmap <- hmap +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(
      high = col.val,
      low = col.val,
      na.value = col.na
    ) +
    coord_equal() +
    ggplot2::theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        size = text_size
      ),
      axis.text.y = element_text(size = text_size),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    ggplot2::xlab("") + ggplot2::ylab("")

  # Save the heatmap as a pdf
  if (save == TRUE) {
    ggplot2::ggsave(paste0(file_name, ".", file_type),
      hmap,
      dpi = dpi,
      width = plot_width,
      height = plot_height
    )
    return(hmap)
  } else {
    return(hmap)
  }
}

# Impute missing data -----------------------------------------------------
#' Impute missing values
#' @author Chathurani Ranathunge
#' @description This function imputes missing values using a user-specified
#' imputation method.
#'
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca
#' @importFrom pcaMethods completeObs
#' @importFrom VIM kNN
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom stats quantile
#' @importFrom stats median
#'
#'
#' @param df A \code{raw_df} object (output of \code{\link{create_df}})
#' containing missing values.
#' @param method Imputation method to use. Default is \code{"minProb"}.
#' Available methods: \code{"minDet", "RF", "kNN", and "SVD"}.
#' @param tune_sigma A scalar used in the \code{"minProb"} method for
#' controlling the standard deviation of the Gaussian distribution
#' from which random values are drawn for imputation.\cr
#' Default is 1.
#' @param q A scalar used in \code{"minProb"} and \code{"minDet"} methods
#' to obtain a low intensity value for imputation. \code{q} should be set to a
#' very low value. Default is 0.01.
#' @param maxiter Maximum number of iterations to be performed when using the
#' \code{"RF"} method. Default is \code{10}.
#' @param ntree Number of trees to grow in each forest when using the
#' \code{"RF"} method. Default is \code{20}.
#' @param n_pcs Number of principal components to calculate when using the
#' \code{"SVD"} method. Default is 2.
#'
#' @details \code{impute_na} function imputes missing values using a
#' user-specified imputation method from the available options, \code{minProb},
#' \code{minDet}, \code{kNN}, \code{RF}, and \code{SVD}.
#'
#' @seealso More information on the available imputation methods can be found
#' in their respective packages.
#' \itemize{\item \code{\link{create_df}}
#' \item For \code{minProb} and
#' \code{minDet} methods, see
#' \code{imputeLCMD} package.
#' \item For Random Forest (\code{RF}) method, see
#'  \code{\link[missForest]{missForest}}.
#' \item For \code{kNN} method, see \code{\link[VIM]{kNN}} from the
#'  \code{\link[VIM]{VIM}} package.
#' \item For \code{SVD} method, see \code{\link[pcaMethods]{pca}} from the
#' \code{\link[pcaMethods]{pcaMethods}} package.}
#'
#' @return An \code{imp.df} object, which is a data frame of intensities with no
#' missing values.
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(prot_groups = "./proteinGroups.txt")
#'
#' ## Impute missing values in the data frame using the default minProb
#' ## method.
#' imp_raw <- impute_na(raw)
#'
#' ## Using the RF method with the number of iterations set at 5
#' ## and number of trees set at 100.
#' imp_raw <- impute_na(raw, method = "RF", maxiter = 5, ntree = 100)
#'
#' ## Using the kNN method.
#' imp_raw <- impute_na(raw, method = "kNN")
#'
#' ## Using the SVD method with n_pcs set to 3.
#' imp_raw <- impute_na(raw, method = "SVD", n_pcs = 3)
#'
#' ## Using the minDet method with q set at 0.001.
#' imp_raw <- impute_na(raw, method = "minDet", q = 0.001)
#' }
#'
#' @export
impute_na <- function(df,
                      method = "minProb",
                      tune_sigma = 1,
                      q = 0.01,
                      maxiter = 10,
                      ntree = 20,
                      n_pcs = 2) {

  # Setting global variables to NULL
  value <- protgroup <- NULL

  # Run the user-specified imputation method

  if (method == "minDet") {
    df_imputed_mindet <- impute.MinDet(df,
                                       q = q)
    return(df_imputed_mindet)
  } else if (method == "RF") {
    df_imp_temp <- missForest::missForest(df,
      maxiter = maxiter,
      ntree = ntree,
      verbose = TRUE
    )
    df_imputed_rf <- df_imp_temp$ximp
    return(df_imputed_rf)
  } else if (method == "kNN") {
    df_imputed_knn <- VIM::kNN(df, imp_var = FALSE)
    rownames(df_imputed_knn) <- rownames(df)
    return(df_imputed_knn)
  } else if (method == "SVD") {
    df <- as.matrix(df)
    df[is.nan(df)] <- NA
    df_imp_temp <- pcaMethods::pca(
      object = df,
      method = "svdImpute",
      n_pcs = n_pcs,
      verbose = TRUE
    )
    df_imputed_svd <- completeObs(df_imp_temp)
    return(df_imputed_svd)
  } else if (method == "minProb") {
    df_imputed_minprob <- impute.Min.Prob(df,
      q = q,
      tune_sigma = tune_sigma
    )
    return(df_imputed_minprob)
  }
}


# Visualize the imputation effects ----------------------------------------
#' Visualize the impact of imputation
#' @author Chathurani Ranathunge
#' @description This function generates density plots to visualize the impact of
#' missing data imputation on the data.
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import viridis
#'
#' @param original A \code{raw_df} object (output of \code{\link{create_df}})
#' containing missing values.
#' @param imputed An \code{imp.df} object obtained from running \code{impute_na}
#'  on the same data frame provided as \code{original}.
#' @param global Logical. If \code{TRUE} ({default}), a global density plot is
#' produced. If \code{FALSE}, sample-wise density plots are produced.
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis: scale_color_viridis]{scale_color_viridis}}
#' for available options.
#' @param nrow Required if \code{global = FALSE} to indicate the number of rows
#' to print the plots.
#' @param ncol Required if \code{global = FALSE} to indicate the number of
#' columns to print the plots.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_name file_name File name to save the density plot/s.
#' Default is \code{"Impute_plot."}
#' @param file_type File type to save the density plot/s.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{7}.
#' @param plot_height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details
#' \itemize{\item Given two data frames, one with missing values
#' (\code{raw_df} object) and the other, an imputed data frame
#' (\code{imp.df} object) of the same data set, \code{impute_plot}
#' generates global or sample-wise density plots to visualize the
#' impact of imputation on the data set.
#' \item Note, when sample-wise option is selected (\code{global = FALSE}),
#' \code{nrow} * \code{ncol} should match the number of samples in the
#' data frame.}
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(prot.groups = "./proteinGroups.txt")
#'
#' ## Impute missing values in the data frame.
#' imp_raw <- impute_na(raw)
#'
#' ## Visualize the impact of missing data imputation with a global density
#' ## plot.
#' impute_plot(original = raw, imputed = imp_raw)
#'
#' ## Make sample-wise density plots for a data set of 25 samples.
#' impute_plot(raw, imp_raw, global = FALSE, nrow = 5, ncol = 5)
#' }
#' @export
impute_plot <- function(original,
                        imputed,
                        global = TRUE,
                        text_size = 10,
                        palette = "viridis",
                        nrow,
                        ncol,
                        save = FALSE,
                        file_name = "Impute_plot",
                        file_type = "pdf",
                        plot_width = 7,
                        plot_height = 7,
                        dpi = 80) {

  # Set global variables to null
  value <- NULL
  # Make necessary changes and combine data frames before plotting
  original <- as.matrix(original)
  orig_data <- reshape2::melt(original, na.rm = FALSE)
  orig_data$mjprot <- sapply(
    strsplit(as.character(orig_data[, 1]), ";"),
    getElement, 1
  )
  orig_data[1] <- NULL
  colnames(orig_data) <- c("sample", "value", "protgroup")
  orig_data$stage <- "Before imputation"

  # do the same for imputed data set
  imputed <- as.matrix(imputed)
  imp_data <- reshape2::melt(imputed, na.rm = FALSE)
  imp_data$mjprot <- sapply(
    strsplit(as.character(imp_data[, 1]), ";"),
    getElement, 1
  )
  imp_data[1] <- NULL
  colnames(imp_data) <- c("sample", "value", "protgroup")
  imp_data$stage <- "After imputation"
  dplot_data <- rbind(orig_data, imp_data)

  # Plot global density plot

  if (global == TRUE) {
    g_densityplot <- ggplot2::ggplot(
      dplot_data,
      ggplot2::aes(x = value)
    ) +
      ggplot2::geom_density(ggplot2::aes(fill = factor(stage,
        levels = c(
          "Before imputation",
          "After imputation"
        )
      )),
      alpha = 0.7,
      lwd = text_size * 0.02
      ) +
      ggplot2::xlab("Intensity") +
      ggplot2::ylab("Density") +
      viridis::scale_fill_viridis(
        discrete = TRUE,
        option = palette,
        begin = 0.3,
        end = 0.7
      ) +
      promor_theme +
      ggplot2::theme(
        axis.title.x = element_text(
          size = text_size,
          face = "bold"
        ),
        axis.title.y = element_text(
          size = text_size,
          face = "bold"
        ),
        axis.text = element_text(size = text_size * 0.7),
        legend.position = "bottom",
        legend.text = element_text(size = text_size)
      )


    if (save == TRUE) {
      ggplot2::ggsave(paste0(file_name, ".", file_type),
        g_densityplot,
        dpi = dpi,
        width = plot_width,
        height = plot_height
      )
      return(g_densityplot)
    } else {
      return(g_densityplot)
    }
  } else {
    s_densplot <- ggplot2::ggplot(
      dplot_data,
      ggplot2::aes(x = value)
    ) +
      ggplot2::geom_density(ggplot2::aes(fill = factor(stage,
        levels = c(
          "Before imputation",
          "After imputation"
        )
      )),
      alpha = 0.7,
      lwd = text_size * 0.02
      ) +
      ggplot2::xlab("Intensity") +
      ggplot2::ylab("Density") +
      viridis::scale_fill_viridis(
        discrete = TRUE,
        option = palette,
        begin = 0.3,
        end = 0.7
      ) +
      promor_facet_theme +
      ggplot2::theme(
        axis.text = element_text(size = text_size * 0.7),
        legend.position = "bottom",
        legend.text = element_text(size = text_size),
      ) +
      ggplot2::facet_wrap(~sample,
        nrow = nrow,
        ncol = ncol,
        strip.position = "top"
      )
    if (save == TRUE) {
      ggplot2::ggsave(paste0(file_name, ".", file_type),
        s_densplot,
        dpi = dpi,
        width = plot_width * ncol,
        height = plot_height * nrow,
        units = "cm"
      )
      return(s_densplot)
    } else {
      return(s_densplot)
    }
  }
}
