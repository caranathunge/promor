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
#' @param raw_df A \code{raw_df} object (output from \code{\link{create_df}}).
#' @param protein_range The range or subset of proteins (rows) to plot. If not
#' provided, all the proteins (rows) in the data frame will be used.
#' @param sample_range The range of samples to plot. If not
#' provided, all the samples (columns) in the data frame will be used.
#' @param reorder_x Logical. If \code{TRUE} samples on the x axis are reordered
#' using the function given in \code{x_fun}. Default is \code{FALSE}.
#' @param reorder_y Logical. If \code{TRUE} proteins in the y axis are reordered
#' using the function given in \code{y_fun}. Default is \code{FALSE}.
#' @param x_fun Function to reorder samples along the x axis. Possible options
#' are \code{mean} and \code{sum}. Default is \code{mean}.
#' @param y_fun Function to reorder proteins along the y axis. Possible options
#' are \code{mean} and \code{sum}. Default is \code{mean}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis:viridis]{viridis}}
#' for available options.
#' @param label_proteins If \code{TRUE} proteins on the y axis
#' will be labeled with their Majority Protein IDs. Default is \code{FALSE}.
#' @param text_size Text size for axis labels. Default is \code{10}.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' directory provided in \code{file_path}.
#' @param file_path A string containing the directory path to save the file.
#' @param file_name File name to save the heatmap. Default is
#' \code{"Missing_data_heatmap"}.
#' @param file_type File type to save the heatmap. Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{15}.
#' @param plot_height Height of the plot. Default is \code{15}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details
#' This function visualizes patterns of missing value occurrence using a
#' heatmap. The user can choose to reorder the axes using the available functions
#' (\code{x_fun}, \code{y_fun}) to better understand the underlying cause of
#' missing data.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @seealso \code{\link{create_df}}
#'
#' @examples
#' ## Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
#' )
#'
#' ## Missing data heatmap with default settings.
#' heatmap_na(raw_df)
#'
#' ## Missing data heatmap with x and y axes reordered by the mean (default) of
#' ## protein intensity.
#' heatmap_na(raw_df,
#'   reorder_x = TRUE, reorder_y = TRUE
#' )
#'
#' ## Missing data heatmap with x and y axes reordered by the sum of
#' ## protein intensity.
#' heatmap_na(raw_df,
#'   reorder_x = TRUE, reorder_y = TRUE, x_fun = sum,
#'   y_fun = sum
#' )
#'
#' ## Missing data heatmap for a subset of the proteins with x and y axes
#' ## reordered by the mean (default) of protein intensity and the y axis
#' ## labeled with protein IDs.
#' heatmap_na(raw_df,
#'   protein_range = 1:30,
#'   reorder_x = TRUE, reorder_y = TRUE,
#'   label_proteins = TRUE
#' )
#'
#' @export
heatmap_na <- function(raw_df,
                       protein_range,
                       sample_range,
                       reorder_x = FALSE,
                       reorder_y = FALSE,
                       x_fun = mean,
                       y_fun = mean,
                       palette = "viridis",
                       label_proteins = FALSE,
                       text_size = 10,
                       save = FALSE,
                       file_type = "pdf",
                       file_path = NULL,
                       file_name = "Missing_data_heatmap",
                       plot_width = 15,
                       plot_height = 15,
                       dpi = 80) {
  # Binding global variables to the local function
  value <- protgroup <- NULL

  # Assign all rows and columns if protein_range and sample_range aren't defined
  if (missing(protein_range)) {
    protein_range <- 1:nrow(raw_df)
  }
  if (missing(sample_range)) {
    sample_range <- 1:ncol(raw_df)
  }

  raw_df <- as.matrix(raw_df[protein_range, sample_range])


  # Convert the data into long format for plotting and make necessary changes
  # i.e. adding column headers etc.
  hmap_data <- reshape2::melt(raw_df, na.rm = FALSE)
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
  col_na <- set_col(palette, n = 1, direction = -1)
  col_val <- set_col(palette, n = 1, direction = 1)

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
      high = col_val,
      low = col_val,
      na.value = col_na
    ) +
    # coord_equal() +
    ggplot2::theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        size = text_size
      ),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    ggplot2::xlab("") + ggplot2::ylab("")

  # add protein labels
  if (label_proteins == TRUE) {
    hmap <- hmap +
      ggplot2::theme(
        axis.text.y = element_text(size = text_size)
      )
  }

  # Set temporary file_path if not specified
  if (is.null(file_path)) {
    file_path <- tempdir()
  }

  # Save the heatmap as a pdf
  if (save == TRUE) {
    ggplot2::ggsave(paste0(file_path, "/", file_name, ".", file_type),
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
#' containing missing values or a \code{norm_df} object after performing
#' normalization.
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
#' @param seed Numerical. Random number seed. Default is \code{NULL}
#'
#' @details \itemize{\item Ideally, you should first remove proteins with
#' high levels of missing data using the \code{filterbygroup_na} function
#' before running \code{impute_na} on the \code{raw_df} object or the
#' \code{norm_df} object.
#' \item \code{impute_na} function imputes missing values using a
#' user-specified imputation method from the available options, \code{minProb},
#' \code{minDet}, \code{kNN}, \code{RF}, and \code{SVD}.
#' \item **Note: Some imputation methods may require that the data be normalized
#' prior to imputation.**
#' \item Make sure to fix the random number seed with \code{seed} for reproducibility}.
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
#' @return An \code{imp_df} object, which is a data frame of protein intensities
#' with no missing values.
#'
#' @examples
#' ## Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
#' )
#'
#' ## Impute missing values in the data frame using the default minProb
#' ## method.
#' imp_df1 <- impute_na(raw_df, seed = 3312)
#'
#' \donttest{
#' ## Impute using the RF method with the number of iterations set at 5
#' ## and number of trees set at 100.
#' imp_df2 <- impute_na(raw_df,
#'   method = "RF",
#'   maxiter = 5, ntree = 100,
#'   seed = 3312
#' )
#'
#'
#' ## Using the kNN method.
#' imp_df3 <- impute_na(raw_df, method = "kNN", seed = 3312)
#' }
#'
#'
#' ## Using the SVD method with n_pcs set to 3.
#' imp_df4 <- impute_na(raw_df, method = "SVD", n_pcs = 3, seed = 3312)
#'
#' ## Using the minDet method with q set at 0.001.
#' imp_df5 <- impute_na(raw_df, method = "minDet", q = 0.001, seed = 3312)
#'
#' ## Impute a normalized data set using the kNN method
#' imp_df6 <- impute_na(ecoli_norm_df, method = "kNN")
#'
#' @references Lazar, Cosmin, et al. "Accounting for the multiple natures of
#' missing values in label-free quantitative proteomics data sets to compare
#' imputation strategies." Journal of proteome research 15.4 (2016): 1116-1125.
#'
#' @export
impute_na <- function(df,
                      method = "minProb",
                      tune_sigma = 1,
                      q = 0.01,
                      maxiter = 10,
                      ntree = 20,
                      n_pcs = 2,
                      seed = NULL) {
  # Setting global variables to NULL
  value <- protgroup <- NULL

  # Run the user-specified imputation method

  if (method == "minDet") {
    set.seed(seed)
    df_imputed_mindet <- impute.MinDet(df,
      q = q
    )
    return(df_imputed_mindet)
  } else if (method == "RF") {
    set.seed(seed)
    df_imp_temp <- missForest::missForest(df,
      maxiter = maxiter,
      ntree = ntree,
      verbose = TRUE
    )
    df_imputed_rf <- df_imp_temp$ximp
    return(df_imputed_rf)
  } else if (method == "kNN") {
    set.seed(seed)
    df_imputed_knn <- VIM::kNN(df, imp_var = FALSE)
    rownames(df_imputed_knn) <- rownames(df)
    return(df_imputed_knn)
  } else if (method == "SVD") {
    set.seed(seed)
    df <- as.matrix(df)
    df[is.nan(df)] <- NA
    df_imp_temp <- pcaMethods::pca(
      object = df,
      method = "svdImpute",
      n_pcs = n_pcs,
      verbose = TRUE
    )
    df_imputed_svd <- as.data.frame(completeObs(df_imp_temp))
    return(df_imputed_svd)
  } else if (method == "minProb") {
    set.seed(seed)
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
#' containing missing values or a \code{norm_df} object containing normalized
#' protein intensity data.
#' @param imputed An \code{imp_df} object obtained from running \code{impute_na}
#'  on the same data frame provided as \code{original}.
#' @param global Logical. If \code{TRUE} ({default}), a global density plot is
#' produced. If \code{FALSE}, sample-wise density plots are produced.
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis:viridis]{viridis}}
#' for available options.
#' @param n_row Used if \code{global = FALSE} to indicate the number of rows
#' to print the plots.
#' @param n_col Used if \code{global = FALSE} to indicate the number of
#' columns to print the plots.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' directory provided in \code{file_path}.
#' @param file_path A string containing the directory path to save the file.
#' @param file_name File name to save the density plot/s.
#' Default is \code{"Impute_plot."}
#' @param file_type File type to save the density plot/s.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{7}.
#' @param plot_height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details
#' \itemize{\item Given two data frames, one with missing values
#' and the other, an imputed data frame (\code{imp_df} object) of the same
#' data set, \code{impute_plot} generates global or sample-wise density plots
#' to visualize the impact of imputation on the data set.
#' \item Note, when sample-wise option is selected (\code{global = FALSE}),
#' \code{n_col} and \code{n_row} can be used to specify the number of columns
#' and rows to print the plots.
#' \item If you choose to specify \code{n_row} and \code{n_col}, make sure that
#' \code{n_row} * \code{n_col} matches the total number of samples in the
#' data frame.}
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#'
#' ## Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
#' )
#'
#' ## Impute missing values in the data frame using the default minProb
#' ## method.
#' imp_df <- impute_na(raw_df)
#'
#' ## Visualize the impact of missing data imputation with a global density
#' ## plot.
#' impute_plot(original = raw_df, imputed = imp_df)
#'
#' ## Make sample-wise density plots
#' impute_plot(raw_df, imp_df, global = FALSE)
#'
#' ## Print plots in user-specified numbers of rows and columns
#' impute_plot(raw_df, imp_df, global = FALSE, n_col = 2, n_row = 3)
#'
#' @export
impute_plot <- function(original,
                        imputed,
                        global = TRUE,
                        text_size = 10,
                        palette = "viridis",
                        n_row,
                        n_col,
                        save = FALSE,
                        file_path = NULL,
                        file_name = "Impute_plot",
                        file_type = "pdf",
                        plot_width = 7,
                        plot_height = 7,
                        dpi = 80) {
  # Set global variables to null
  value <- NULL

  # Assign n_row and n_col if not defined
  if (missing(n_row) & missing(n_col)) {
    n_row <- ncol(imputed)
    n_col <- 1
  }

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
      ggplot2::geom_density(
        ggplot2::aes(fill = factor(stage,
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
      promor_theme() +
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
      ggplot2::geom_density(
        ggplot2::aes(fill = factor(stage,
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
      promor_facet_theme() +
      ggplot2::theme(
        axis.text = element_text(size = text_size * 0.7),
        legend.position = "bottom",
        legend.text = element_text(size = text_size),
      ) +
      ggplot2::facet_wrap(~sample,
        nrow = n_row,
        ncol = n_col,
        strip.position = "top"
      )

    # Set temporary file_path if not specified
    if (is.null(file_path)) {
      file_path <- tempdir()
    }

    if (save == TRUE) {
      ggplot2::ggsave(paste0(file_path, "/", file_name, ".", file_type),
        s_densplot,
        dpi = dpi,
        width = plot_width * n_col,
        height = plot_height * n_row,
        units = "cm"
      )
      return(s_densplot)
    } else {
      return(s_densplot)
    }
  }
}
