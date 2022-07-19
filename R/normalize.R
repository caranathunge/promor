# Normalize intensities ----------------------------------------------------
#' Normalize intensity data
#' @author Chathurani Ranathunge
#' @description This function normalizes data using a user-specified
#' normalization method.
#' @import limma
#'
#' @param imp_df An \code{imp_df} object with missing values imputed using
#' \code{impute_na}.
#' @param method Name of the normalization method to use. Choices are
#' \code{"none", "scale", "quantile" or "cyclicloess."}
#' Default is \code{"quantile."}
#'
#' @details \itemize{\item\code{normalize_data} is a wrapper function around
#' the \code{\link[limma]{normalizeBetweenArrays}} function from the
#' \code{limma} package. \item This function normalizes
#' intensity values to achieve consistency among samples.
#' \item It assumes that the intensities in the
#' data frame have been log-transformed, therefore, it is important to make sure
#' that \code{create_df} was run with \code{log_tr = TRUE}(default) when
#' creating the \code{raw_df} object.}
#'
#' @return A \code{norm_df} object, which is a data frame of
#' normalized protein intensities.
#'
#' @seealso \itemize{\item \code{create_df}
#' \item \code{impute_na}
#' \item See \code{\link[limma]{normalizeBetweenArrays}} in the R package
#' \code{limma} for more information on the different normalization methods
#' available.}
#'
#'
#' @examples
#' ## Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#' prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#' exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
#' )
#'
#' ## Impute missing values in the data frame using the default minProb
#' ## method.
#' imp_df <- impute_na(raw_df)
#'
#' ## Normalize the imp_df object using the default quantile method
#' norm_df1 <- normalize_data(imp_df)
#'
#' ## Use the cyclicloess method
#' norm_df2 <- normalize_data(imp_df, method = "cyclicloess")
#'
#' @export

normalize_data <- function(imp_df,
                           method = "quantile") {
  norm_df <- limma::normalizeBetweenArrays(imp_df,
    method = method
  )

  return(norm_df)
}

# Visualize normalization effects -----------------------------------------
#' Visualize the effect of normalization
#' @author Chathurani Ranathunge
#' @description This function visualizes the impact of normalization on
#' the data
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import viridis


#' @param original A \code{raw_df} object (output of \code{\link{create_df}})
#' containing missing values, or ideally, an \code{imp_df} object after
#' imputing the missing values with \code{impute_na}.
#' @param normalized A \code{norm_df} object after normalizing the data frame
#' provided as \code{original} using \code{normalize_data}.
#' @param type Type of plot to generate. Choices are "box" or "density." Default
#' is \code{"box."}
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' directory provided in \code{file_path}.
#' @param file_path A string containing the directory path to save the file.
#' @param file_name File name to save the plot.
#' Default is \code{"Norm_plot."}
#' @param file_type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{10}.
#' @param plot_height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details Given two data frames, one with data prior to normalization
#' (\code{original}), and the other, after normalization (\code{normalized}),
#' \code{norm_plot} generates side-by-side plots to visualize the effect of
#' normalization on the protein intensity data.
#'
#' @seealso \itemize{\item \code{\link{normalize_data}}
#' \item \code{create_df}
#' \item \code{impute_na}
#' }
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \donttest{
#' ## Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#' prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#' exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
#' )
#'
#' ## Impute missing values in the data frame using the default minProb
#' ## method.
#' imp_df <- impute_na(raw_df)
#'
#' ## Normalize the imp_df object using the default quantile method
#' norm_df <- normalize_data(imp_df)
#'
#' ## Visualize normalization using box plots
#' norm_plot(original = imp_df, normalized = norm_df)
#'
#' ## Visualize normalization using density plots
#' norm_plot(imp_df, norm_df, type = "density")
#' }
#'
#' @export
norm_plot <- function(original,
                      normalized,
                      type = "box",
                      text_size = 10,
                      palette = "viridis",
                      save = FALSE,
                      file_path = NULL,
                      file_name = "Norm_plot",
                      file_type = "pdf",
                      dpi = 80,
                      plot_width = 10,
                      plot_height = 7) {

  # Set global variables to null
  intensity <- value <- group <- NULL

  # Pre-prossesing data for plotting
  normalized <- as.matrix(normalized)
  norm1 <- reshape2::melt(normalized, na.rm = FALSE)
  norm1$normstage <- "After normalization"

  original <- as.matrix(original)
  orig1 <- reshape2::melt(original, na.rm = FALSE)
  orig1$normstage <- "Before normalization"

  # combine the two data sets
  plot_data <- rbind(orig1, norm1)
  colnames(plot_data) <- c("prot", "sample", "intensity", "normstage")
  plot_data$group <- factor(sapply(
    strsplit(
      as.character(plot_data[, "sample"]),
      "_"
    ),
    getElement, 1
  ))
  plot_data$normstage <- factor(plot_data$normstage,
    levels = c(
      "Before normalization",
      "After normalization"
    )
  )
  # Make density plots
  if (type == "density") {
    n_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = intensity,
        color = sample
      )
    ) +
      ggplot2::geom_density(lwd = text_size * 0.02) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      viridis::scale_color_viridis(
        discrete = TRUE,
        option = palette,
        begin = 0.3,
        end = 0.7,
        direction = -1
      ) +
      promor_facet_theme() +
      ggplot2::theme(
        legend.position = "none",
        axis.text = element_text(size = text_size * 0.7)
      ) +
      ggplot2::facet_wrap(~normstage)


    # Default: Boxplots
  } else {
    n_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = stats::reorder(
          sample,
          group
        ),
        y = intensity,
        fill = group
      )
    ) +
      geom_boxplot(
        color = "grey30",
        alpha = 0.7,
        outlier.shape = 1,
        outlier.stroke = 0.1,
        outlier.size = text_size * 0.04,
        outlier.color = "grey30",
        lwd = text_size * 0.02
      ) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~normstage) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      viridis::scale_fill_viridis(
        discrete = TRUE,
        option = palette,
        begin = 0.3,
        end = 0.7
      ) +
      promor_facet_theme() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text = element_text(size = text_size * 0.7)
      )
  }

  #Set temporary file_path if not specified
  if(is.null(file_path)){
    file_path <- tempdir()
  }

  if (save == TRUE) {
    ggsave(paste0(file_path, "/", file_name, ".", file_type),
      n_plot,
      dpi = dpi,
      width = plot_width,
      height = plot_height
    )
    return(n_plot)
  } else {
    return(n_plot)
  }
}
