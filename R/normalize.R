# Normalize intensities ----------------------------------------------------
#' Normalize intensity data
#' @author Chathurani Ranathunge
#' @description This function normalizes data using a user-specified
#' normalization method.
#' @import limma
#'
#' @param df An \code{imp.df} object with missing values imputed using
#' \code{impute_na}.
#' @param method Name of the normalization method to use. Choices are
#' \code{"none", "scale", "quantile" or "cyclicloess."}
#' Default is \code{"quantile."}
#'
#' @details \code{normalize_data} normalizes intensity values to achieve
#' consistency among samples. The function assumes that the intensities in the
#' data frame have been log-transformed, therefore, it is important to make sure
#' that \code{create_df} was run with \code{log.tr = TRUE} when creating the
#' \code{raw_df} object.
#'
#' @return A \code{norm.df} object, which is a data frame with
#' normalized intensities.
#'
#' @seealso \itemize{\item \code{create_df}
#' \item \code{impute_na}
#' \item See \code{\link[limma]{normalizeBetweenArrays}} in the R package
#' \code{limma} for more information on the different normalization methods
#' available.}
#'
#'
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(file.path = "./proteinGroups.txt", log.tr = TRUE)
#'
#' ## Normalize a data set.
#' raw_nm <- normalize_data(raw, method = "cyclicloess")
#'
#' ## Normalize an imputed data set.
#' raw_imp <- impute_na(raw)
#' raw_nm <- normalize_data(raw_imp)
#' }
#'
#' @export

normalize_data <- function(df,
                           method = "quantile") {
  norm_df <- limma::normalizeBetweenArrays(df,
    method = method
  )

  return(norm_df)
}

# Visualize normalization effects -----------------------------------------
#' Visualize the effect of normalization
#' @author Chathurani Ranathunge
#' @description This function visualizes the impact of normalization on
#' the data with plots.
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import viridis


#' @param original A \code{raw_df} object (output of \code{\link{create_df}})
#' containing missing values or an \code{imp.df} object after imputing the
#' missing values with \code{impute_na}.
#' @param normalized A \code{norm.df} object after normalizing the data frame
#' provided as \code{original}.
#' @param type Type of plot to generate. Choices are "box" or "density." Default
#' is \code{"box."}
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_name file_name File name to save the plot.
#' Default is \code{"Norm_plot."}
#' @param file_type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot_width Width of the plot. Default is \code{7}.
#' @param plot_height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details Given two data frames, one with data prior to normalization
#' (\code{original}) and one after normalization (\code{normalized}),
#' \code{norm_plot} generates side-by-side plots to visualize the effect of
#' normalization on the intensity data.
#'
#' @seealso \itemize{\item \code{\link{normalize_data}}
#' \item \code{create_df}
#' \item \code{impute_na}
#' }
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(file.path = "./proteinGroups.txt", log.tr = TRUE)
#'
#' ## Normalize an imputed data set.
#' raw_imp <- impute_na(raw)
#' raw_nm <- normalize_data(raw_imp)
#'
#' ## Visualize normalization with box plots.
#' norm_plot(raw_imp, raw_nm)
#'
#' ## Visualize normalization with density plots.
#' norm_plot(raw_imp, raw_nm, type = "density")
#' }
#'
#' @export
norm_plot <- function(original,
                      normalized,
                      type = "box",
                      text_size = 10,
                      palette = "viridis",
                      save = FALSE,
                      file_name = "Norm_plot",
                      file_type = "pdf",
                      dpi = 80,
                      plot_width = 7,
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
      promor_facet_theme +
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
      promor_facet_theme +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text = element_text(size = text_size * 0.7)
      )
  }
  if (save == TRUE) {
    ggsave(paste0(file_name, ".", file_type),
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
