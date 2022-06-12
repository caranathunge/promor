# Scatter plots : correlation among technical replicates ----------------
#' Correlation between technical replicates
#' @author Chathurani Ranathunge
#' @description This function generates scatter plots to visualize the
#' correlation between a given pair of technical replicates (Eg: 1 vs 2)
#' for each sample.
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import gridExtra
#' @import viridis
#'
#' @param df A \code{raw_df} object (output of \code{\link{create_df}})
#' containing technical replicates.
#' @param rep1 Numerical. Technical replicate number.
#' @param rep2 Numerical. Number of the second technical replicate to compare
#' to \code{rep1}.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' working directory.
#' @param file_type File type to save the scatter plots.
#' Default is \code{"pdf"}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param nrow Numerical. Number of plots to print in a row in a single page.
#' Default is \code{4}.
#' @param ncol Numerical. Number of plots to print in a column in a single
#' page. Default is \code{4}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details
#' \itemize{\item Given a data frame of log-transformed intensities
#' (a \code{raw_df} object) and a pair of numbers referring to the technical
#' replicates, \code{corr_plot} produces a list of scatter plots showing
#' correlation between the given pair of technical replicates for all the
#' samples provided in the data frame.
#' \item Note: \code{corr_plot} assumes that the sample names (columns) in the
#' data frame follow the "Group_UniqueSampleID_TechnicalReplicate" notation.
#' For example,
#' "technical replicate 01" of "sample 10" belonging to "group A" should be
#' labeled as "A_10_01" in the data frame.
#' \item Note: \code{nrow} * \code{ncol} should be equal to the number of
#' samples to display in a single page.
#' }
#'
#' @seealso
#' \code{create_df}
#'
#' @examples
#' \dontrun{
#'
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(file.path = "./proteinGroups.txt", log.tr = TRUE)
#'
#' ## Compare technical replicates 1 vs. 2 for all samples.
#' corr_plot(raw, 1, 2, nrow = 4, ncol = 6)
#' }
#'
#' @return
#' A list of \code{ggplot2} plot objects.
#'
#'
#'
#'
#' @export
corr_plot <- function(df,
                      rep1,
                      rep2,
                      save = FALSE,
                      file_type = "pdf",
                      palette = "viridis",
                      text_size = 5,
                      nrow = 4,
                      ncol = 4,
                      dpi = 80) {

  # Separate the sample id from the column names, remove duplicates, paste _ and
  # create a list. This is now a list of unique sample names used in the
  # dataframe. Eg: _12_ instead of 12 because '12' can match 120, 121, 212 etc.
  sample_name <- unlist(paste0("_", unique(sapply(
    strsplit(colnames(df), "_"),
    "[", 2
  )),
  "_",
  sep = ""
  ))

  # split the dataframe by each sample and output a list of dataframes for
  # plotting. 'as.data.frame' forces even those with just one column or
  #replicate into a dataframe rather than a row
  sub_df <- lapply(
    sample_name,
    function(y) as.data.frame(df[, grepl(y, names(df))])
  )

  # Only keep the samples with at least 2 replicates for plotting TR1 v TR2
  if (rep1 == 3 || rep2 == 3) {
    plot_data <- Filter(function(c) ncol(c) == 3, sub_df)
  } else {
    plot_data <- Filter(function(c) ncol(c) >= 2, sub_df)
  }

  # Set colors
  pal_col <- set_col(palette, 1)

  # Create a list of scatter plots and print/save
  plot_list <- lapply(seq_along(plot_data), function(t) {
    ggplot2::ggplot(
      plot_data[[t]],
      ggplot2::aes(x = plot_data[[t]][, rep1], y = plot_data[[t]][, rep2])
    ) +
      ggplot2::geom_point(col = pal_col, size = text_size * 0.1) +
      ggplot2::geom_text(
        label = sapply(
          strsplit(rownames(plot_data[[t]]), ";"),
          getElement, 1
        ),
        hjust = "inward",
        vjust = "inward",
        size = text_size * 0.2,
        check_overlap = TRUE
      ) +
      promor::promor_facet_theme +
      ggplot2::theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(
          size = text_size,
          face = "bold",
          colour = "grey20"
        )
      ) +
      ggplot2::ggtitle(gsub("\\_\\d+$", "", colnames(plot_data[[t]][1])))
  })
  if (save == TRUE) {
    ggplot2::ggsave(paste0("TR", rep1, "vs", "TR", rep2, ".", file_type),
      marrangeGrob(
        grobs = plot_list,
        nrow = nrow,
        ncol = ncol,
        top = ""
      ),
      dpi = dpi
    )
  } else {
    grid.arrange(grobs = plot_list, newpage = TRUE)
  }
}


# Remove user-specified samples -------------------------------------------
#' Remove user-specified samples
#' @author Chathurani Ranathunge
#' @description This function removes user-specified samples from the
#' data frame.
#'
#' @param df A \code{raw_df} object.
#' @param rem Name of the sample to remove.
#'
#' @details \itemize{\item\code{rem_sample} assumes that sample names follow the
#' "Group_Sample_TechnicalReplicate" notation.
#' \item If all the technical replicates representing a sample needs to be
#' removed, provide "Group_UniqueSampleID" as \code{rem}.
#' \item If a specific technical replicate needs to be removed in case it
#' shows weak correlation with other technical replicates for example, you can
#' remove that particular technical replicate by providing
#' "Group_UniqueSampleID_TechnicalReplicate" as \code{rem}.
#' }
#' @return A \code{raw_df} object.
#'
#' @seealso \code{\link{corr_plot}}, \code{\link{create_df}}
#'
#' @examples
#' \dontrun{
#'
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(file.path = "./proteinGroups.txt", log.tr = TRUE)
#'
#' ## Remove all technical replicates of "A_10"
#' raw_1 <- rem_sample(raw, "A_10")
#'
#' ## Remove only technical replicate number 2 of "A_10"
#' raw_1 <- rem_sample(raw, "A_10_2")
#' }
#'
#' @export
rem_sample <- function(df, rem) {
  df_rem <- df[, -grep(rem, colnames(df), fixed = TRUE)]
  return(df_rem)
  message(paste0("Sample ", rem, "has been removed."))
}


# Compute average intensity across tech.replicates for each sample --------
#' Compute average intensity
#' @author Chathurani Ranathunge
#'
#' @description This function computes average intensities across
#' technical replicates for each sample.
#'
#' @import limma
#'
#' @param df A \code{raw_df} object containing technical replicates.
#'
#' @details
#' \code{aver_techreps} assumes that sample names in the data frame
#' follow the "Group_UniqueSampleID_TechnicalReplicate" notation.
#'
#' @return A \code{raw_df} object of averaged intensities.
#'
#' @seealso \code{\link[limma]{avearrays}} from \code{\link[limma]{limma}}
#' package.
#'
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(file.path = "./proteinGroups.txt", log.tr = TRUE)
#'
#' # Compute average intensities across technical replicates.
#' raw_1 <- aver_techreps(raw)
#' }
#'
#' @export
aver_techreps <- function(df) {
  # Convert dataframe back to a matrix as dataframes don't allow multiple
  # columns with the same name. We need that feature to average over TRs.
  df_mat <- as.matrix(df)

  # substitute technical replicates with the sample name in the vector.
  # Now in each sample, all technical replicates are labelled the same way.
  colnames(df_mat) <- gsub("\\_\\d+$", "", colnames(df_mat))

  # Average across technical replicates
  df_ave <- limma::avearrays(df_mat,
    ID = colnames(df_mat),
    weights = NULL
  )
  df_ave <- as.data.frame(df_ave)
  return(df_ave)
}
