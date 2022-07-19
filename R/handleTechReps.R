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
#' @param raw_df A \code{raw_df} object (output of \code{\link{create_df}})
#' containing technical replicates.
#' @param rep_1 Numerical. Technical replicate number.
#' @param rep_2 Numerical. Number of the second technical replicate to compare
#' to \code{rep1}.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' directory provided in \code{file_path}.
#' @param file_type File type to save the scatter plots.
#' Default is \code{"pdf"}.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridisLite:viridis]{viridis}}
#' for available options.
#' @param text_size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param n_row Numerical. Number of plots to print in a row in a single page.
#' Default is \code{4}.
#' @param n_col Numerical. Number of plots to print in a column in a single
#' page. Default is \code{4}.
#' @param dpi Plot resolution. Default is \code{80}.
#' @param file_path A string containing the directory path to save the file.
#'
#' @details
#' \itemize{\item Given a data frame of log-transformed intensities
#' (a \code{raw_df} object) and a pair of numbers referring to the technical
#' replicates, \code{corr_plot} produces a list of scatter plots showing
#' correlation between the given pair of technical replicates for all the
#' samples provided in the data frame.
#' \item Note: \code{n_row} * \code{n_col} should be equal to the number of
#' samples to display in a single page.
#' }
#'
#' @seealso
#' \code{create_df}
#'
#' @examples
#' \donttest{
#'
#' ## Use a data set containing technical replicates to create a raw_df object
#' raw_df <- create_df(
#' prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg2.txt",
#' exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed2.txt",
#' tech_reps = TRUE
#' )
#' }
#' \dontrun{
#' ## Compare technical replicates 1 vs. 2 for all samples and save plots in a pdf in
#' ## the working directory.
#' corr_plot(raw_df, rep_1 = 1, rep_2 = 2, n_row = 3, n_col = 2, save = TRUE, file_path = ".")
#' }
#'
#' @return
#' A list of \code{ggplot2} plot objects.
#'
#'
#'
#'
#' @export
corr_plot <- function(raw_df,
                      rep_1,
                      rep_2,
                      save = FALSE,
                      file_type = "pdf",
                      palette = "viridis",
                      text_size = 5,
                      n_row = 4,
                      n_col = 4,
                      dpi = 80,
                      file_path = NULL) {

  # Separate the sample id from the column names, remove duplicates, paste _ and
  # create a list. This is now a list of unique sample names used in the
  # dataframe. Eg: _12_ instead of 12 because '12' can match 120, 121, 212 etc.
  sample_name <- unlist(paste0("_", unique(sapply(
    strsplit(colnames(raw_df), "_"),
    "[", 2
  )),
  "_",
  sep = ""
  ))

  # split the dataframe by each sample and output a list of dataframes for
  # plotting. 'as.data.frame' forces even those with just one column or
  # replicate into a dataframe
  sub_df <- lapply(
    sample_name,
    function(y) as.data.frame(raw_df[, grepl(y, names(raw_df))])
  )

  # Only keep the samples with at least 2 replicates for plotting TR1 v TR2
  if (rep_1 == 3 || rep_2 == 3) {
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
      ggplot2::aes(x = plot_data[[t]][, rep_1], y = plot_data[[t]][, rep_2])
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
      promor_facet_theme() +
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

  #Set temporary file_path if not specified
  if(is.null(file_path)){
    file_path <- tempdir()
  }

  if (save == TRUE) {
    ggplot2::ggsave(paste0(file_path, "/TR", rep_1, "vs", "TR", rep_2, ".", file_type),
      marrangeGrob(
        grobs = plot_list,
        nrow = n_row,
        ncol = n_col,
        top = paste0("TR", rep_1, " vs ", "TR", rep_2)
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
#' @param raw_df A \code{raw_df} object.
#' @param rem Name of the sample to remove.
#'
#' @details \itemize{\item\code{rem_sample} assumes that sample names follow the
#' "Group_UniqueSampleID_TechnicalReplicate" notation (Use \code{head(raw_df)}
#' to see the structure of the \code{raw_df} object.)
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
#' \donttest{
#'
#' ## Use a data set containing technical replicates to create a raw_df object
#' raw_df <- create_df(
#' prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg2.txt",
#' exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed2.txt",
#' tech_reps = TRUE
#' )
#' # Check the first few rows of the raw_df object
#' head(raw_df)
#'
#' ## Remove all technical replicates of "WT_4"
#' raw_df1 <- rem_sample(raw_df, "WT_4")
#'
#' ## Remove only technical replicate number 2 of "WT_4"
#' raw_df2 <- rem_sample(raw_df, "WT_4_2")
#' }
#'
#' @export
rem_sample <- function(raw_df, rem) {
  df_rem <- raw_df[, -grep(rem, colnames(raw_df), fixed = TRUE)]
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
#' @param raw_df A \code{raw_df} object containing technical replicates.
#'
#' @details
#' \code{aver_techreps} assumes that column names in the data frame
#' follow the "Group_UniqueSampleID_TechnicalReplicate" notation. (Use
#' \code{head(raw_df)} to see the structure of the \code{raw_df} object.)
#'
#' @return A \code{raw_df} object of averaged intensities.
#'
#' @seealso \itemize{\item\code{\link[limma]{avearrays}} from \code{\link[limma]{limma}}
#' package.\item \code{\link{create_df}}
#' }
#'
#' @examples
#' \donttest{
#'
#' ## Use a data set containing technical replicates to create a raw_df object
#' raw_df <- create_df(
#' prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg2.txt",
#' exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed2.txt",
#' tech_reps = TRUE
#' )
#'
#' # Compute average intensities across technical replicates.
#' rawdf_ave <- aver_techreps(raw_df)
#' }
#'
#' @export
aver_techreps <- function(raw_df) {
  # Convert dataframe back to a matrix as dataframes don't allow multiple
  # columns with the same name. We need that feature to average over TRs.
  df_mat <- as.matrix(raw_df)

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
