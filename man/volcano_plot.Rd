% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findDE.R
\name{volcano_plot}
\alias{volcano_plot}
\title{Volcano plot}
\usage{
volcano_plot(
  fit_df,
  adj_method = "BH",
  sig = "adjP",
  cutoff = 0.05,
  lfc = 1,
  line_fc = TRUE,
  line_p = TRUE,
  palette = "viridis",
  text_size = 10,
  label_top = FALSE,
  n_top = 10,
  save = FALSE,
  file_path = NULL,
  file_name = "Volcano_plot",
  file_type = "pdf",
  plot_height = 7,
  plot_width = 7,
  dpi = 80
)
}
\arguments{
\item{fit_df}{A \code{fit_df} object from performing \code{find_dep}.}

\item{adj_method}{Method used for adjusting the p-values for multiple
testing. Default is \code{"BH"}.}

\item{sig}{Criteria to denote significance. Choices are \code{"adjP"}
(default) for adjusted p-value or \code{"P"} for p-value.}

\item{cutoff}{Cutoff value for p-values and adjusted p-values. Default is
0.05.}

\item{lfc}{Minimum absolute log2-fold change to use as threshold for
differential expression.}

\item{line_fc}{Logical. If \code{TRUE}(default), a dotted line will be shown
to indicate the \code{lfc} threshold in the plot.}

\item{line_p}{Logical. If \code{TRUE}(default), a dotted line will be shown
to indicate the p-value or adjusted p-value \code{cutoff.}}

\item{palette}{Viridis color palette option for plots. Default is
\code{"viridis"}. See
\code{\link[viridis:viridis]{viridis}}
for available options.}

\item{text_size}{Text size for axis text, labels etc.}

\item{label_top}{Logical. If \code{TRUE} (default), labels are added to the
dots to indicate protein names.}

\item{n_top}{The number of top hits to label with protein name when
\code{label_top = TRUE.} Default is \code{10}.}

\item{save}{Logical. If \code{TRUE} saves a copy of the plot in the
directory provided in \code{file_path}.}

\item{file_path}{A string containing the directory path to save the file.}

\item{file_name}{File name to save the plot. Default is "Volcano_plot."}

\item{file_type}{File type to save the plot. Default is \code{"pdf".}}

\item{plot_height}{Height of the plot. Default is 7.}

\item{plot_width}{Width of the plot. Default is 7.}

\item{dpi}{Plot resolution. Default is \code{80.}}
}
\value{
A \code{ggplot2} plot object.
}
\description{
This function generates volcano plots to visualize
differentially expressed proteins between groups.
}
\details{
\itemize{
\item Volcano plots show log-2-fold change on the x-axis,
and based on the significance criteria chosen, either -log10(p-value) or
-log10(adjusted p-value) on the y-axis.
\item \code{volcano_plot} requires a \code{fit_df} object from performing
differential expression analysis with \code{find_dep.}
\item User has the option to choose criteria that denote significance.
}
}
\examples{

## Create a volcano plot with default settings.
volcano_plot(ecoli_fit_df)

## Change significance criteria and cutoff
volcano_plot(ecoli_fit_df, cutoff = 0.1, sig = "P")

## Label top 30 differentially expressed proteins and
## change the color palette of the plot
volcano_plot(ecoli_fit_df, label_top = TRUE, n_top = 30, palette = "mako")

}
\seealso{
\itemize{
\item \code{find_dep}
\item \code{\link[limma:topTable]{topTable}} and
\code{\link[limma:lmFit]{lmFit}} functions from the
\code{\link[limma]{limma}} package.
}
}
\author{
Chathurani Ranathunge
}
