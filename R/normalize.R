# Normalize intensities ----------------------------------------------------
#' Normalize intensity data
#' @author Chathurani Ranathunge
#' @description This function normalizes data using a user-specified
#' normalization method.
#' @import limma
#'
#' @param df A \code{raw.df} object, or ideally, an \code{imp.df} object with
#' missing values imputed.
#' @param method Name of the normalization method to use. Choices are
#' \code{"none", "scale", "quantile" or "cyclicloess."}
#' Default is \code{"quantile."}
#'
#' @details \code{normalize.data} normalizes intensity values to achieve
#' consistency among samples. The function assumes that the intensities in the
#' data frame have been log-transformed, therefore, it is important to make sure
#' that \code{create.df} was run with \code{log.tr = TRUE} when creating the
#' \code{raw.df} object.
#'
#' @return A \code{norm.df} object, which is a data frame with
#' normalized intensities.
#'
#' @seealso \itemize{\item \code{create.df}
#' \item \code{impute.NA}
#' \item See \code{\link[limma]{normalizeBetweenArrays}} in the R package
#' \code{limma} for more information on the different normalization methods
#' available.}
#'
#'
#' @examples
#' \dontrun{
#' ## Create a raw.df object from a proteinGroups.txt file.
#' raw <- create.df(file.path = "./proteinGroups.txt", log.tr = TRUE)
#'
#' ## Normalize a data set.
#' raw_nm <- normalize.data(raw, method = "cyclicloess")
#'
#' ## Normalize an imputed data set.
#' raw_imp <- impute.NA(raw)
#' raw_nm <- normalize.data(raw_imp)
#'
#' }
#'
#'

#' @export

normalize.data <- function(df,
                      method = "quantile"){
  norm_df <- limma::normalizeBetweenArrays(df,
                                           method = method)
  return(norm_df)
}

# Visualize normalization effects -----------------------------------------
#' Visualize the effect of normalization
#' @author Chathurani Ranathunge
#' @description This function visualizes the impact of normalization on
#' the data with plots.
#' @importFrom reshape2 melt
#' @import ggplot2

#' @param original A \code{raw.df} object (output of \code{\link{create.df}})
#' containing missing values or an \code{imp.df} object after imputing the
#' missing values with \code{impute.NA}.
#' @param normalized A \code{norm.df} object after normalizing the data frame
#' provided as \code{original}.
#' @param type Type of plot to generate. Choices are "box" or "density." Default
#' is \code{"box."}
#' @param text.size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param palette Color palette for box plots. Default is \code{"YlGnBu."}
#' @param save Logical. If \code{TRUE} (default) saves a copy of the plot in the
#' working directory.
#' @param file.name file.name File name to save the plot.
#' Default is \code{"Norm_plot."}
#' @param file.type File type to save the plot.
#' Default is \code{"pdf"}.
#' @param plot.width Width of the plot. Default is \code{7}.
#' @param plot.height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details Given two data frames, one with data prior to normalization
#' (\code{original}) and one after normalization (\code{normalized}),
#' \code{norm.plot} generates side-by-side plots to visualize the effect of
#' normalization on the intensity data.
#'
#' @seealso \itemize{\item \code{\link{normalize.data}}
#' \item \code{create.df}
#' \item \code{impute.NA}
#' }
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' ## Create a raw.df object from a proteinGroups.txt file.
#' raw <- create.df(file.path = "./proteinGroups.txt", log.tr = TRUE)
#'
#' ## Normalize an imputed data set.
#' raw_imp <- impute.NA(raw)
#' raw_nm <- normalize.data(raw_imp)
#'
#' ## Visualize normalization with box plots.
#' norm.plot(raw_imp, raw_nm)
#'
#' ## Visualize normalization with density plots.
#' norm.plot(raw_imp, raw_nm, type = "density")
#' }
#'
#' @export
norm.plot <- function(original,
                      normalized,
                      type = "box",
                      text.size = 10,
                      palette = "YlGnBu",
                      save = TRUE,
                      file.name = "Norm_plot",
                      file.type = "pdf",
                      dpi = 80,
                      plot.width = 7,
                      plot.height = 7){

  #Set global variables to null
  intensity <- value <- group <- NULL

  #Pre-prossesing data for plotting
  normalized <- as.matrix(normalized)
  norm1 <- reshape2::melt(normalized, na.rm = FALSE)
  norm1$normstage <- "After normalization"

  original <- as.matrix(original)
  orig1 <- reshape2::melt(original, na.rm =FALSE)
  orig1$normstage <- "Before normalization"

  #combine the two data sets
  plot_data<- rbind(orig1, norm1)
  colnames(plot_data) <- c("prot", "sample", "intensity", "normstage")
  plot_data$group<- sapply(strsplit(as.character(plot_data[,"sample"]),'_'),
                           getElement,1)
  plot_data$normstage <- factor(plot_data$normstage,
                                levels = c("Before normalization",
                                           "After normalization"))

  if(type == "density"){
    norm_plot <- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = intensity,
                                         color = sample)) +
      ggplot2::geom_density(lwd = text.size * 0.02)+
      ggplot2::xlab("") +
      ggplot2::ylab("")+
      #ggplot2::scale_fill_brewer(palette = palette)+
      ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "none",
            panel.border = element_rect(fill = NA,
                                        colour = "grey",
                                        size = 0.5),
            text= element_text(size = text.size),
            axis.line.x = element_line(size = 0.1),
            axis.line.y = element_line(size = 0.1),
            axis.ticks.x = element_line(size = 0.1),
            axis.ticks.y = element_line(size = 0.1),
            axis.title.x = element_text(size = text.size),
            axis.title.y= element_text(size = text.size),
            axis.line = element_line(colour = "grey",
                                     size = 0.5),
            strip.text = element_text(size= text.size),
            strip.background = element_rect(fill = "grey95",
                                            colour = "grey",
                                            size = 0.5))+
      ggplot2::facet_wrap( ~normstage)


#Default: Boxplots

  }else{
    norm_plot <- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = sample,
                                         y = intensity,
                                         fill = group)) +
      geom_boxplot(color = "grey30",
                   alpha = 0.9,
                   outlier.shape = 1,
                   outlier.stroke = 0.1,
                   outlier.size = text.size * 0.04,
                   outlier.color = "grey30",
                   lwd = text.size * 0.02)+
      ggplot2::coord_flip()+
      ggplot2::facet_wrap( ~normstage)+
      ggplot2::xlab("") +
      ggplot2::ylab("")+
      ggplot2::scale_fill_brewer(palette = palette)+
      ggplot2::theme_classic()+
      ggplot2::theme(panel.border = element_rect(fill = NA,
                                                 colour = "grey",
                                                 size = 0.5),
            text = element_text(size = text.size),
            legend.title = element_blank(),
            axis.line.x = element_line(size = 0.1),
            axis.line.y = element_line(size = 0.1),
            axis.ticks.x = element_line(size = 0.1),
            axis.ticks.y = element_line(size = 0.1),
            axis.title.x = element_text(size = text.size),
            axis.title.y= element_text(size = text.size),
            axis.line = element_line(colour = "grey",
                                     size = 0.5),
            strip.text = element_text(size= text.size),
            strip.background = element_rect(fill = "grey95",
                                            colour = "grey",
                                            size = 0.5))

  }
  if(save == TRUE){
    ggsave(paste0(file.name,".", file.type),
           norm_plot,
           dpi = dpi,
           width = plot.width,
           height = plot.height)
    return(norm_plot)
    }else{
      return(norm_plot)
    }
  }
