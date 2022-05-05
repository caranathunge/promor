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
#'
#' @param df A \code{raw.df} object (output from \code{\link{create.df}}).
#' @param protein.range The range of proteins to plot. Default is \code{ALL},
#' meaning all the proteins in the data frame.
#' @param sample.range The range of samples to plot. Default is \code{ALL},
#' meaning all the samples in the data frame.
#' @param reorder.x Logical. If \code{TRUE} samples on the x axis are reordered
#' using the function given in \code{x.FUN}. Default is \code{FALSE}.
#' @param reorder.y Logical. If \code{TRUE} proteins in the y axis are reordered
#' using the function given in \code{y.FUN}. Default is \code{FALSE}.
#' @param x.FUN Function to reorder samples along the x axis. Possible options
#' include \code{mean} and \code{sum}. Default is \code{mean}.
#' @param y.FUN Function to reorder proteins along the y axis. Possible options
#' include \code{mean} and \code{sum}. Default is \code{mean}.
#' @param col.na Color assigned for missing values (NAs). Default is "Black."
#' @param col.val Color assigned for valid values. Default is "Grey."
#' @param text.size Text size for axis labels. Default is \code{10}.
#' @param save Logical. If \code{TRUE} (default) saves a copy of the plot in the
#' working directory.
#' @param file.type File type to save the heatmap. Default is \code{"pdf"}.
#' @param file.name File name to save the heatmap.
#' Default is \code{"MissingData_heatmap"}.
#' @param plot.width Width of the plot. Default is \code{15}.
#' @param plot.height Height of the plot. Default is \code{15}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details
#' This function visualizes patterns of missing value occurrence using a
#' heatmap. User can choose to reorder the axes using the available functions
#' (\code{x.FUN}, \code{y.FUN}) to better understand the underlying cause of
#' missing data.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @seealso \code{\link{create.df}}
#'
#' @examples
#' \dontrun{
#' ## Create a raw.df object from a proteinGroups.txt file.
#' raw <- create.df(file.path = "./proteinGroups.txt")
#'
#' ## Missing data heatmap with default settings.
#' heatmap.NA(raw)
#'
#' ## Missing data heatmap with x and y axes reordered by sum of intensity.
#' heatmap.NA(raw, reorder.x = TRUE, reorder.y = TRUE, x.FUN = sum,
#' y.FUN = sum)
#' }
#'
#' @export
heatmap.NA <- function(df,
                       protein.range = "ALL",
                       sample.range = "ALL",
                       reorder.x = FALSE,
                       reorder.y = FALSE,
                       x.FUN = mean,
                       y.FUN = mean,
                       col.na = "black",
                       col.val = "grey",
                       text.size = 10,
                       save = TRUE,
                       file.type = "pdf",
                       file.name = "MissingData_heatmap",
                       plot.width = 15,
                       plot.height = 15,
                       dpi = 80){

 #Binding global variables to the local function
  value <- protgroup <- NULL

  #select the range of data to plot
  if(protein.range == "ALL" && sample.range != "ALL"){
    df <- as.matrix(df[,sample.range])
  }else if(protein.range != "ALL" && sample.range == "ALL"){
    df <- as.matrix(df[protein.range,])
  }else if(protein.range != "ALL" && sample.range != "ALL"){
    df <- as.matrix(df[protein.range,sample.range])
  }else{
    df <- as.matrix(df)
  }

  #Convert the data into long format for plotting and make necessary changes
  #i.e. adding column headers etc.
  hmap_data <- reshape2::melt(df, na.rm = FALSE)
  hmap_data$mjprot <- sapply(strsplit(as.character(hmap_data[,1]),';'),
                             getElement,1)
  hmap_data[1] <- NULL
  hmap_df <- as.data.frame(hmap_data)
  colnames(hmap_df) <- c( "sample",
                          "value",
                          "protgroup")

  #Options for arranging the rows and the columns of the heat map
  if(reorder.x == TRUE & reorder.y == TRUE){
    hmap <- ggplot2::ggplot(hmap_df,
                            ggplot2::aes(x = stats::reorder(sample,
                                                     value,
                                                     na.rm = TRUE,
                                                     FUN = x.FUN),
                                         y = stats::reorder(protgroup,
                                                     value,
                                                     na.rm = TRUE,
                                                     FUN = y.FUN),
                                         fill = value))

    }else if(reorder.x==FALSE & reorder.y == TRUE){
    hmap <- ggplot2::ggplot(hmap_df,
                            ggplot2::aes(x = sample,
                                         y = stats::reorder(protgroup,
                                                     value,
                                                     na.rm = TRUE,
                                                     FUN = y.FUN),
                                         fill = value))

    }else if(reorder.x== TRUE & reorder.y == FALSE){
    hmap <- ggplot2::ggplot(hmap_df,
                            ggplot2::aes(x = stats::reorder(sample,
                                                     value,
                                                     na.rm = TRUE,
                                                     FUN = x.FUN),
                                         y = protgroup,
                                         fill = value))

    }else{

    hmap <- ggplot2::ggplot(hmap_df,
                            ggplot2::aes(x = sample,
                                         y = protgroup,
                                         fill = value))
    }

  #Create heat map
  hmap <- hmap +
    ggplot2::geom_tile()+
    ggplot2::scale_fill_gradient(high = col.val,
                                 low = col.val,
                                 na.value = col.na)+
    ggplot2::theme(plot.background = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90,
                                              size = text.size),
                   axis.text.y = element_text(size = text.size),
                   axis.ticks.x = element_blank(),
                   legend.position = "none")+
    ggplot2::xlab("") + ggplot2::ylab("")

#Save the heatmap as a pdf
if(save == TRUE){
  ggplot2::ggsave(paste0(file.name,".",file.type),
                  hmap,
                  dpi = dpi,
                  width = plot.width,
                  height = plot.height)
  return(hmap)
  }else{
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
#' @param df A \code{raw.df} object (output of \code{\link{create.df}})
#' containing missing values.
#' @param method Imputation method to use. Default is \code{"minProb"}.
#' Available methods: \code{"minDet", "RF", "kNN", and "SVD"}.
#' @param tune.sigma A scalar used in the \code{"minProb"} method for
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
#' @param nPcs Number of principal components to calculate when using the
#' \code{"SVD"} method. Default is 2.
#'
#' @details \code{impute.NA} function imputes missing values using a
#' user-specified imputation method from the available options, \code{minProb},
#' \code{minDet}, \code{kNN}, \code{RF}, and \code{SVD}.
#'
#' @seealso More information on the available imputation methods can be found
#' in their respective packages.
#' \itemize{\item \code{\link{create.df}}
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
#' ## Create a raw.df object from a proteinGroups.txt file.
#' raw <- create.df(file.path = "./proteinGroups.txt")
#'
#' ## Impute missing values in the data frame using the default minProb
#' method.
#' imp_raw <- impute.NA(raw)
#'
#' ## Using the RF method with the number of iterations set at 5,
#' and the number of trees set at 100.
#' imp_raw <- impute.NA(raw, method = "RF", maxiter = 5, ntree = 100 )
#'
#' ## Using the kNN method.
#' imp_raw <- impute.NA(raw, method = "kNN")
#'
#' ## Using the SVD method with nPCs set to 3.
#' imp_raw <- impute.NA(raw, method = "SVD", nPCs = 3)
#'
#' ## Using the minDet method with q set at 0.001.
#' imp_raw <- impute.NA(raw, method = "minDet", q = 0.001)
#'
#' }
#'
#' @export
impute.NA <- function(df,
                      method = "minProb",
                      tune.sigma = 1,
                      q = 0.01,
                      maxiter = 10,
                      ntree = 20,
                      nPcs = 2){

# Setting functions to NULL
impute.Min.Prob <- NULL
impute.MinDet <- NULL
value <- protgroup <- NULL


#Run imputeLCMD function imputeMinProb
  impute.Min.Prob <<- function (dataSet.mvs, q = 0.01, tune.sigma = 1){
    nSamples = dim(dataSet.mvs)[2]
    nFeatures = dim(dataSet.mvs)[1]
    dataSet.imputed = dataSet.mvs
    min.samples = apply(dataSet.imputed, 2, quantile, prob = q,
                        na.rm = T)
    count.NAs = apply(!is.na(dataSet.mvs), 1, sum)
    count.NAs = count.NAs/nSamples
    dataSet.filtered = dataSet.mvs[which(count.NAs > 0.5), ]
    protSD = apply(dataSet.filtered, 1, sd)
    sd.temp = median(protSD, na.rm = T) * tune.sigma
    #print(sd.temp)
    for (i in 1:(nSamples)) {
      dataSet.to.impute.temp = rnorm(nFeatures,
                                     mean = min.samples[i],
                                     sd = sd.temp)
      dataSet.imputed[
        which(
          is.na(
            dataSet.mvs[, i])), i] = dataSet.to.impute.temp[which(
              is.na(dataSet.mvs[,i]))]
    }
    return(dataSet.imputed)
  }


#Run imputeLCMD function imputeMinDet
  impute.MinDet <<- function (dataSet.mvs, q = 0.01){
    nSamples = dim(dataSet.mvs)[2]
    dataSet.imputed = dataSet.mvs
    lowQuantile.samples = apply(dataSet.imputed, 2, quantile,
                                prob = q, na.rm = T)
    for (i in 1:(nSamples)) {
      dataSet.imputed[which(
        is.na(dataSet.mvs[, i])), i] = lowQuantile.samples[i]
    }
    return(dataSet.imputed)
  }

#Run the user-specified imputation method

  if (method == "minDet"){
    df_imputed_mindet <- impute.MinDet(df,
                                       q = q)
    return(df_imputed_mindet)

    }else if(method == "RF"){
    df_imp_temp <- missForest::missForest(df,
                                          maxiter = maxiter,
                                          ntree = ntree,
                                          verbose= TRUE)
    df_imputed_RF <- df_imp_temp$ximp
    return(df_imputed_RF)

    }else if (method == "kNN"){
    df_imputed_knn <- VIM::kNN(df, imp_var = FALSE)
    rownames(df_imputed_knn) <- rownames(df)
    return(df_imputed_knn)

    }else if (method == "SVD") {
    df[is.nan(df)] <- NA
    df_imp_temp <- pcaMethods::pca(object = df,
                                   method = "svdImpute",
                                   nPcs = nPcs,
                                   verbose= TRUE)
    df_imputed_svd <- completeObs(df_imp_temp)
    return(df_imputed_svd)

    }else if (method == "minProb"){
    df_imputed_minprob <- impute.Min.Prob(df,
                                         q = q,
                                         tune.sigma = tune.sigma)
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
#'
#' @param original A \code{raw.df} object (output of \code{\link{create.df}})
#' containing missing values.
#' @param imputed An \code{imp.df} object obtained from running \code{impute.NA}
#'  on the same data frame provided as \code{original}.
#' @param global Logical. If \code{TRUE} ({default}), a global density plot is
#' produced. If \code{FALSE}, sample-wise density plots are produced.
#' @param text.size Text size for plot labels, axis labels etc. Default is
#' \code{10}.
#' @param orig.col Fill color for the \code{original} data set.\cr
#' Default is "Blue."
#' @param imp.col Fill color color for the \code{imputed} data set.\cr
#' Default is "Red."
#' @param nrow Required if \code{global = FALSE} to indicate the number of rows
#' to print the plots.
#' @param ncol Required if \code{global = FALSE} to indicate the number of
#' columns to print the plots.
#' @param save Logical. If \code{TRUE} (default) saves a copy of the plot in the
#' working directory.
#' @param file.name file.name File name to save the density plot/s.
#' Default is \code{"Impute_plot."}
#' @param file.type File type to save the density plot/s.
#' Default is \code{"pdf"}.
#' @param plot.width Width of the plot. Default is \code{7}.
#' @param plot.height Height of the plot. Default is \code{7}.
#' @param dpi Plot resolution. Default is \code{80}.
#'
#' @details
#' \itemize{\item Given two data frames, one with missing values
#' (\code{raw.df} object) and the other, an imputed data frame
#' (\code{imp.df} object) of the same data set, \code{impute.plot}
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
#' ## Create a raw.df object from a proteinGroups.txt file.
#' raw <- create.df(file.path = "./proteinGroups.txt")
#'
#' ## Impute missing values in the data frame.
#' imp_raw <- impute.NA(raw)
#'
#' ## Visualize the impact of missing data imputation with a global density plot.
#' impute.plot(original = raw, imputed = imp_raw)
#'
#' ## Make sample-wise density plots for a data set of 25 samples.
#' impute.plot(raw, imp_raw, global = FALSE, nrow = 5, ncol = 5)
#'
#' }
#' @export
impute.plot <- function(original,
                            imputed,
                            global = TRUE,
                            text.size = 10,
                            orig.col = "blue",
                            imp.col = "red",
                            nrow,
                            ncol,
                            save = TRUE,
                            file.name = "Impute_plot",
                            file.type ="pdf",
                            plot.width = 7,
                            plot.height = 7,
                            dpi = 80){

  #Set global variables to null
  value <- NULL
  #Make necessary changes and combine data frames before plotting
  original <- as.matrix(original)
  orig_data <- reshape2::melt(original, na.rm = FALSE)
  orig_data$mjprot <- sapply(strsplit(as.character(orig_data[,1]),';'),
                              getElement,1)
  orig_data[1] <- NULL
  colnames(orig_data) <- c( "sample", "value", "protgroup")
  orig_data$stage <- "Before imputation"

  #do the same for imputed data set
  imputed <- as.matrix(imputed)
  imp_data <- reshape2::melt(imputed, na.rm = FALSE)
  imp_data$mjprot <- sapply(strsplit(as.character(imp_data[,1]),';'),
                              getElement,1)
  imp_data[1] <- NULL
  colnames(imp_data) <- c( "sample", "value", "protgroup")
  imp_data$stage <- "After imputation"
  dplot_data <- rbind(orig_data, imp_data)

  #Plot global density plot

  if(global ==TRUE){
    g_densityplot <- ggplot2::ggplot(dplot_data,
                                     ggplot2::aes(x=value)) +
      ggplot2::geom_density(ggplot2::aes(fill = factor(stage,
                                                       levels = c(
                                                         "Before imputation",
                                                         "After imputation"))),
                            alpha = 0.25,
                            lwd = 0.1 )+
      ggplot2::xlab("Intensity") +
      ggplot2::ylab("Density")+
      ggplot2::scale_fill_manual(values=c(orig.col,imp.col))+
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
                     legend.text = element_text(size = text.size))


    if(save == TRUE){
      ggplot2::ggsave(paste0(file.name,".",file.type),
                      g_densityplot,
                      dpi = dpi,
                      width = plot.width,
                      height = plot.height)
      return(g_densityplot)
      }else{
        return(g_densityplot)
        }
  }else{
    s_densplot <- ggplot2::ggplot(dplot_data,
                                  ggplot2::aes(x = value)) +
      ggplot2::geom_density(ggplot2::aes(fill = factor(stage,
                                                     levels=c(
                                                       "Before imputation",
                                                       "After imputation"))),
                            alpha =   0.25,
                            lwd = 0.1 )+
      ggplot2::xlab("Intensity") + ggplot2::ylab("Density")+
      ggplot2::scale_fill_manual(values = c(orig.col, imp.col))+
      ggplot2::theme_classic()+
      ggplot2::theme(panel.border = element_rect(fill = NA,
                                                 colour = "grey",
                                                 size = 0.5),
                     legend.title = element_blank(),
                     axis.ticks = element_line(colour = "grey"),
                     axis.title.x = element_blank(),
                     axis.title.y= element_blank(),
                     axis.text = element_text(size = text.size/2),
                     axis.line = element_line(colour = "grey",
                                              size = 0.5),
                     legend.text = element_text(size = text.size),
                     strip.text = element_text(size = text.size),
                     strip.background = element_rect(fill = "grey95",
                                                     colour = "grey",
                                                     size = 0.5))+
      ggplot2::facet_wrap(~sample,
                          nrow = nrow,
                          ncol = ncol,
                          strip.position = "top")
    if(save == TRUE){
      ggplot2::ggsave(paste0(file.name,".", file.type),
                      s_densplot,
                      dpi = dpi,
                      width = plot.width,
                      height = plot.height)
      return(s_densplot)
      }else{
        return(s_densplot)
      }
  }
  }
