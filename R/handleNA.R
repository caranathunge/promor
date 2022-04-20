# Visualize missing data -------------------------------------
#' Visualize missing data
#' @description This function visualizes the patterns of missing value
#' occurrence using a heatmap.
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @param df A \code{raw.df} object (output from \code{create.df}).
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


  #Convert the data into long format for plotting and make necessary changes
  #i.e. adding column headers etc.
  df <- as.matrix(df)
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
                            ggplot2::aes(x = reorder(sample,
                                                     value,
                                                     na.rm = TRUE,
                                                     FUN = x.FUN),
                                         y = reorder(protgroup,
                                                     value,
                                                     na.rm = TRUE,
                                                     FUN = y.FUN),
                                         fill = value))

    }else if(reorder.x==FALSE & reorder.y == TRUE){
    hmap <- ggplot2::ggplot(hmap_df,
                            ggplot2::aes(x = sample,
                                         y = reorder(protgroup,
                                                     value,
                                                     na.rm = TRUE,
                                                     FUN = y.FUN),
                                         fill = value))

    }else if(reorder.x== TRUE & reorder.y == FALSE){
    hmap <- ggplot2::ggplot(hmap_df,
                            ggplot2::aes(x = reorder(sample,
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
                   axis.text.x = element_text(angle = 90,
                                              size = text.size),
                   axis.text.y=element_text(size = text.size),
                   axis.ticks.y = element_blank(),
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
  }else{
    return(hmap)
  }
  }

# Impute missing data -----------------------------------------------------
#' Impute missing data
#' @author Chathurani Ranathunge
#' @description This function imputes missing data using a user-specified
#' imputation method.
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca
#' @importFrom VIM kNN
#' @export



impute.NA <- function(x,
                      method = "minProb",
                      tune.sigma = 1,
                      q = 0.001,
                      maxiter = 6,
                      ntree = 20,
                      verbose= TRUE,
                      nPcs = 2,
                      center = TRUE,
                      imp_var = FALSE){


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
    print(sd.temp)
    for (i in 1:(nSamples)) {
      dataSet.to.impute.temp = rnorm(nFeatures,
                                     mean = min.samples[i],
                                     sd = sd.temp)
      dataSet.imputed[which(is.na(dataSet.mvs[, i])), i] =
        dataSet.to.impute.temp[which(is.na(dataSet.mvs[,i]))]
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
      dataSet.imputed[which(is.na(dataSet.mvs[, i])), i] = lowQuantile.samples[i]
    }
    return(dataSet.imputed)
  }

#Run the user-specified imputation method

  if (method == "minDet"){
    df_imputed_mindet <- impute.MinDet(x,
                                       q = q)
    return(df_imputed_mindet)

    }else if(method == "RF"){
    df_imp_temp <- missForest::missForest(x,
                                          maxiter = maxiter,
                                          ntree = ntree,
                                          verbose= verbose)
    df_imputed_RF <- df_imp_temp$ximp
    return(df_imputed_RF)

    }else if (method == "kNN"){
    df_imputed_knn <- VIM::kNN(x, imp_var = imp_var)
    rownames(df_imputed_knn) <- rownames(x)
    return(df_imputed_knn)

    }else if (method == "SVD") {
    x[is.nan(x)] <- NA
    df_imp_temp <- pcaMethods::pca(object = x,
                                   method = "svdImpute",
                                   nPcs = nPcs,
                                   center = center,
                                   verbose= verbose)
    df_imputed_svd <- completeObs(df_imp_temp)
    return(df_imputed_svd)

    }else if (method == "minProb"){
    df_imputed_minprob <- impute.Min.Prob(x,
                                         q = q,
                                         tune.sigma = tune.sigma)
    return(df_imputed_minprob)
    }
  }


# Visualize the imputation effects ----------------------------------------
#' Visualize the imputation effects
#' @author Chathurani Ranathunge
#' @description This function compares imputed data to original data with a
#' user-defined plot.
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export

#Options: Global density plot, Sample-wise density plot
impute.densplot <- function(original,
                            imputed,
                            global = TRUE,
                            alpha=0.25,
                            lwd = 0.1,
                            xlabel = "Intensity",
                            ylabel="Density",
                            xlab.size = 7,
                            ylab.size =7,
                            col1 ="blue",
                            col2="red",
                            strip.txt.size = 5,
                            strip.ln.col = "grey",
                            strip.ln.size = 0.5,
                            nrow = 15, ncol = 5,
                            dpi = 80,
                            filename= "Impute_densityplot",
                            filetype="pdf",
                            width = 7,
                            height= 7,
                            save = TRUE){

  #Make necessary changes and combine data frames before plotting
  orig_data <- reshape2::melt(original, na.rm = FALSE)
  orig_data$mjprot <- sapply(strsplit(as.character(orig_data[,1]),';'),
                              getElement,1)
  orig_data[1] <- NULL
  orig_df1 <- as.data.frame(orig_data)
  colnames(orig_df1) <- c( "sample", "value", "protgroup")
  orig_df1$stage <- "Before imputation"

  imp_data <- reshape2::melt(imputed, na.rm = FALSE)
  imp_data$mjprot <- sapply(strsplit(as.character(imp_data[,1]),';'),
                              getElement,1)
  imp_data[1] <- NULL
  imp_df1 <- as.data.frame(imp_data)
  colnames(imp_df1) <- c( "sample", "value", "protgroup")
  imp_df1$stage <- "After imputation"
  dplot_data <- rbind(orig_df1, imp_df1)

  #Plot global density plot
  if(global ==TRUE){
    g_densityplot <- ggplot2::ggplot(dplot_data,
                                     ggplot2::aes(x=value)) +
      ggplot2::geom_density(ggplot2::aes(fill = factor(stage,
                                                       levels = c(
                                                         "Before imputation",
                                                         "After imputation"))),
                            alpha=alpha,
                            lwd = lwd )+
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel)+
      ggplot2::scale_fill_manual(values=c(col1,col2))+
      ggplot2::theme_classic()+
      ggplot2::theme(legend.title = element_blank(),
                     axis.title.x = element_text(size = xlab.size),
                     axis.title.y = element_text(size = ylab.size))


    if(save == TRUE){
      ggplot2::ggsave(paste0(filename,".",filetype),
                      g_densityplot,
                      dpi = dpi,
                      width = width,
                      height = height)
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
                            alpha = alpha,
                            lwd = lwd )+
      ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)+
      ggplot2::scale_fill_manual(values = c(col1, col2))+
      ggplot2::theme_classic()+
      ggplot2::theme(legend.title = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y= element_blank(),
                     strip.text = element_text(size = strip.txt.size),
                     strip.background = element_rect(colour = strip.ln.col,
                                                     size = strip.ln.size))+
      ggplot2::facet_wrap(~sample,
                          nrow = nrow,
                          ncol = ncol,
                          strip.position = "top")
    if(save == TRUE){
      ggplot2::ggsave(paste0(filename,".", filetype),
                      s_densplot,
                      dpi = dpi,
                      width = width,
                      height = height)
      }else{
        return(s_densplot)
      }
  }
  }
