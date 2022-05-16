# Pre-process protein intensity data for modeling -----------------------------
#' Pre-process protein intensity data for modeling
#' @description This function pre-processes protein intensity data from
#' the top differentially expressed proteins identified with `find_dep` for
#' modeling.
#'
#' @author Chathurani Ranathunge
#'
#' @import caret
#' @importFrom stats cor
#'
#' @param fit.df A \code{fit.df} object from performing \code{find_dep}.
#' @param norm.df The \code{norm.df} object from which the \code{fit.df} object
#' was obtained.
#' @param sig Criteria to denote significance in differential expression.
#' Choices are \code{"adjP"} (default) for adjusted p-value or \code{"P"}
#' for p-value.
#' @param sig.cutoff Cutoff value for p-values and adjusted p-values in
#' differential expression. Default is 0.05.
#' @param FC Minimum absolute log-fold change to use as threshold for
#' differential expression. Default is 1.
#' @param n.top The number of top hits from \code{find_dep} to be used in
#' modeling. Default is 20.
#' @param find.highcorr Logical. If \code{TRUE} (default), finds highly
#' correlated proteins.
#' @param corr.cutoff A numeric value specifying the correlation cutoff.
#' Default is 0.90.
#' @param rem.highcorr Logical. If \code{TRUE} (default), removes highly correlated
#' proteins (predictors or features).
#'
#' @details This function creates a dataframe that contains protein intensities
#' for a user-specified number of top differentially expressed proteins.
#'  \itemize{\item Using \code{find.highcorr = TRUE}, highly correlated
#'  proteins can be identified, and can be removed with
#'  \code{rem.highcorr = TRUE}.
#'  \item Note: Most models will benefit from reducing correlation between
#'  proteins (predictors or features), therefore we recommend removing those
#'  proteins at this stage to reduce pairwise-correlation.
#'  \item If no proteins meet the significance threshold for differential
#'  expression, you may adjust \code{sig}, \code{FC}, and \code{sig.cutoff}
#'  accordingly to obtain a set of proteins for modeling.}
#'
#' @return A \code{model.df} object, which is a data frame of protein
#' intensities with proteins indicated by columns.
#' @seealso
#' \itemize{
#' \item \code{find_dep}, \code{normalize_data}
#' \item \code{\link[caret:findCorrelation]{findCorrelation}}}
#' @examples
#' \dontrun{
#'
#'  ## Create a model.df object with default settings.
#'  model_df <- pre_process(fit_df, norm_df)
#'
#'  ## Change the correlation cutoff.
#'  model_df <- pre_process(fit_df, norm_df, corr.cutoff = 0.95)
#'  }
#'
#'
#' @export

pre_process <- function(fit.df,
                        norm.df,
                        sig = "adjP",
                        sig.cutoff = 0.05,
                        FC = 1,
                        n.top = 20,
                        find.highcorr = TRUE,
                        corr.cutoff = 0.90,
                        rem.highcorr = TRUE
                        ){

  #Extract the results from the differential expression analysis.
  exp_DE <- limma::topTable(fit.df,
                            coef = colnames(fit.df)[2],
                            n = length(fit.df$df.total),
                            adjust.method = "BH")

  #Subset results by logFC and p-value cutoff
  if(sig == "P"){
    top_proteins <- rownames(subset(exp_DE,
                                    abs(logFC) > FC & P.Value < sig.cutoff,
                                    drop = FALSE))

  #Or default: based on adj.P value
  }else{
    top_proteins <- rownames(subset(exp_DE,
                                    abs(logFC) > FC & adj.P.Val < sig.cutoff,
                                    drop = FALSE))
  }

  #If the total number of DE proteins < n.top, replace n.top with that number.
  if (length(top_proteins) < n.top){
    n.top = length(top_proteins)
    message(paste0("Total number of differentially expressed proteins (",
                 n.top, ") ", "is less than n.top."))
  }

  #Extract the top n.top hits from the top hit list
  top_proteins <- top_proteins[1:n.top]

  #Check if there are sig. proteins before moving on to pre-processing
  if (identical(top_proteins, character(0))){
    stop(message
         (paste0("No significant proteins found at ",
                 sig,
                 " < ",
                 sig.cutoff,
                 "." )))
  }else{
    #Extract intensity values for top proteins based on FC and sig.cutoff
    top_intensity <- subset(norm.df,
                            rownames(norm.df) %in% top_proteins,
                            drop = FALSE)
  }
  #Extract group or condition information from sample names in the data frame
  group <- factor(c(sapply(strsplit(colnames(top_intensity), "_"),
                           getElement, 1)))

  #Transpose the data frame. Columns are now proteins and rows are samples.
  topint_trans <- as.data.frame(t(top_intensity))

  #Remove sample names.
  rownames(topint_trans) <- NULL

  #Add a new column with the group or condition information.
  #Condition column is now the rightmost column in the data frame.
  topint_trans$Condition <- group

  #For correlation calculations, make a matrix without the Condition column
  topint_cor <- topint_trans[, 1: ncol(topint_trans)-1]

  #Create a correlation matrix
  cor_matrix <- cor(topint_cor)


  if (find.highcorr == TRUE){
  #Identify protein columns with high pairwise-correlation to remove
  highcor <- findCorrelation(cor_matrix, cutoff = corr.cutoff,  names = TRUE)
  message(
    "Following proteins show high pariwise-correlation")
  message(paste0(highcor, collapse = "\n"))

  if (rem.highcorr == TRUE){
    topint_trans <- topint_trans[,!(colnames(topint_trans) %in% highcor)]

    message("Proteins with high pairwise-correlation have been removed.")

    }else{
      warning("Proteins with high pairwise-correlation have NOT been removed.",
              call. = FALSE)

      }
  }else{
    warning("Your data could have proteins with high pairwise-correlation.",
            call. = FALSE )
  }


  return(topint_trans)
}


