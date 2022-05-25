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
#' @param save.corrmatrix Logical. If \code{TRUE}, saves a copy of the
#' protein correlation matrix in a tab-delimited text file labeled
#' "Protein_correlation.txt" in the working directory.
#' @param rem.highcorr Logical. If \code{TRUE} (default), removes highly
#' correlated proteins (predictors or features).
#'
#' @details This function creates a data frame that contains protein intensities
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
                        save.corrmatrix = FALSE,
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

  if (save.corrmatrix == TRUE){
    write.table(cor_matrix,
                file="./Protein_correlation.txt",
                row.names = TRUE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE)
  }

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


# Remove user-specified predictors -------------------------------------------
#' Remove user-specified proteins (predictors) from the data frame
#' @author Chathurani Ranathunge
#' @description This function removes user-specified proteins from the
#' data frame.
#'
#' @param model.df A \code{model.df} object.
#' @param rem.protein Name of the protein to remove.
#'
#' @details \itemize{\item After visualizing protein intensity variation
#' among conditions with \code{predictor_plot}, you can choose to remove
#' specific proteins (predictors) from the data frame before modeling if they
#' do not show distinct patterns of variation among conditions. For
#' example, these proteins could have overlapping density distributions.
#' }
#' @return A \code{model.df} object.
#'
#' @seealso \code{\link{predictor_plot}}, \code{\link{pre_process}}
#'
#' @examples
#' \dontrun{
#'
#' ## Remove protein "A113403" from model_df data frame.
#' model_df <- rem_predictor(model_df, "A113403")
#'
#' }
#'
#' @export
rem_predictor <- function(model.df, rem.protein){
  df_rem <- df[, -grep(rem.protein, colnames(df))]
  message(paste0("Protein ", rem.protein, "has been removed."))
  return(df_rem)

}

# Split data frame -------------------------------------------------------------
#' Split the data frame to create training and test data
#' @description This function can be used to create balanced splits of the
#' protein intensity data to create training and test data
#'
#' @author Chathurani Ranathunge
#'
#' @import caret
#'
#' @param model.df A \code{model.df} object from performing \code{pre_process}.
#' @param train.size The size of the training data set as a percentage of the
#' complete data set. Default is 0.8.
#'
#' @details This function splits the \code{model.df} object in to training and
#' test data sets using random sampling while preserving the original
#' class distribution of the data.
#'
#' @return A list of data frames.
#' @seealso
#' \itemize{
#' \item \code{pre_process}
#' \item \code{\link[caret:createDataPartition]{createDataPartition}}}
#' @examples
#' \dontrun{
#'
#'  ## Split the data frame into training and test data sets with 70% of the
#'  data in training and 30% in test data sets
#'  split_data <- split_df(model-df, train.size = 0.7)
#'
#'  ## Access training data set
#'  split_data$training
#'
#'  ## Access test data set
#'  split_data$test
#'
#'  }
#'
#'
#' @export
split_df <- function(model.df,
                     train.size = 0.80
                     ){
  set.seed(8314)
  train_index <- createDataPartition(model.df$Condition,
                                     p = train.size,
                                     list = FALSE)

  #Use the train_index to subset the data frame
  train_df <- model.df[train_index,]
  test_df <- model.df[-train_index,]

  #Remove rownames
  rownames(train_df)  <- NULL
  rownames(test_df)  <- NULL

  #Create a list with test and training data frames
  split_dataframes <- list()
  split_dataframes[[1]] <- train_df
  split_dataframes[[2]] <- test_df

  #Rename the items of the list
  names(split_dataframes) <- c("training", "test")

  return(split_dataframes)

}


# Train the model -------------------------------------------------------------
#' Train machine learning models on protein intensity data
#' @description This function can be used to train models on protein intensity
#' data using a number of different machine learning algorithms
#'
#' @author Chathurani Ranathunge
#'
#' @import caret
#'
#' @param split.df A \code{split.df} object from performing \code{split_df}.
#' @param resample.method The resampling method to use. Default is "repeatedcv"
#' for repeated cross validation.
#' See \code{\link[caret:trainControl]{trainControl}} for
#' details on other available methods.
#' @param resample.iterations Number of resampling iterations. Default is 10.
#' @param num.repeats The number of complete sets of folds to compute (For
#' \code{resampling method = "repeatedcv"} only).
#' @param algorithm.list A list of classification or regression algorithms to use.
#' A full list of machine learning algorithms available through
#' the \code{\link[caret]{caret package}} can be found here:
#' \url{http://topepo.github.io/caret/train-models-by-tag.html}. see below for
#' default options.
#'
#' @details \itemize{\item \code{train_model} function can be used to first
#' define the control parameters to be used in training models, calculate
#' resampling-based performance measures for models based on a given set of
#' machine-learning algorithms, and output the best model for each algorithm.
#' \item In the event that \code{algorithm.list} is not provided, a default
#' list of five classification-based machine-learning algorithms will be used
#' for building and training models. Default \code{algorithm.list}:
#' "rf", "ranger", "earth", "knn", "svmLinear."
#' \item Note: Models that failed to build are removed from the output. }
#'
#' @return
#' A list of class \code{train} for each machine-learning algorithm.
#' See \code{\link[caret:train]{train}} for more information on accessing
#' different elements of this list.
#'
#' @seealso
#' \itemize{
#' \item \code{pre_process}
#' \item \code{\link[caret:trainControl]{caret: trainControl}}
#' \item \code{\link[caret:train]{caret: train}}
#' }
#' @examples
#' \dontrun{
#' #Fit models based on the default list of ML algorithms.
#' model_fit <- train_model(split.df = split_df)
#'
#' #Fit models using a user-specified list of ML algorithms.
#' model_fit <- train_model(split.df = split_df,
#'                          algorithm.list = c("svmRadial", "rpart"))
#'
#' #Change resampling method and resampling iterations.
#' model_fit <- train_model(split.df = split_df,
#'                          resample.method = "cv",
#'                          resample.iterations = 50)
#'  }
#'
#'
#' @export
train_model <- function(split.df,
                        resample.method = "repeatedcv",
                        resample.iterations = 10,
                        num.repeats = 3,
                        algorithm.list,
                        ...
                        ){

  #If algorithm.list is not provided, use the default list of algorithms.
  if(missing(algorithm.list)) {
    algorithm.list = c("rf",  "earth", "knn", "svmLinear", "ranger")
  }

  #Set trainControl parameters for resampling
  set.seed(351)
  fit_control <- trainControl(method = resample.method,
                              number = resample.iterations,
                              repeats = num.repeats)

  #Extract the training data set from the split.df object
  training_data <- split.df$training


  #Train models using ML algorithms from the algorithm.list.
  model_list <- lapply(setNames(algorithm.list, algorithm.list),
                       function (x) tryCatch({set.seed(351);
                         message(paste0("\n", "Running ", x, "...", "\n"));
                         train(Condition ~ .,
                               data = training_data,
                               trControl = fit_control,
                               method = x,
                               importance = TRUE,
                               ...)},
                         error = function(e) {message(paste0(x,
                                                           " failed."))}))

  message(paste0("Done!"))

  #Drop models that failed to build from the list
  model_list <- Filter(Negate(is.null), model_list)
  return(model_list)
}
