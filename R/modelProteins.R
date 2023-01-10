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
#' @param fit_df A \code{fit_df} object from performing \code{find_dep}.
#' @param norm_df The \code{norm_df} object from which the \code{fit_df} object
#' was obtained.
#' @param sig Criteria to denote significance in differential expression.
#' Choices are \code{"adjP"} (default) for adjusted p-value or \code{"P"}
#' for p-value.
#' @param sig_cutoff Cutoff value for p-values and adjusted p-values in
#' differential expression. Default is \code{0.05}.
#' @param fc Minimum absolute log-fold change to use as threshold for
#' differential expression. Default is \code{1}.
#' @param n_top The number of top hits from \code{find_dep} to be used in
#' modeling. Default is \code{20}.
#' @param find_highcorr Logical. If \code{TRUE} (default), finds highly
#' correlated proteins.
#' @param corr_cutoff A numeric value specifying the correlation cutoff.
#' Default is \code{0.90}.
#' @param save_corrmatrix Logical. If \code{TRUE}, saves a copy of the
#' protein correlation matrix in a tab-delimited text file labeled
#' "Protein_correlation.txt" in the directory specified by \code{file_path}.
#' @param file_path A string containing the directory path to save the file.
#' @param rem_highcorr Logical. If \code{TRUE} (default), removes highly
#' correlated proteins (predictors or features).
#'
#' @details This function creates a data frame that contains protein intensities
#' for a user-specified number of top differentially expressed proteins.
#'  \itemize{\item Using \code{find_highcorr = TRUE}, highly correlated
#'  proteins can be identified, and can be removed with
#'  \code{rem_highcorr = TRUE}.
#'  \item Note: Most models will benefit from reducing correlation between
#'  proteins (predictors or features), therefore we recommend removing those
#'  proteins at this stage to reduce pairwise-correlation.
#'  \item If no or few proteins meet the significance threshold for differential
#'  expression, you may adjust \code{sig}, \code{fc}, and/or \code{sig_cutoff}
#'  accordingly to obtain more proteins for modeling.}
#'
#' @return A \code{model_df} object, which is a data frame of protein
#' intensities with proteins indicated by columns.
#' @seealso
#' \itemize{
#' \item \code{find_dep}, \code{normalize_data}
#' \item \code{\link[caret:findCorrelation]{caret: findCorrelation}}}
#'
#' @examples
#'
#' ## Create a model_df object with default settings.
#' covid_model_df1 <- pre_process(fit_df = covid_fit_df, norm_df = covid_norm_df)
#'
#' ## Change the correlation cutoff.
#' covid_model_df2 <- pre_process(covid_fit_df, covid_norm_df, corr_cutoff = 0.95)
#'
#' ## Change the significance criteria to include more proteins
#' covid_model_df3 <- pre_process(covid_fit_df, covid_norm_df, sig = "P")
#'
#' ## Change the number of top differentially expressed proteins to include
#' covid_model_df4 <- pre_process(covid_fit_df, covid_norm_df, sig = "P", n_top = 24)
#'
#' @export
pre_process <- function(fit_df,
                        norm_df,
                        sig = "adjP",
                        sig_cutoff = 0.05,
                        fc = 1,
                        n_top = 20,
                        find_highcorr = TRUE,
                        corr_cutoff = 0.90,
                        save_corrmatrix = FALSE,
                        file_path = NULL,
                        rem_highcorr = TRUE) {

  #binding for global variable
  logFC <- P.Value <- adj.P.Val <-  NULL

  # Extract the results from the differential expression analysis.
  exp_de <- limma::topTable(fit_df,
    coef = colnames(fit_df)[2],
    n = length(fit_df$df.total),
    adjust.method = "BH"
  )

  # Subset results by logFC and p-value cutoff
  if (sig == "P") {
    top_proteins <- rownames(subset(exp_de,
      abs(logFC) > fc & P.Value < sig_cutoff,
      drop = FALSE
    ))

    # Or default: based on adj.P value
  } else {
    top_proteins <- rownames(subset(exp_de,
      abs(logFC) > fc & adj.P.Val < sig_cutoff,
      drop = FALSE
    ))
  }

  # If the total number of DE proteins < n_top, replace n_top with that number.
  if (length(top_proteins) < n_top) {
    n_top <- length(top_proteins)
    message(paste0(
      "Total number of differentially expressed proteins (",
      n_top, ") ", "is less than n_top."
    ))
  }

  # Extract the top n_top hits from the top hit list
  top_proteins <- top_proteins[1:n_top]

  # Check if there are sig. proteins before moving on to pre-processing
  if (identical(top_proteins, character(0))) {
    stop(message
    (paste0(
        "No significant proteins found at ",
        sig,
        " < ",
        sig_cutoff,
        "."
    )))
  } else {
    # Extract intensity values for top proteins based on fc and sig_cutoff
    top_intensity <- subset(norm_df,
      rownames(norm_df) %in% top_proteins,
      drop = FALSE
    )
  }
  # Extract group or condition information from sample names in the data frame
  group <- factor(c(sapply(
    strsplit(colnames(top_intensity), "_"),
    getElement, 1
  )))

  # Transpose the data frame. Columns are now proteins and rows are samples.
  topint_trans <- as.data.frame(t(top_intensity))

  # Remove sample names.
  rownames(topint_trans) <- NULL

  # Add a new column with the group or condition information.
  # condition column is now the rightmost column in the data frame.
  topint_trans$condition <- group

  # For correlation calculations, make a matrix without the condition column
  topint_cor <- topint_trans[, seq_len(ncol(topint_trans)) - 1]

  # Create a correlation matrix
  cor_matrix <- cor(topint_cor)

  #Set temporary file_path if not specified
  if(is.null(file_path)){
    file_path <- tempdir()
  }

  if (save_corrmatrix == TRUE) {
    write.table(cor_matrix,
      file = paste0(file_path,"/Protein_correlation.txt"),
      row.names = TRUE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
  }

  if (find_highcorr == TRUE) {
    # Identify protein columns with high pairwise-correlation to remove
    highcor <- findCorrelation(cor_matrix, cutoff = corr_cutoff, names = TRUE)
    if (length(highcor != 0)) {
      message(
        "Following protein(s) show high pariwise-correlation"
      )
    } else {
      message(
        "None of the proteins show high pair-wise correlation."
      )
    }
    message(paste0(highcor, collapse = "\n"))

    if (rem_highcorr == TRUE) {
      topint_trans_1 <- topint_trans[, !(colnames(topint_trans) %in% highcor)]
      if (ncol(topint_trans_1) == ncol(topint_trans)) {
        message("No highly correlated proteins to be removed.")
      } else {
        message("Proteins with high pairwise-correlation have been removed.")
      }
    } else {
      warning("Proteins with high pairwise-correlation have NOT been removed.",
        call. = FALSE
      )
      topint_trans_1 <- topint_trans
    }
  } else {
    warning("Your data could have proteins with high pairwise-correlation.",
      call. = FALSE
    )
    topint_trans_1 <- topint_trans
  }

  # Convert condition names to R compatible names
  topint_trans_1$condition <- make.names(topint_trans_1$condition)

  # Convert condition to a factor (important for varimp calculations)
  topint_trans_1$condition <- factor(topint_trans_1$condition)
  return(topint_trans_1)
}


# Remove user-specified features -------------------------------------------
#' Remove user-specified proteins (features) from a data frame
#' @author Chathurani Ranathunge
#' @description This function removes user-specified proteins from a `model_df`
#' object
#'
#' @param model_df A \code{model_df} object.
#' @param rem_protein Name of the protein to remove.
#'
#' @details \itemize{\item After visualizing protein intensity variation
#' among conditions with \code{feature_plot} or after assessing the importance
#' of each protein in models using \code{varimp_plot}, you can choose to remove
#' specific proteins (features) from the data frame. \item For example, you can
#' choose to remove a protein from the \code{model_df} object if the protein
#' does not show distinct patterns of variation among conditions. This protein
#' may show mostly overlapping distributions in the feature plots.
#' \item Another incidence would be removing a protein that is very low in
#' variable importance in the models built using \code{train_models}. You can
#' visualize variable importance using \code{varimp_plot}.
#' }
#' @return A \code{model_df} object.
#'
#' @seealso \code{\link{feature_plot}}, \code{\link{pre_process}}
#'
#' @examples
#'
## Create a model_df object with default settings.
#' covid_model_df <- pre_process(fit_df = covid_fit_df, norm_df = covid_norm_df)
#'
#' ## Remove sp|P22352|GPX3_HUMAN protein from the model_df object
#' covid_model_df1 <- rem_feature(covid_model_df, rem_protein = "sp|P22352|GPX3_HUMAN")
#'
#' @export
rem_feature <- function(model_df,
                        rem_protein) {
  df_rem <- model_df[, -grep(rem_protein, colnames(model_df), fixed = TRUE)]
  message(paste0("Protein ", rem_protein, " has been removed."))
  return(df_rem)
}

# Split data frame -------------------------------------------------------------
#' Split the data frame to create training and test data
#' @description This function can be used to create balanced splits of the
#' protein intensity data in a `model_df` object to create training and test data
#'
#' @author Chathurani Ranathunge
#'
#' @import caret
#'
#' @param model_df A \code{model_df} object from performing \code{pre_process}.
#' @param train_size The size of the training data set as a proportion of the
#' complete data set. Default is 0.8.
#' @param seed Numerical. Random number seed. Default is \code{NULL}
#'
#' @details This function splits the \code{model_df} object in to training and
#' test data sets using random sampling while preserving the original
#' class distribution of the data. Make sure to fix the random number seed with
#' \code{seed} for reproducibility
#'
#' @return A list of data frames.
#' @seealso
#' \itemize{
#' \item \code{pre_process}
#' \item \code{\link[caret:createDataPartition]{createDataPartition}}}
#'
#' @examples
#'
#' ## Create a model_df object
#' covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#'
#' ## Split the data frame into training and test data sets using default settings
#' covid_split_df1 <- split_data(covid_model_df, seed = 8314)
#'
#' ## Split the data frame into training and test data sets with 70% of the
#' ## data in training and 30% in test data sets
#' covid_split_df2 <- split_data(covid_model_df, train_size = 0.7, seed = 8314)
#'
#' ## Access training data set
#' covid_split_df1$training
#'
#' ## Access test data set
#' covid_split_df1$test
#'
#' @export
split_data <- function(model_df,
                       train_size = 0.80,
                       seed = NULL) {
  set.seed(seed)
  train_index <- createDataPartition(model_df$condition,
    p = train_size,
    list = FALSE
  )

  # Use the train_index to subset the data frame
  train_df <- model_df[train_index, ]
  test_df <- model_df[-train_index, ]

  # Remove rownames
  rownames(train_df) <- NULL
  rownames(test_df) <- NULL

  # Create a list with test and training data frames
  split_dataframes <- list()
  split_dataframes[[1]] <- train_df
  split_dataframes[[2]] <- test_df

  # Rename the items of the list
  names(split_dataframes) <- c("training", "test")


  return(split_dataframes)
}


# Train models -------------------------------------------------------------
#' Train machine learning models on training data
#' @description This function can be used to train models on protein intensity
#' data using different machine learning algorithms
#'
#' @author Chathurani Ranathunge
#'
#' @import caret
#' @importFrom kernlab sigest
#' @import xgboost
#' @import naivebayes
#'
#' @param split_df A \code{split_df} object from performing \code{split_data}.
#' @param resample_method The resampling method to use. Default is
#' \code{"repeatedcv"} for repeated cross validation.
#' See \code{\link[caret:trainControl]{trainControl}} for
#' details on other available methods.
#' @param resample_iterations Number of resampling iterations. Default is
#' \code{10}.
#' @param num_repeats The number of complete sets of folds to compute (For
#' \code{resampling method = "repeatedcv"} only).
#' @param algorithm_list A list of classification or regression algorithms to
#' use.
#' A full list of machine learning algorithms available through
#' the \code{caret} package can be found here:
#' \url{http://topepo.github.io/caret/train-models-by-tag.html}. See below for
#' default options.
#' @param seed Numerical. Random number seed. Default is \code{NULL}
#' @param ... Additional arguments to be passed on to
#' \code{\link[caret:train]{train}} function in the \code{caret} package.
#'
#' @details \itemize{\item \code{train_models} function can be used to first
#' define the control parameters to be used in training models, calculate
#' resampling-based performance measures for models based on a given set of
#' machine-learning algorithms, and output the best model for each algorithm.
#' \item In the event that \code{algorithm_list} is not provided, a default
#' list of four classification-based machine-learning algorithms will be used
#' for building and training models. Default \code{algorithm_list}:
#'  "svmRadial", "rf", "glm", "xgbLinear."
#' \item Note: Models that fail to build are removed from the output.
#' \item Make sure to fix the random number seed with
#' \code{seed} for reproducibility}
#'
#' @return
#' A list of class \code{train} for each machine-learning algorithm.
#' See \code{\link[caret:train]{train}} for more information on accessing
#' different elements of this list.
#'
#' @seealso
#' \itemize{
#' \item \code{pre_process}
#' \item \code{\link[caret:trainControl]{trainControl}}
#' \item \code{\link[caret:train]{train}}
#' }
#'
#' @references Kuhn, Max. "Building predictive models in R using the caret
#' package." Journal of statistical software 28 (2008): 1-26.
#'
#' @examples
#' \donttest{
#'
#' ## Create a model_df object
#' covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#'
#' ## Split the data frame into training and test data sets
#' covid_split_df <- split_data(covid_model_df, seed = 8314)
#'
#' ## Fit models based on the default list of machine learning (ML) algorithms
#' covid_model_list1 <- train_models(split_df = covid_split_df, seed = 351)
#'
#' ## Fit models using a user-specified list of ML algorithms.
#' covid_model_list2 <- train_models(
#'   covid_split_df,
#'   algorithm_list = c("naive_bayes", "glmboost"),
#'   seed = 351
#' )
#'
#' ## Change resampling method and resampling iterations.
#' covid_model_list3 <- train_models(
#'   covid_split_df,
#'   resample_method = "cv",
#'   resample_iterations = 50,
#'   seed = 351
#' )
#' }
#'
#' @export
train_models <- function(split_df,
                         resample_method = "repeatedcv",
                         resample_iterations = 10,
                         num_repeats = 3,
                         algorithm_list,
                         seed = NULL,
                         ...) {

  # If algorithm_list is not provided, use the default list of algorithms.
  if (missing(algorithm_list)) {
    algorithm_list <- c("svmRadial", "rf", "glm", "xgbLinear")
  }

  # Set trainControl parameters for resampling
  set.seed(seed)
  fit_control <- trainControl(
    method = resample_method,
    number = resample_iterations,
    repeats = num_repeats,
    classProbs = TRUE
  )

  # Extract the training data set from the split_df object
  training_data <- split_df$training


  # Train models using ML algorithms from the algorithm_list.
  model_list <- lapply(
    setNames(algorithm_list, algorithm_list),
    function(x) {
      tryCatch(
        {
          set.seed(seed)
          message(paste0("\n", "Running ", x, "...", "\n"))
          train(condition ~ .,
            data = training_data,
            trControl = fit_control,
            method = x,
            ...
          )
        },
        error = function(e) {
          message(paste0(
            x,
            " failed."
          ))
        }
      )
    }
  )

  message(paste0("Done!"))

  # Drop models that failed to build from the list
  model_list <- Filter(Negate(is.null), model_list)
  return(model_list)
}

# Test the models -------------------------------------------------------------
#' Test machine learning models on test data
#' @description This function can be used to predict test data using models
#' generated by different machine learning algorithms
#'
#' @author Chathurani Ranathunge
#'
#' @import stats
#' @import caret
#' @importFrom reshape2 melt
#'
#' @param model_list A \code{model_list} object from performing
#' \code{train_models}.
#' @param split_df A \code{split_df} object from performing \code{split_data}.
#' @param type Type of output. Set \code{type} as "prob" (default) to output
#' class probabilities, and "raw" to output class predictions.
#' @param save_confusionmatrix Logical. If \code{TRUE}, a tab-delimited
#' text file ("Confusion_matrices.txt") with confusion matrices in the
#' long-form data format will be saved in the directory specified by
#' \code{file_path}.
#' See below for more details.
#' @param file_path A string containing the directory path to save the file.
#' @param ... Additional arguments to be passed on to
#' \code{\link[stats:predict]{predict}}.
#'
#' @details \itemize{\item \code{test_models} function uses
#' models obtained from \code{train_models} to predict a given test data set.
#' \item Setting \code{type = "raw"} is required to obtain confusion matrices.
#' \item Setting \code{type = "prob"} (default) will output a list of
#' probabilities that can be used to generate ROC curves using \code{roc_plot}.}
#'
#' @return
#' \itemize{\item \code{probability_list}: If \code{type = "prob"}, a list of
#' data frames containing class probabilities for each method in the
#' \code{model_list} will be returned.
#' \item \code{prediction_list}: If \code{type = "raw"}, a list of factors
#' containing class predictions for each method will be returned.}
#'
#' @seealso
#' \itemize{
#' \item \code{split_df}
#' \item \code{train_models}
#' \item \code{\link[stats:predict]{predict}}
#' \item \code{\link[caret:confusionMatrix]{confusionMatrix}}
#' }
#'
#' @examples
#' \donttest{
#' ## Create a model_df object
#' covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#'
#' ## Split the data frame into training and test data sets
#' covid_split_df <- split_data(covid_model_df)
#'
#' ## Fit models using the default list of machine learning (ML) algorithms
#' covid_model_list <- train_models(covid_split_df)
#'
#' # Test a list of models on a test data set and output class probabilities,
#' covid_prob_list <- test_models(model_list = covid_model_list, split_df = covid_split_df)
#' }
#'
#' \dontrun{
#' # Save confusion matrices in the working directory and output class predictions
#' covid_pred_list <- test_models(
#'   model_list = covid_model_list,
#'   split_df = covid_split_df,
#'   type = "raw",
#'   save_confusionmatrix = TRUE,
#'   file_path = "."
#' )
#' }
#'
#' @export
test_models <- function(model_list,
                        split_df,
                        type = "prob",
                        save_confusionmatrix = FALSE,
                        file_path = NULL,
                        ...) {

  # Extract test data from the split_df object
  test_data <- split_df$test

  # Predict test data
  pred_list <- lapply(
    model_list,
    function(x) {
      message(paste0(
        "\n",
        "Testing ",
        x$method,
        "...",
        "\n"
      ))
      predict(x,
        test_data,
        type = type
      )
    }
  )
  message(paste0("\n", "Done!"))

  # Get confusion matrices and associated statistics
  if (type == "raw") {
    cm_list <- lapply(
      pred_list,
      function(x) {
        confusionMatrix(
          x,
          test_data$condition
        )
      }
    )

    #Set temporary file_path if not specified
    if(is.null(file_path)){
      file_path <- tempdir()
    }

    if (save_confusionmatrix == TRUE) {

      # Convert c.matrices to long-form data frames
      cm_df <- lapply(
        cm_list,
        function(x) reshape2::melt(as.table(x))
      )
      # Get the list of methods
      method_list <- names(cm_df)

      # Add method names to the data frames
      cm_dfm <- lapply(
        seq_along(method_list),
        function(x) {
          cm_df[[x]]["method"] <- method_list[x]
          cm_df[[x]]
        }
      )

      # Combine all data frames into one
      cm_dfm_long <- do.call(
        "rbind",
        cm_dfm
      )

      # Add column names before saving
      colnames(cm_dfm_long) <- c("Prediction", "Reference", "Value", "Method")

      # Save data in a text file
      write.table(cm_dfm_long,
        file = paste0(file_path, "/Confusion_matrices.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  return(pred_list)
}
