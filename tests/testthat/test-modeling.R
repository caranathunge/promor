#Testing modeling functions
#Load test norm_df object (subset of 50 proteins from original ecoli data set)
test_norm_df <- read.csv( "./testdata/test_ecoli_norm.txt",
                          sep = "\t",
                          stringsAsFactors = TRUE)
test_norm_df <- as.matrix(test_norm_df)


#Load test fit_df object.
test_fit_df <- readRDS(file = "./testdata/test_ecoli_fit.RDS")

#Load model_df
test_model_df <- readRDS( "./testdata/test_ecoli_model.RDS")

#Testing pre_process------------------------------------------------------------
#Test if the function works - default
testthat::test_that("pre_process default works", {
  suppressMessages(
    model_df <- pre_process(test_fit_df,
                               test_norm_df))
    testthat::expect_equal(model_df, test_model_df)

})

#Load model_df - highcorr proteins included
test_model_df_hc <- readRDS( "./testdata/test_ecoli_model_hc.RDS")

#Test if function works when rem_highcorr = FALSE
testthat::test_that("pre_process works when rem_highcorr = FALSE", {

  suppressWarnings(
    model_df <- pre_process(test_fit_df,
                           test_norm_df,
                           rem_highcorr = FALSE))
  testthat::expect_equal(model_df, test_model_df_hc)

})

# Testing rem_feature function--------------------------------------------------
test_remf <- readRDS(file = "./testdata/test_ecoli_remf.RDS")

#Test if function works
testthat::test_that("rem_feature works", {
  suppressMessages(
    rem_p <- rem_feature(test_model_df,
                               "A6H8Z3;Q15042;C9J837"))
    testthat::expect_equal(rem_p, test_remf)

})

# Testing split_data------------------------------------------------------------
#load split_df object
test_splitdf <- readRDS( "./testdata/test_ecoli_splitdf.RDS")

#Test if function works
testthat::test_that("split_data works", {
  suppressMessages(testthat::expect_equal(split_data(test_model_df,
                           train_size = 0.50), test_splitdf))

})

# Testing train_models----------------------------------------------------------
#load model_list object
test_model_list <- readRDS( "./testdata/test_ecoli_modellist.RDS")

#rf - confusion_matrix dimensions
rf_conf_mat <- length(test_model_list$rf$finalModel$confusion)

#glm - coefficients
glm_coef <- length(test_model_list$glm$finalModel$coefficients)

#Test if function works
testthat::test_that("train_models works", {
  suppressWarnings(suppressMessages(model_list <- train_models(test_splitdf,
                                             algorithm_list = c("rf", "glm"))))

  rf_len <- length(model_list$rf$finalModel$confusion)
  glm_len <- length(model_list$glm$finalModel$coefficients)

  testthat::expect_equal(rf_len, rf_conf_mat)
  testthat::expect_equal(glm_len, glm_coef)

})

#Testing if test_models works---------------------------------------------------
#load prob_list object
test_prob_list <- readRDS( "./testdata/test_ecoli_problist.RDS")

#Test if function works
testthat::test_that("test_models works", {
  prob_list <- test_models(test_model_list,
                           test_splitdf)
  testthat::expect_equal(prob_list, test_prob_list)

})


