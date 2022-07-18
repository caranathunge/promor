# Testing modeling functions----------------------------------------------------
## Input files------------------------------------------------------------------
# create fake norm_df
set.seed(581)
x <- data.frame(
  Control_1 = rnbinom(100, size = 0.5, prob = 1e-5),
  Control_2 = rnbinom(100, size = 0.5, prob = 1e-5),
  Case_1 = rnbinom(100, size = 0.5, prob = 1e-5),
  Case_2 = rnbinom(100, size = 0.5, prob = 1e-5),
  row.names = paste0("protein_", 1:100)
)
x[x == 0] <- NA
fake_norm_df <- log(x, base = 2)

# Load fake fit_df
fake_fit_df <- readRDS("./testdata/fake_fit.RDS")

# Create test_model_df
test_model_df <- structure(list(protein_18 = c(
  17.316272691676, 16.9847967072082,
  8.47167521439204, 15.5617652636277
), protein_24 = c(
  3.8073549220576,
  12.4553272203046, 13.5992154950177, 17.5539742084579
), protein_78 = c(
  15.9073131982696,
  9.67771964164101, 16.6606784030219, 17.4781376069187
), condition = structure(c(
  2L,
  2L, 1L, 1L
), levels = c("Case", "Control"), class = "factor")), row.names = c(
  NA,
  -4L
), class = "data.frame")


# create test_model_df with highcorr proteins
test_model_df_hc <- structure(list(protein_18 = c(
  17.316272691676, 16.9847967072082,
  8.47167521439204, 15.5617652636277
), protein_24 = c(
  3.8073549220576,
  12.4553272203046, 13.5992154950177, 17.5539742084579
), protein_35 = c(
  13.5634346357102,
  10.0953970227926, 16.5253219843842, 16.8754613979164
), protein_61 = c(
  6.7279204545632,
  10.2070143201775, 14.5470752676215, 16.495355392236
), protein_71 = c(
  14.7317432907078,
  15.0738486783698, 18.2667728177542, 16.5871182396906
), protein_78 = c(
  15.9073131982696,
  9.67771964164101, 16.6606784030219, 17.4781376069187
), protein_83 = c(
  17.289496610897,
  17.8753954004763, 13.9438883698759, 16.3635530265969
), condition = structure(c(
  2L,
  2L, 1L, 1L
), levels = c("Case", "Control"), class = "factor")), row.names = c(
  NA,
  -4L
), class = "data.frame")

# Testing pre_process------------------------------------------------------------
# Test if the function works - default
test_that("pre_process default works", {
  suppressMessages(
    expect_equal(
      pre_process(fake_fit_df,
        fake_norm_df,
        sig = "P",
        sig_cutoff = 0.1
      ),
      test_model_df
    )
  )
})



## Test if function works when rem_highcorr = FALSE-----------------------------
test_that("pre_process with highcorr works", {
  suppressWarnings(suppressMessages(
    expect_equal(
      pre_process(fake_fit_df,
        fake_norm_df,
        sig = "P",
        sig_cutoff = 0.1,
        rem_highcorr = FALSE
      ),
      test_model_df_hc
    )
  ))
})

# Testing rem_feature------------------------------------------------------------
# Create model_df_rem
model_df_rem <- structure(list(protein_18 = c(
  17.316272691676, 16.9847967072082,
  8.47167521439204, 15.5617652636277
), protein_24 = c(
  3.8073549220576,
  12.4553272203046, 13.5992154950177, 17.5539742084579
), condition = structure(c(
  2L,
  2L, 1L, 1L
), levels = c("Case", "Control"), class = "factor")),
class = "data.frame", row.names = c(
  NA,
  -4L
)
)


# Test if function works
test_that("rem_feature works", {
  suppressMessages(
    rem_p <- rem_feature(
      test_model_df,
      "protein_78"
    )
  )
  expect_equal(rem_p, model_df_rem)
})

# Testing split_data-------------------------------------------------------------
# Create a fake model_df with a subset of iris data
s1 <- iris[1:10, ]
s2 <- iris[51:61, ]
fake_model_df <- as.data.frame(rbind(s1, s2))
colnames(fake_model_df) <- c(
  "protein_1",
  "protein_2",
  "protein_3",
  "protein_4",
  "condition"
)
fake_model_df$condition <- gsub("versicolor", "Case", fake_model_df$condition)
fake_model_df$condition <- gsub("setosa", "Control", fake_model_df$condition)

# Create test_split_df
test_split_df <- list(training = structure(list(protein_1 = c(
  5.1, 4.7, 4.6, 5,
  5.4, 5, 4.9, 7, 6.4, 6.9, 5.5, 5.7, 4.9, 6.6, 5.2
), protein_2 = c(
  3.5,
  3.2, 3.1, 3.6, 3.9, 3.4, 3.1, 3.2, 3.2, 3.1, 2.3, 2.8, 2.4, 2.9,
  2.7
), protein_3 = c(
  1.4, 1.3, 1.5, 1.4, 1.7, 1.5, 1.5, 4.7, 4.5,
  4.9, 4, 4.5, 3.3, 4.6, 3.9
), protein_4 = c(
  0.2, 0.2, 0.2, 0.2,
  0.4, 0.2, 0.1, 1.4, 1.5, 1.5, 1.3, 1.3, 1, 1.3, 1.4
), condition = c(
  "Control",
  "Control", "Control", "Control", "Control", "Control", "Control",
  "Case", "Case", "Case", "Case", "Case", "Case", "Case", "Case"
)), row.names = c(NA, -15L), class = "data.frame"), test = structure(list(
  protein_1 = c(4.9, 4.6, 4.4, 6.5, 6.3, 5), protein_2 = c(
    3,
    3.4, 2.9, 2.8, 3.3, 2
  ), protein_3 = c(
    1.4, 1.4, 1.4, 4.6,
    4.7, 3.5
  ), protein_4 = c(0.2, 0.3, 0.2, 1.5, 1.6, 1), condition = c(
    "Control",
    "Control", "Control", "Case", "Case", "Case"
  )
), row.names = c(
  NA,
  -6L
), class = "data.frame"))

# Test if function works
test_that("split_data works", {
  suppressMessages(expect_equal(
    split_data(fake_model_df, train_size = 0.7, seed = 8314), test_split_df
  ))
})

# Testing train_models----------------------------------------------------------
rf_conf <- 8
glm_coef_length <- 5
svm_acc <- round(0.962963, digits = 4)


# Test if function works
test_that("train_models works", {
  suppressWarnings(suppressMessages(model_list <- train_models(test_split_df,
    algorithm_list = c("rf", "glm", "svmRadial"), seed = 351
  )))

  rf <- model_list$rf$finalModel$confusion[1, 1]
  glm <- length(model_list$glm$finalModel$coefficients)
  svm <- round(model_list$svmRadial$results$Accuracy[1], digits = 4)

  expect_equal(rf, rf_conf)
  expect_equal(glm, glm_coef_length)
  expect_equal(svm, svm_acc)
})

# Testing test_models------------------------------------------------------------
# Load test_model_list
test_model_list <- readRDS("./testdata/test_model_list.RDS")

# Create a test_prob_list
test_prob_list <- list(glm = structure(list(Case = c(
  4.25085500133093e-09, 4.41258141137268e-12,
  5.50266587850956e-09, 1, 0.999999999999469, 1
), Control = c(
  0.999999995749145,
  0.999999999995587, 0.999999994497334, 2.22044604925031e-16,
  5.31303986722389e-13,
  2.22044604925031e-16
)), class = "data.frame", row.names = c(
  "1",
  "2", "3", "4", "5", "6"
)), rf = structure(list(Case = c(
  0.058,
  0, 0.07, 1, 0.986, 0.912
), Control = c(
  0.942, 1, 0.93, 0, 0.014,
  0.088
)), class = "data.frame", row.names = c(
  "1", "2", "3", "4",
  "5", "6"
)), svmRadial = structure(list(Case = c(
  0.208990430665199,
  0.147462988884377, 0.395826109932927, 0.903842832194967, 0.891552398156344,
  0.806926899565157
), Control = c(
  0.791009569334801, 0.852537011115623,
  0.604173890067073, 0.0961571678050328, 0.108447601843656, 0.193073100434843
)), class = "data.frame", row.names = c(NA, -6L)))


# Test if function works
testthat::test_that("test_models works", {
  prob_list <- test_models(
    test_model_list,
    test_split_df
  )
  testthat::expect_equal(prob_list, test_prob_list)
})
