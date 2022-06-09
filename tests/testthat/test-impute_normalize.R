#Missing data imputation
#Load data with NA
test_rawdf_na <- read.csv( "./testdata/test_rawdf_na.txt",
                        sep = "\t",
                        stringsAsFactors = TRUE)
#Load imputed data (minprob)
test_imp_minprob <- read.csv( "./testdata/test_imputed_minprob.txt",
                           sep = "\t",
                           stringsAsFactors = TRUE)


test_that("impute_na works for minprob", {
  expect_equal(impute_na(test_rawdf_na), test_imp_minprob)

})

#Load imputed data (mindet)
test_imp_mindet <- read.csv( "./testdata/test_imputed_mindet.txt",
                              sep = "\t",
                              stringsAsFactors = TRUE)


test_that("impute_na works for mindet", {
  expect_equal(impute_na(test_rawdf_na, method = "minDet"), test_imp_mindet)

})

#Load imputed data (rf)
test_imp_rf <- read.csv( "./testdata/test_imputed_rf.txt",
                              sep = "\t",
                              stringsAsFactors = TRUE)


test_that("impute_na works for rf", {
  expect_equal(impute_na(test_rawdf_na, method = "RF"), test_imp_rf)

})

#Load imputed data (knn)
test_imp_knn <- read.csv( "./testdata/test_imputed_knn.txt",
                              sep = "\t",
                              stringsAsFactors = TRUE)


test_that("impute_na works for knn", {
  expect_equal(impute_na(test_rawdf_na, method = "kNN"), test_imp_knn)

})

#Load imputed data (svd)
test_imp_svd <- read.csv( "./testdata/test_imputed_svd.txt",
                              sep = "\t",
                              stringsAsFactors = TRUE)


test_that("impute_na works for svd", {
  expect_equal(impute_na(test_rawdf_na, method = "SVD"), test_imp_svd)

})

#Normalization
#Load normalized data - quantile method
test_norm_q <- read.csv( "./testdata/test_norm_q.txt",
                           sep = "\t",
                           stringsAsFactors = TRUE)

test_that("normalization works for quantile", {
  expect_equal(normalize_data(test_imp_minprob),
               test_norm_q)

})

#Load normalized data - scale method
test_norm_s <- read.csv( "./testdata/test_norm_s.txt",
                         sep = "\t",
                         stringsAsFactors = TRUE)

test_that("normalization works for sacle", {
  expect_equal(normalize_data(test_imp_minprob,
                              method = "scale"),
               test_norm_s)

})

#Load normalized data - cyclicloess method
test_norm_c <- read.csv( "./testdata/test_norm_c.txt",
                         sep = "\t",
                         stringsAsFactors = TRUE)

test_that("normalization works for cyclicloess", {
  expect_equal(normalize_data(test_imp_minprob,
                              method = "cyclicloess"),
               test_norm_c)

})
