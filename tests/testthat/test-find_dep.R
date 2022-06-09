#Testing find_dep function
#Load test norm_df object (subset of 50 proteins from original ecoli data set)
test_norm_df <- read.csv( "./testdata/test_ecoli_norm.txt",
                         sep = "\t",
                         stringsAsFactors = TRUE)
test_norm_df <- as.matrix(test_norm_df)


#Load test fit_df object.
test_fit_df <- readRDS(file = "./testdata/test_ecoli_fit.RDS")


#Test if the function works
test_that("find_dep works", {
  expect_equal(find_dep(test_norm_df), test_fit_df)

})
