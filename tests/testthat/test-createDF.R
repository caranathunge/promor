#No tech-reps
#Load test_rawdf
test_rawdf <- read.csv( "./testdata/test_rawdf.txt",
                             sep = "\t",
                             stringsAsFactors = TRUE)


test_that("create_df works for data without tech reps", {
  expect_equal(create_df("./testdata/test_protgroups.txt",
                         "./testdata/test_expdesign.txt"),
               test_rawdf)
})


#With tech-reps
#Load test_rawdf_tr
test_rawdf_tr <- read.csv( "./testdata/test_rawdf_tr.txt",
                        sep = "\t",
                        stringsAsFactors = TRUE)


test_that("create_df works for data without tech reps", {
  expect_equal(create_df("./testdata/test_protgroups_tr.txt",
                         "./testdata/test_expdesign_tr.txt",
                         tech.reps = TRUE),
               test_rawdf_tr)
})
