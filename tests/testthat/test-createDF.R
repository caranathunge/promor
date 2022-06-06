#No tech-reps
#Load test_rawdf
test_rawdf <- read.csv( "./testdata/test_rawdf.txt",
                             sep = "\t",
                             stringsAsFactors = TRUE)

#Check if the function works
test_that("create_df works for data without tech reps", {
  expect_equal(create_df("./testdata/test_protgroups.txt",
                         "./testdata/test_expdesign.txt"),
               test_rawdf)
})


#Filter out rows with high levels of NA in each group
#Load test_rawdf_filtered
test_rawdf_filtered <- read.csv( "./testdata/test_rawdf_filtered.txt",
                        sep = "\t",
                        stringsAsFactors = TRUE)

#Check if the function works
test_that("filterbygroup_na works ", {
  expect_equal(filterbygroup_na(test_rawdf),
               test_rawdf_filtered)
})


#Data with tech-reps
#Load test_rawdf_tr
test_rawdf_tr <- read.csv( "./testdata/test_rawdf_tr.txt",
                           sep = "\t",
                           stringsAsFactors = TRUE)

#Check if the function works
test_that("create_df works for data with tech reps", {
  expect_equal(create_df("./testdata/test_protgroups_tr.txt",
                         "./testdata/test_expdesign_tr.txt",
                         tech_reps = TRUE),
               test_rawdf_tr)
})

#Average across tech reps
#Load test_rawdf_avg_tr
test_rawdf_avg_tr <- read.csv( "./testdata/test_rawdf_avg_tr.txt",
                               sep = "\t",
                               stringsAsFactors = TRUE)

#Check if the function works
test_that("aver_techreps works",{
  expect_equal(aver_techreps(test_rawdf_tr),
               test_rawdf_avg_tr)

})

#Removing samples
#Load test_rawdf_avg_tr
test_rem_sample <- read.csv( "./testdata/test_rem_sample.txt",
                               sep = "\t",
                               stringsAsFactors = TRUE)

#Check if the function works
test_that("rem_sample works",{
  expect_equal(rem_sample(test_rawdf, "L_6"),
               test_rem_sample)

})
