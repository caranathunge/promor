# Pre-processing----------------------------------------------------------------
## No tech-reps-----------------------------------------------------------------
# Create test_rawdf
test_raw_df <- structure(
  list(H_1 = c(
    27.3134575527802, 25.0475589619685, 25.7489360071628,
    22.3217409656277, 23.9379563661029, 25.768729719086, 27.6466382487042,
    NA, 26.2290168009068, 20.7899474944062, 20.6427110188685, NA,
    23.3404162709489, 23.9196169668074, 23.9141515162468
  ), H_2 = c(
    27.3965634965896,
    25.4638231809381, 25.7198105464424, 23.7216107357952, 23.7859632470773,
    25.1020545097735, 27.8039664255518, 21.3147336941717, 26.6504781237034,
    20.3212455062761, NA, NA, 23.5160500915308, 24.1067737781079,
    23.7869603409319
  ), H_3 = c(
    27.0682521470214, 24.5609251894038,
    25.8093128192732, 23.0995044870265, 23.4582634148661, 24.3433337543312,
    27.7661660273555, 20.9903006674839, 26.9776011840784, 21.7031583967283,
    NA, NA, 23.4146749142913, NA, 23.9672804455828
  ), L_4 = c(
    27.5967623420646,
    26.1093667962031, 26.3575984475757, 21.2096722968756, 24.2172666952823,
    24.7096194080498, 27.915816116319, 22.5738542523838, 27.3090394121738,
    23.0277123297888, 21.4292575700939, NA, 24.4428468344072, 25.0422241309514,
    24.9093979862761
  ), L_5 = c(
    27.8685009058466, 26.1012931825391,
    26.5974729242791, 23.7620267328653, 24.9532259575982, 25.0495544486124,
    27.9875559091121, 21.4550803438114, 27.2600639899901, 22.5863936035805,
    21.7253396278684, 23.3327894314844, 24.0624994391506, 25.141762588201,
    24.1083700330184
  ), L_6 = c(
    27.9516039906164, 25.7985259878485,
    26.452446499214, 22.6047857618104, 24.2072461353296, 24.6898717464561,
    27.8163496237955, 21.4407576690532, 27.1474412722278, 22.4058228851568,
    21.2955245730356, NA, 24.2565231471014, 24.651792303055, 24.2982404854359
  )),
  row.names = c(
    "A0AV96;B7Z8Z7;A0AV96-2;D6R9D6", "A0AVT1;A0AVT1-2",
    "H7BXI1;A0FGR8-6;A0FGR8-2;A0FGR8;C9JGI7;A0FGR8-4", "A0JLT2;A0JLT2-2",
    "A0MZ66-4;A0MZ66;A0MZ66-3;B7Z7Z9;A0MZ66-5;A0MZ66-6;A0MZ66-2",
    "A1A528;O43264", "A1L0T0;E9PJS0",
    "A1L491;Q86XK2;Q86XK2-5;Q86XK2-2;H0YAV3;Q86XK2-3",
    "Q99798;A2A274;F5H2A5;B4DW08;B4DEC3",
    "Q86X10-2;A2A2F0;Q86X10-3;A2A2E9;Q86X10;B4E2E8;Q86X10-4",
    "Q9Y672;A2A2G4", "Q9Y312;A2A2Q9", "A2A2V2;P42696;P42696-2;Q5TCT4",
    "P35611-2;E7EV99;Q86XM2;A2A3N8;E7ENY0;P35611;P35611-3;H0Y9H2;Q96D30",
    "A2RRF3;Q9UBC2;Q9UBC2-2"
  ), class = "data.frame"
)

# Check if the function works
test_that("create_df works for data without tech reps", {
  expect_equal(
    create_df(
      "./testdata/test_protgroups.txt",
      "./testdata/test_expdesign.txt"
    ),
    test_raw_df
  )
})


# Filter out rows with high levels of NA in each group
# Create test_rawdf_filtered
test_raw_df_filt <- structure(list(H_1 = c(
  27.3134575527802, 25.0475589619685, 25.7489360071628,
  22.3217409656277, 23.9379563661029, 25.768729719086, 27.6466382487042,
  26.2290168009068, 20.7899474944062, 23.3404162709489, 23.9141515162468
), H_2 = c(
  27.3965634965896, 25.4638231809381, 25.7198105464424,
  23.7216107357952, 23.7859632470773, 25.1020545097735, 27.8039664255518,
  26.6504781237034, 20.3212455062761, 23.5160500915308, 23.7869603409319
), H_3 = c(
  27.0682521470214, 24.5609251894038, 25.8093128192732,
  23.0995044870265, 23.4582634148661, 24.3433337543312, 27.7661660273555,
  26.9776011840784, 21.7031583967283, 23.4146749142913, 23.9672804455828
), L_4 = c(
  27.5967623420646, 26.1093667962031, 26.3575984475757,
  21.2096722968756, 24.2172666952823, 24.7096194080498, 27.915816116319,
  27.3090394121738, 23.0277123297888, 24.4428468344072, 24.9093979862761
), L_5 = c(
  27.8685009058466, 26.1012931825391, 26.5974729242791,
  23.7620267328653, 24.9532259575982, 25.0495544486124, 27.9875559091121,
  27.2600639899901, 22.5863936035805, 24.0624994391506, 24.1083700330184
), L_6 = c(
  27.9516039906164, 25.7985259878485, 26.452446499214,
  22.6047857618104, 24.2072461353296, 24.6898717464561, 27.8163496237955,
  27.1474412722278, 22.4058228851568, 24.2565231471014, 24.2982404854359
)), row.names = c(
  "A0AV96;B7Z8Z7;A0AV96-2;D6R9D6", "A0AVT1;A0AVT1-2",
  "H7BXI1;A0FGR8-6;A0FGR8-2;A0FGR8;C9JGI7;A0FGR8-4", "A0JLT2;A0JLT2-2",
  "A0MZ66-4;A0MZ66;A0MZ66-3;B7Z7Z9;A0MZ66-5;A0MZ66-6;A0MZ66-2",
  "A1A528;O43264", "A1L0T0;E9PJS0", "Q99798;A2A274;F5H2A5;B4DW08;B4DEC3",
  "Q86X10-2;A2A2F0;Q86X10-3;A2A2E9;Q86X10;B4E2E8;Q86X10-4",
  "A2A2V2;P42696;P42696-2;Q5TCT4",
  "A2RRF3;Q9UBC2;Q9UBC2-2"
), class = "data.frame")

# Check if the function works
test_that("filterbygroup_na works ", {
  expect_equal(
    filterbygroup_na(test_raw_df),
    test_raw_df_filt
  )
})


## Data with tech-reps----------------------------------------------------------
# Create test_rawdf_tr
test_raw_df_tr <- structure(list(WT_4_1 = c(
  NA, NA, 33.4466523714569, NA, NA, NA,
  NA, NA, NA, NA
), WT_4_2 = c(
  NA, NA, 33.4031167803233, NA, NA,
  NA, NA, NA, NA, NA
), WT_4_3 = c(
  NA, NA, 33.3963695509796, NA,
  NA, NA, NA, NA, NA, NA
), WT_5_1 = c(
  NA, NA, 33.6075256022659,
  NA, 27.7291004088961, 29.8983623872195, 28.7191446901938, NA,
  NA, NA
), WT_5_2 = c(
  NA, NA, 33.5700044833146, NA, 27.9750880629869,
  29.8886270681255, 29.7229268435248, NA, NA, NA
), WT_5_3 = c(
  NA,
  NA, 33.558874449839, NA, 28.0128123150503, 30.071759272516, 29.2723088499544,
  NA, NA, NA
), WT_6_1 = c(
  NA, NA, 33.3597978675293, NA, NA, NA,
  NA, NA, NA, NA
), WT_6_2 = c(
  NA, NA, 33.340296349835, NA, NA,
  NA, NA, NA, NA, NA
), WT_6_3 = c(
  NA, NA, 33.1261253775022, NA,
  NA, NA, NA, NA, NA, NA
), D8_1_1 = c(
  NA, NA, 33.6397902993625,
  NA, NA, NA, NA, NA, NA, NA
), D8_1_2 = c(
  NA, NA, 33.6331999260124,
  NA, NA, NA, NA, NA, NA, NA
), D8_1_3 = c(
  NA, NA, 33.6451698828979,
  NA, NA, NA, NA, NA, NA, NA
), D8_2_1 = c(
  NA, NA, 33.732064069791,
  NA, NA, NA, NA, NA, NA, NA
), D8_2_2 = c(
  NA, NA, 33.8470745315561,
  NA, NA, NA, NA, NA, NA, NA
), D8_2_3 = c(
  NA, NA, 34.034364466028,
  NA, NA, NA, NA, NA, NA, NA
), D8_3_1 = c(
  32.651534866107, 29.9718582903499,
  33.7256797804011, 25.8938399462174, NA, NA, 28.0640425936895,
  29.8058276204678, 28.1101836284433, 32.291235781689
), D8_3_2 = c(
  32.7083716870042,
  29.9744590980922, 33.6990288832398, 25.8779251409578, NA, NA,
  27.7160725522537, 30.0867662850199, 28.1139629227287, 32.3195320217604
), D8_3_3 = c(
  32.7525053859095, 29.8938139141624, 33.6982009971023,
  26.1088485652935, 26.9912915696221, NA, 27.8317712225526, 29.8307290910728,
  28.2692796614525, 32.2906592143389
)), row.names = c(
  "CON__O43790",
  "CON__O76009;CON__Q6NTB9", "CON__P00761", "CON__P48668;CON__P02538;CON__P04259",
  "CON__P04264", "CON__P13645", "CON__P35908", "CON__P78386", "CON__Q14525",
  "CON__Q9UE12;CON__Q15323"
), class = "data.frame")

# Check if the function works
test_that("create_df works for data with tech reps", {
  expect_equal(
    create_df("./testdata/test_protgroups_tr.txt",
      "./testdata/test_expdesign_tr.txt",
      tech_reps = TRUE
    ),
    test_raw_df_tr
  )
})

# Average across tech reps
# Create test_rawdf_avg_tr

test_raw_df_tr_avg <- structure(list(WT_4 = c(
  NaN, NaN, 33.4153795675866, NaN, NaN,
  NaN, NaN, NaN, NaN, NaN
), WT_5 = c(
  NaN, NaN, 33.5788015118065,
  NaN, 27.9056669289778, 29.9529162426203, 29.2381267945577, NaN,
  NaN, NaN
), WT_6 = c(
  NaN, NaN, 33.2754065316221, NaN, NaN, NaN,
  NaN, NaN, NaN, NaN
), D8_1 = c(
  NaN, NaN, 33.6393867027576, NaN,
  NaN, NaN, NaN, NaN, NaN, NaN
), D8_2 = c(
  NaN, NaN, 33.871167689125,
  NaN, NaN, NaN, NaN, NaN, NaN, NaN
), D8_3 = c(
  32.7041373130069,
  29.9467104342015, 33.7076365535811, 25.9602045508229, 26.9912915696221,
  NaN, 27.8706287894986, 29.9077743321868, 28.1644754042081, 32.3004756725961
)), class = "data.frame", row.names = c(
  "CON__O43790", "CON__O76009;CON__Q6NTB9",
  "CON__P00761", "CON__P48668;CON__P02538;CON__P04259", "CON__P04264",
  "CON__P13645", "CON__P35908", "CON__P78386", "CON__Q14525",
  "CON__Q9UE12;CON__Q15323"
))

# Check if the function works
test_that("aver_techreps works", {
  expect_equal(
    aver_techreps(test_raw_df_tr),
    test_raw_df_tr_avg
  )
})

## Removing samples-------------------------------------------------------------
# Make a test data frame
df <- data.frame(
  a = c(12, 34, 11, 83, NA, 67),
  b = c(NA, NA, 39, 29, 20, NA),
  c = c(82, 40, 22, NA, NA, 70)
)
# Create df_rem
df_rem <- structure(list(a = c(12, 34, 11, 83, NA, 67), c = c(
  82, 40, 22,
  NA, NA, 70
)), class = "data.frame", row.names = c(NA, -6L))


# Check if the function works
test_that("rem_sample works", {
  expect_equal(
    rem_sample(df, "b"),
    df_rem
  )
})

# Missing data imputation--------------------------------------------------------

# minprob results
df_mp <- structure(
  list(
    a = c(12, 34, 11, 83, 9.80950652319248, 67),
    b = c(12.6260927833319, 29.8391137373841, 39, 29, 20, -1.37591651207944),
    c = c(82, 40, 22, 14.7098046749082, 22.299106363566, 70)
  ),
  row.names = c(NA, -6L),
  class = "data.frame"
)

# rf - results
df_rf <- structure(
  list(
    a = c(12, 34, 11, 83, 62.1633333333333, 67),
    b = c(26.6, 32.15, 39, 29, 20, 22.3),
    c = c(82, 40, 22, 53.125, 66.025, 70)
  ),
  row.names = c(NA, -6L),
  class = "data.frame"
)

# knn - results
df_knn <- structure(
  list(
    a = c(12, 34, 11, 83, 34, 67),
    b = c(39, 39, 39, 29, 20, 39),
    c = c(82, 40, 22, 70, 70, 70)
  ),
  row.names = c("1", "2", "3", "4", "5", "6"),
  class = "data.frame"
)

# svd - results
df_svd <- structure(
  list(
    a = c(12, 34, 11, 83, 40.5316879411941, 67),
    b = c(15.8432402677882, 34.3000432486029, 39, 29, 20, 24.3678401876791),
    c = c(82, 40, 22, 59.6698946596492, 73.7190158521945, 70)
  ),
  class = "data.frame",
  row.names = c(NA, -6L)
)

# mindet - results
df_md <- structure(
  list(
    a = c(12, 34, 11, 83, 11.04, 67),
    b = c(20.18, 20.18, 39, 29, 20, 20.18),
    c = c(82, 40, 22, 22.54, 22.54, 70)
  ),
  row.names = c(NA, -6L),
  class = "data.frame"
)

# testing minprob
test_that("impute_na works for minprob", {
  expect_equal(impute_na(df), df_mp)
})

# testing rf
test_that("impute_na works for rf", {
  suppressWarnings(expect_equal(impute_na(df,
    method = "RF"
  ), df_rf))
})

# testing knn
test_that("impute_na works for knn", {
  expect_equal(impute_na(df,
    method = "kNN"
  ), df_knn)
})

# testing svd
test_that("impute_na works for svd", {
  suppressMessages(expect_equal(impute_na(df,
    method = "SVD"
  ), df_svd))
})

# testing minDet
test_that("impute_na works for minDet", {
  expect_equal(impute_na(df,
    method = "minDet"
  ), df_md)
})

# Normalization------------------------------------------------------------------
# quantile - results
df_norm_q <- structure(c(
  18.0997021211887, 34.3333333333333, 15.208697594444,
  68, 7.71446489534041, 55.6130379124614, 15.208697594444,
  55.6130379124614,
  68, 34.3333333333333, 18.0997021211887, 7.71446489534041, 68,
  34.3333333333333, 15.208697594444, 7.71446489534041,
  18.0997021211887,
  55.6130379124614
),
dim = c(6L, 3L),
dimnames = list(NULL, c("a", "b", "c"))
)

# scale -results
df_norm_s <- structure(c(
  13.5592392536598, 38.4178445520362, 12.4293026491882,
  93.7847381711472, 11.0841204923586, 75.7057524996008,
  13.3932139334231,
  31.6520431717492, 41.3695156820846, 30.7619475584732,
  21.2151362472229,
  -1.45951281342845, 68.4138364231764, 33.3726031332568,
  18.3549317232912,
  12.2726118395859, 18.6044806724393, 58.4020554831994
),
dim = c(6L, 3L),
dimnames = list(NULL, c("a", "b", "c"))
)

# cyclicloess - results
df_norm_c <- structure(c(
  8.87594722638867, 35.9180013309211, 24.1612907491603,
  43.5415649767626, 17.3396297815539, 44.825964199607,
  21.9075705937958,
  38.5487769581834, 24.1612907491603, 43.5415649767626,
  17.3396297815539,
  44.825964199607, 81.9288096106512, 35.9180013309212,
  24.1612907491603,
  43.5415649767626, 17.3396297815539, 44.825964199607
),
dim = c(6L, 3L),
dimnames = list(NULL, c("a", "b", "c"))
)


## testing - quantile method----------------------------------------------------

test_that("normalization works for quantile", {
  expect_equal(
    normalize_data(df_mp),
    df_norm_q
  )
})

## testing - scale method-------------------------------------------------------

test_that("normalization works for sacle", {
  expect_equal(
    normalize_data(df_mp,
      method = "scale"
    ),
    df_norm_s
  )
})

## testing- cyclicloess method--------------------------------------------------

test_that("normalization works for cyclicloess", {
  expect_equal(
    normalize_data(df_mp,
      method = "cyclicloess"
    ),
    df_norm_c
  )
})

# Differential expression -------------------------------------------------------
# make a fake - norm_df
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

# Load fit_df
fake_fit_df <- readRDS("./testdata/fake_fit.RDS")

# Test find_dep works
test_that("find_dep works", {
  expect_equal(
    find_dep(fake_norm_df),
    fake_fit_df
  )
})
