# Run imputeLCMD function imputeMinProb
impute.Min.Prob <- function(dataSet.mvs, q = 0.01, tune_sigma = 1) {
  nSamples <- dim(dataSet.mvs)[2]
  nFeatures <- dim(dataSet.mvs)[1]
  dataSet.imputed <- dataSet.mvs
  min.samples <- apply(dataSet.imputed, 2, quantile,
    prob = q,
    na.rm = T
  )
  count.NAs <- apply(!is.na(dataSet.mvs), 1, sum)
  count.NAs <- count.NAs / nSamples
  dataSet.filtered <- dataSet.mvs[which(count.NAs > 0.5), ]
  protSD <- apply(dataSet.filtered, 1, sd)
  sd.temp <- median(protSD, na.rm = T) * tune_sigma

  for (i in 1:(nSamples)) {
    dataSet.to.impute.temp <- rnorm(nFeatures,
      mean = min.samples[i],
      sd = sd.temp
    )
    dataSet.imputed[
      which(
        is.na(
          dataSet.mvs[, i]
        )
      ), i
    ] <- dataSet.to.impute.temp[which(
      is.na(dataSet.mvs[, i])
    )]
  }
  return(dataSet.imputed)
}


# Run imputeLCMD function imputeMinDet
impute.MinDet <- function(dataSet.mvs, q = 0.01) {
  nSamples <- dim(dataSet.mvs)[2]
  dataSet.imputed <- dataSet.mvs
  lowQuantile.samples <- apply(dataSet.imputed, 2, quantile,
    prob = q, na.rm = T
  )
  for (i in 1:(nSamples)) {
    dataSet.imputed[which(
      is.na(dataSet.mvs[, i])
    ), i] <- lowQuantile.samples[i]
  }
  return(dataSet.imputed)
}
