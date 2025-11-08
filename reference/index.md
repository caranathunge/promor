# Package index

## Proteomics Data Analysis

### Pre-processing

- [`create_df()`](https://caranathunge.github.io/promor/reference/create_df.md)
  : Create a data frame of protein intensities
- [`aver_techreps()`](https://caranathunge.github.io/promor/reference/aver_techreps.md)
  : Compute average intensity

### Quality Control

- [`filterbygroup_na()`](https://caranathunge.github.io/promor/reference/filterbygroup_na.md)
  : Filter proteins by group level missing data
- [`rem_sample()`](https://caranathunge.github.io/promor/reference/rem_sample.md)
  : Remove user-specified samples
- [`impute_na()`](https://caranathunge.github.io/promor/reference/impute_na.md)
  : Impute missing values
- [`normalize_data()`](https://caranathunge.github.io/promor/reference/normalize_data.md)
  : Normalize intensity data

### Differential Expression Analysis

- [`find_dep()`](https://caranathunge.github.io/promor/reference/find_dep.md)
  : Identify differentially expressed proteins between groups

### Visualization

- [`corr_plot()`](https://caranathunge.github.io/promor/reference/corr_plot.md)
  : Correlation between technical replicates
- [`heatmap_na()`](https://caranathunge.github.io/promor/reference/heatmap_na.md)
  : Visualize missing data
- [`impute_plot()`](https://caranathunge.github.io/promor/reference/impute_plot.md)
  : Visualize the impact of imputation
- [`norm_plot()`](https://caranathunge.github.io/promor/reference/norm_plot.md)
  : Visualize the effect of normalization
- [`volcano_plot()`](https://caranathunge.github.io/promor/reference/volcano_plot.md)
  : Volcano plot
- [`heatmap_de()`](https://caranathunge.github.io/promor/reference/heatmap_de.md)
  : Heatmap of differentially expressed proteins

### Miscellaneous

- [`onegroup_only()`](https://caranathunge.github.io/promor/reference/onegroup_only.md)
  : Proteins that are only expressed in a given group

## Building Models

### Pre-processing

- [`pre_process()`](https://caranathunge.github.io/promor/reference/pre_process.md)
  : Pre-process protein intensity data for modeling
- [`rem_feature()`](https://caranathunge.github.io/promor/reference/rem_feature.md)
  : Remove user-specified proteins (features) from a data frame
- [`split_data()`](https://caranathunge.github.io/promor/reference/split_data.md)
  : Split the data frame to create training and test data

### Modeling

- [`train_models()`](https://caranathunge.github.io/promor/reference/train_models.md)
  : Train machine learning models on training data
- [`test_models()`](https://caranathunge.github.io/promor/reference/test_models.md)
  : Test machine learning models on test data

### Visualization

- [`feature_plot()`](https://caranathunge.github.io/promor/reference/feature_plot.md)
  : Visualize feature (protein) variation among conditions
- [`varimp_plot()`](https://caranathunge.github.io/promor/reference/varimp_plot.md)
  : Variable importance plot
- [`performance_plot()`](https://caranathunge.github.io/promor/reference/performance_plot.md)
  : Model performance plot
- [`roc_plot()`](https://caranathunge.github.io/promor/reference/roc_plot.md)
  : ROC plot

## Data

- [`covid_fit_df`](https://caranathunge.github.io/promor/reference/covid_fit_df.md)
  : Suvarna et al 2021 LFQ data (fit object)
- [`covid_norm_df`](https://caranathunge.github.io/promor/reference/covid_norm_df.md)
  : Suvarna et al 2021 LFQ data (normalized)
- [`ecoli_norm_df`](https://caranathunge.github.io/promor/reference/ecoli_norm_df.md)
  : Cox et al 2014 LFQ data (normalized)
- [`ecoli_fit_df`](https://caranathunge.github.io/promor/reference/ecoli_fit_df.md)
  : Cox et al 2014 LFQ data (fit object)
