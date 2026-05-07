# Model performance plot

This function generates plots to visualize model performance

## Usage

``` r
performance_plot(
  model_list,
  type = "box",
  text_size = 10,
  palette = "viridis",
  save = FALSE,
  file_path = NULL,
  file_name = "Performance_plot",
  file_type = "pdf",
  plot_width = 7,
  plot_height = 7,
  dpi = 80
)
```

## Arguments

- model_list:

  A `model_list` object from performing `train_models`.

- type:

  Type of plot to generate. Choices are "box" or "dot." Default is
  `"box."` for boxplots.

- text_size:

  Text size for plot labels, axis labels etc. Default is `10`.

- palette:

  Viridis color palette option for plots. Default is `"viridis"`. See
  [`viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html)
  for available options.

- save:

  Logical. If `TRUE` saves a copy of the plot in the directory provided
  in `file_path`.

- file_path:

  A string containing the directory path to save the file.

- file_name:

  File name to save the plot. Default is `"Performance_plot."`

- file_type:

  File type to save the plot. Default is `"pdf"`.

- plot_width:

  Width of the plot. Default is `7`.

- plot_height:

  Height of the plot. Default is `7`.

- dpi:

  Plot resolution. Default is `80`.

## Value

A `ggplot2` object.

## Details

- The default metrics used for classification based models are
  "Accuracy" and "Kappa."

- These metric types can be changed by providing additional arguments to
  the `train_models` function. See
  [`train`](https://rdrr.io/pkg/caret/man/train.html) and
  [`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html) for
  more information.

## See also

- `train_models`

- [`resamples`](https://rdrr.io/pkg/caret/man/resamples.html)

- [`train`](https://rdrr.io/pkg/caret/man/train.html)

- [`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html)

## Author

Chathurani Ranathunge

## Examples

``` r
# \donttest{

## Create a model_df object
covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#> Total number of differentially expressed proteins (8) is less than n_top.
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.

## Split the data frame into training and test data sets
covid_split_df <- split_data(covid_model_df)

## Fit models based on the default list of machine learning (ML) algorithms
covid_model_list <- train_models(covid_split_df)
#> 
#> Running svmRadial...
#> Loading required package: ggplot2
#> Loading required package: lattice
#> 
#> Running rf...
#> 
#> Running glm...
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> 
#> Running xgbLinear...
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep1: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep2: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10.Rep3: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: There were missing values in resampled performance measures.
#> Something is wrong; all the Accuracy metric values are missing:
#>     Accuracy       Kappa    
#>  Min.   : NA   Min.   : NA  
#>  1st Qu.: NA   1st Qu.: NA  
#>  Median : NA   Median : NA  
#>  Mean   :NaN   Mean   :NaN  
#>  3rd Qu.: NA   3rd Qu.: NA  
#>  Max.   : NA   Max.   : NA  
#>  NAs    :27    NAs    :27   
#> xgbLinear failed.
#> 
#> Running naive_bayes...
#> Done!

## Generate box plots to visualize performance of different ML algorithms
performance_plot(covid_model_list)
#> Using Resample as id variables


## Generate dot plots
performance_plot(covid_model_list, type = "dot")
#> Using Resample as id variables
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_segment()`).


## Change color palette
performance_plot(covid_model_list, type = "dot", palette = "inferno")
#> Using Resample as id variables
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_segment()`).

# }
```
