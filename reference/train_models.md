# Train machine learning models on training data

This function can be used to train models on protein intensity data
using different machine learning algorithms

## Usage

``` r
train_models(
  split_df,
  resample_method = "repeatedcv",
  resample_iterations = 10,
  num_repeats = 3,
  algorithm_list,
  seed = NULL,
  ...
)
```

## Arguments

- split_df:

  A `split_df` object from performing `split_data`.

- resample_method:

  The resampling method to use. Default is `"repeatedcv"` for repeated
  cross validation. See
  [`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html) for
  details on other available methods.

- resample_iterations:

  Number of resampling iterations. Default is `10`.

- num_repeats:

  The number of complete sets of folds to compute (For
  `resampling method = "repeatedcv"` only).

- algorithm_list:

  A list of classification or regression algorithms to use. A full list
  of machine learning algorithms available through the `caret` package
  can be found here:
  [http://topepo.github.io/caret/train-models-by-tag.html](http://topepo.github.io/caret/train-models-by-tag.md).
  See below for default options.

- seed:

  Numerical. Random number seed. Default is `NULL`

- ...:

  Additional arguments to be passed on to
  [`train`](https://rdrr.io/pkg/caret/man/train.html) function in the
  `caret` package.

## Value

A list of class `train` for each machine-learning algorithm. See
[`train`](https://rdrr.io/pkg/caret/man/train.html) for more information
on accessing different elements of this list.

## Details

- In the event that `algorithm_list` is not provided, a default list of
  four classification-based machine-learning algorithms will be used for
  building and training models. Default `algorithm_list`: "svmRadial",
  "rf", "glm", "xgbLinear, and "naive_bayes."

- Note: Models that fail to build are removed from the output.

- Make sure to fix the random number seed with `seed` for
  reproducibility

## References

Kuhn, Max. "Building predictive models in R using the caret package."
Journal of statistical software 28 (2008): 1-26.

## See also

- `pre_process`

- [`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html)

- [`train`](https://rdrr.io/pkg/caret/man/train.html)

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
covid_split_df <- split_data(covid_model_df, seed = 8314)

## Fit models based on the default list of machine learning (ML) algorithms
covid_model_list1 <- train_models(split_df = covid_split_df, seed = 351)
#> 
#> Running svmRadial...
#> 
#> Running rf...
#> 
#> Running glm...
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
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
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

## Fit models using a user-specified list of ML algorithms.
covid_model_list2 <- train_models(
  covid_split_df,
  algorithm_list = c("svmRadial", "glmboost"),
  seed = 351
)
#> 
#> Running svmRadial...
#> 
#> Running glmboost...
#> glmboost failed.
#> Done!

## Change resampling method and resampling iterations.
covid_model_list3 <- train_models(
  covid_split_df,
  resample_method = "cv",
  resample_iterations = 50,
  seed = 351
)
#> Warning: `repeats` has no meaning for this resampling method.
#> 
#> Running svmRadial...
#> Warning: There were missing values in resampled performance measures.
#> 
#> Running rf...
#> Warning: There were missing values in resampled performance measures.
#> 
#> Running glm...
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: There were missing values in resampled performance measures.
#> 
#> Running xgbLinear...
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold01: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold02: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold03: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold04: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold05: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold06: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold07: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold08: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold09: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold10: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold11: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold12: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold13: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold14: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold15: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold16: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold17: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold18: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold19: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold20: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold21: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold22: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold23: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold24: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold25: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold26: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold27: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold28: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=0e+00, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=1e-01, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=1e-04, nrounds= 50, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=0e+00, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=1e-01, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=1e-04, nrounds=100, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=0e+00, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=1e-01, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=0e+00, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-01, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
#>   ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost]
#> Warning: Argument 'objective' is only for custom objectives. For built-in objectives, pass the objective under 'params'. This warning will become an error in a future version.
#> Warning: model fit failed for Fold29: lambda=1e-04, alpha=1e-04, nrounds=150, eta=0.3 Error in modelFit$xNames <- colnames(x) : 
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
#> Warning: There were missing values in resampled performance measures.
#> Done!
# }
```
