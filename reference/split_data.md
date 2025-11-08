# Split the data frame to create training and test data

This function can be used to create balanced splits of the protein
intensity data in a `model_df` object to create training and test data

## Usage

``` r
split_data(model_df, train_size = 0.8, seed = NULL)
```

## Arguments

- model_df:

  A `model_df` object from performing `pre_process`.

- train_size:

  The size of the training data set as a proportion of the complete data
  set. Default is 0.8.

- seed:

  Numerical. Random number seed. Default is `NULL`

## Value

A list of data frames.

## Details

This function splits the `model_df` object in to training and test data
sets using random sampling while preserving the original class
distribution of the data. Make sure to fix the random number seed with
`seed` for reproducibility

## See also

- `pre_process`

- [`createDataPartition`](https://rdrr.io/pkg/caret/man/createDataPartition.html)

## Author

Chathurani Ranathunge

## Examples

``` r
## Create a model_df object
covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#> Total number of differentially expressed proteins (8) is less than n_top.
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.

## Split the data frame into training and test data sets using default settings
covid_split_df1 <- split_data(covid_model_df, seed = 8314)

## Split the data frame into training and test data sets with 70% of the
## data in training and 30% in test data sets
covid_split_df2 <- split_data(covid_model_df, train_size = 0.7, seed = 8314)

## Access training data set
covid_split_df1$training
#>    sp|O00391|QSOX1_HUMAN sp|P00746|CFAD_HUMAN sp|P02652|APOA2_HUMAN
#> 1               22.40059             24.73720              28.66292
#> 2               19.97970             26.80882              30.67527
#> 3               21.12426             22.82051              30.16657
#> 4               20.21675             27.64256              30.67527
#> 5               20.92041             25.60917              31.48626
#> 6               23.37177             22.71174              27.64256
#> 7               22.27765             26.38049              29.94831
#> 8               21.56619             24.97819              30.51976
#> 9               24.73720             23.57248              30.00715
#> 10              22.10170             23.73261              28.73758
#> 11              22.30244             23.87346              29.41927
#> 12              23.73261             26.60611              28.51444
#> 13              22.62693             21.62540              29.33404
#> 14              21.78974             24.47092              29.49655
#> 15              23.15140             24.32863              27.82097
#> 16              19.37365             24.06332              29.41927
#> 17              21.81120             22.19998              27.76262
#> 18              20.49465             25.23133              30.58588
#> 19              20.80441             24.32863              28.14984
#> 20              19.02912             24.32863              30.58588
#> 21              22.11671             24.42514              28.42881
#> 22              22.11671             25.32291              30.00715
#> 23              22.19998             24.79335              31.38459
#> 24              20.99572             25.73887              31.68780
#> 25              22.34632             25.78791              30.44035
#> 26              21.07625             26.97260              29.41927
#> 27              22.11671             23.62285              29.33404
#> 28              22.36332             25.50697              30.16657
#> 29              20.17148             25.45922              31.38459
#>    sp|P02765|FETUA_HUMAN sp|P13796|PLSL_HUMAN sp|P22352|GPX3_HUMAN
#> 1               30.67527             22.57517             23.96220
#> 2               32.70218             17.63880             22.86492
#> 3               29.19797             25.13590             23.92030
#> 4               31.68780             18.80999             21.43144
#> 5               31.96064             19.37365             23.31727
#> 6               28.84603             26.09616             25.73887
#> 7               31.08519             21.74014             23.26651
#> 8               30.73758             19.02912             23.82412
#> 9               30.07485             23.31727             25.55824
#> 10              28.84603             21.90714             24.32863
#> 11              30.73758             19.22486             24.37980
#> 12              30.58588             22.86492             22.95628
#> 13              30.35700             23.77481             25.45922
#> 14              31.68780             19.97970             23.52849
#> 15              29.65531             25.83668             24.47092
#> 16              31.48626             20.59215             21.90714
#> 17              29.03320             24.61318             24.94399
#> 18              31.96064             20.54634             23.38986
#> 19              30.44035             24.16344             22.71174
#> 20              31.08519             20.17148             19.37365
#> 21              30.35700             20.80441             24.32863
#> 22              31.59221             21.12426             23.92030
#> 23              32.17057             21.66673             22.38032
#> 24              31.96064             20.99572             23.31727
#> 25              30.26855             15.86089             23.35368
#> 26              30.87044             22.30244             20.54634
#> 27              30.87044             22.57517             24.20810
#> 28              30.87044             20.08515             22.08669
#> 29              31.28625             20.87460             23.62285
#>    sp|P25311|ZA2G_HUMAN sp|Q16610|ECM1_HUMAN  condition
#> 1              27.70882             25.18199     Severe
#> 2              25.23133             27.22918 Non.Severe
#> 3              27.87499             24.97819     Severe
#> 4              24.84503             28.42881 Non.Severe
#> 5              26.49553             27.08827     Severe
#> 6              28.51444             23.87346     Severe
#> 7              25.93514             26.60611 Non.Severe
#> 8              26.04292             26.19652 Non.Severe
#> 9              27.03502             26.74563     Severe
#> 10             25.69315             24.24620     Severe
#> 11             26.54609             26.04292     Severe
#> 12             27.43537             26.04292     Severe
#> 13             26.49553             25.02850     Severe
#> 14             24.79335             27.39454     Severe
#> 15             28.51444             24.84503     Severe
#> 16             24.94399             26.29915 Non.Severe
#> 17             27.22918             24.79335     Severe
#> 18             25.98335             26.97260 Non.Severe
#> 19             27.03502             26.60611     Severe
#> 20             25.32291             25.78791 Non.Severe
#> 21             26.34214             26.13770     Severe
#> 22             26.92763             26.74563     Severe
#> 23             24.56294             27.76262 Non.Severe
#> 24             25.50697             27.33912 Non.Severe
#> 25             26.43611             27.12807 Non.Severe
#> 26             27.18300             27.54357 Non.Severe
#> 27             25.88645             25.23133 Non.Severe
#> 28             26.04292             27.48114 Non.Severe
#> 29             26.19652             26.65823 Non.Severe

## Access test data set
covid_split_df1$test
#>   sp|O00391|QSOX1_HUMAN sp|P00746|CFAD_HUMAN sp|P02652|APOA2_HUMAN
#> 1              21.28517             25.93514              31.68780
#> 2              21.46057             25.45922              29.41927
#> 3              21.81120             25.36447              29.94831
#> 4              22.88773             23.67558              28.84603
#> 5              20.72842             26.09616              30.73758
#> 6              22.01501             24.84503              29.25736
#>   sp|P02765|FETUA_HUMAN sp|P13796|PLSL_HUMAN sp|P22352|GPX3_HUMAN
#> 1              31.79321             18.19499             24.32863
#> 2              30.16657             21.07625             22.95628
#> 3              30.96996             19.37365             22.98936
#> 4              29.33404             23.20021             25.93514
#> 5              32.09807             21.35759             21.19299
#> 6              29.71094             19.86997             24.61318
#>   sp|P25311|ZA2G_HUMAN sp|Q16610|ECM1_HUMAN  condition
#> 1             25.78791             26.65823 Non.Severe
#> 2             28.73758             26.65823     Severe
#> 3             27.08827             26.87214     Severe
#> 4             27.92312             25.45922     Severe
#> 5             24.51383             27.70882 Non.Severe
#> 6             26.29915             26.24771 Non.Severe
```
