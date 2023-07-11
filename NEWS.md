# promor 0.2.1

#### Documentation updates
* In previous versions, the functionality of the `uniq_pep` parameter 
(see `create_df` help page) was incorrectly described. We corrected the 
description to match its functionality as follows: "Numerical. Proteins that are
identified by this number or fewer number of unique peptides are filtered out 
(default is 2).Only applies when `input_type = "MaxQuant"`.

#### Bug fixes
* In previous versions, `volcano_plot` function plotted "-log~10~(P-value)" on 
the y-axis by default even when `sig = "adjP"` and only colored the dots by 
their  significance based on "adjP". In this version, `sig = "adjP"` shows 
"-log~10~.adj. P-value" on the y-axis and colors the dots by their "adjP"
significance.


# promor 0.2.0

#### New data types allowed
* A new argument (`data_type`) added to the `create_df` function to 
accommodate other types of LFQ data (raw intensity, iBAQ).

#### A new input type added
* A new argument `input_type` added to the `create_df` function to 
allow users to input data from a standard quantitative matrix.

#### Workflow changes
* To allow for missing data imputation prior to or after data normalization 
step (depending on the imputation method used), the following changes were made:
  * `norm_df` and `imp_df` arguments replaced with a generic `df` argument in 
  the functions, `find_dep`, `impute_na`, `normalize_data`, and `heatmap_de`
  * A note was added to the tutorials to clarify that for some imputation 
  methods, such as the kNN method, data normalization should be performed prior
  to imputation.
  
#### Default machine learning algorithms
* `naive_bayes` added to the default `algorithm_list` argument in the 
`train_models` function.

# promor 0.1.1

#### Bug fixes
* Fixes a minor issue with `create_df` when removing potential contaminants. 
The number of potential contaminants removed is now shown in the console.
* Fixes an issue with `find_dep` that previously used a fixed value for the
`adj_method` argument.
* Fixes an issue with the `file_path` argument for saving the "TopHits.txt" 
file produced by the `find_dep` function.

#### Other changes
* Citation file updated with the biorxiv preprint details.
* Readme file updated with information on the Shiny App.
* Help pages updated.

# promor 0.1.0

* First official release!
