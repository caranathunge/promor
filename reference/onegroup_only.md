# Proteins that are only expressed in a given group

This function outputs a list of proteins that are only expressed
(present) in one user-specified group while not expressed (completely
absent) in another user-specified group.

## Usage

``` r
onegroup_only(
  raw_df,
  abs_group,
  pres_group,
  set_na = 0.34,
  save = FALSE,
  file_path = NULL
)
```

## Arguments

- raw_df:

  A `raw_df` object (output of
  [`create_df`](https://caranathunge.github.io/promor/reference/create_df.md))

- abs_group:

  Name of the group in which proteins are not expressed.

- pres_group:

  Name of the group in which proteins are expressed.

- set_na:

  The percentage of missing data allowed in `pres_group`. Default is
  0.34 (one third of the samples in the group).

- save:

  Logical. If `TRUE` (default), it saves the output in a text file named
  "Group\_`pres_group`\_only.txt."

- file_path:

  A string containing the directory path to save the file.

## Value

A list of majority protein IDs.

## Details

Note: `onegroup_only` function assumes that column names in the `raw_df`
object provided as `df` follow "Group_UniqueSampleID" notation. (Use
`head(raw_df)` to check the structure of your `raw_df` object.)

- A text file containing majority protein IDs will be saved in a
  temporary directory if `file_path` is not specified.

## Author

Chathurani Ranathunge

## Examples

``` r
# Generate a raw_df object with default settings. No technical replicates.
raw_df <- create_df(
prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
)
#> 0 empty row(s) removed.
#> 0 empty column(s) removed.
#> 80 protein(s) (rows) only identified by site removed.
#> 65 reverse protein(s) (rows) removed.
#> 42 protein potential contaminant(s) (rows) removed.
#> 1923 protein(s) identified by 2 or fewer unique peptides removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.

## Find the proteins only expressed in group L, but absent in group H.
onegroup_only(raw_df, abs_group = "H",
pres_group = "L")
#>  [1] "H7C5A7;E9PJN0;A6NDT2;O14734;H0Y698;Q9BR14;E9PMC4;G5E995;E9PIS4"                                                                                                
#>  [2] "A8MXV4"                                                                                                                                                        
#>  [3] "A8MYT4;Q8NEB9"                                                                                                                                                 
#>  [4] "B4DXC8;Q9Y5B8-2;Q9Y5B8;E9PNU1"                                                                                                                                 
#>  [5] "B7Z468;P20810-3;E7ES10;P20810-4;E9PCH5;G5E946;E7EVY3;B7Z574;P20810-2;P20810;G5E9D3;P20810-5;E7ESM9;P20810-7;P20810-6;E7EQ12;E7EQA0;H0Y7F0;H0YD33;H0Y9H6;E9PDE4"
#>  [6] "C9J3L8;C9J5W0;P43307;E9PAL7;C9IZQ1;P43307-2;F5H5Y2"                                                                                                            
#>  [7] "Q9HAW4-2;E7ESG2;Q9HAW4-3;Q9HAW4"                                                                                                                               
#>  [8] "G3V158;Q9Y315;E9PMH9;E9PML7;E9PPM8"                                                                                                                            
#>  [9] "H3BRK5;H3BMW1;H3BNF0;Q9Y5Y2;H3BNS4"                                                                                                                            
#> [10] "P53007"                                                                                                                                                        
#> [11] "P61764;F5H5Q4;P61764-2;B7Z1V5"                                                                                                                                 
#> [12] "P81605;P81605-2"                                                                                                                                               
#> [13] "Q13228-2;Q13228;B4E1F3"                                                                                                                                        
#> [14] "Q5SW02;Q5SVZ6"                                                                                                                                                 
#> [15] "Q6Q0C0-2;Q6Q0C0"                                                                                                                                               
#> [16] "Q6UW63"                                                                                                                                                        
#> [17] "Q8IYI6"                                                                                                                                                        
#> [18] "Q92572;F5H459"                                                                                                                                                 
#> [19] "Q9C0B7"                                                                                                                                                        
```
