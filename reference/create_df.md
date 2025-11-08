# Create a data frame of protein intensities

This function creates a data frame of protein intensities

## Usage

``` r
create_df(
  prot_groups,
  exp_design,
  input_type = "MaxQuant",
  data_type = "LFQ",
  filter_na = TRUE,
  filter_prot = TRUE,
  uniq_pep = 2,
  tech_reps = FALSE,
  zero_na = TRUE,
  log_tr = TRUE,
  base = 2
)
```

## Arguments

- prot_groups:

  File path to a proteinGroups.txt file produced by MaxQuant or a
  standard input file containing a quantitative matrix where the
  proteins or protein groups are indicated by rows and the samples by
  columns.

- exp_design:

  File path to a text file containing the experimental design.

- input_type:

  Type of input file indicated by `prot_groups`. Available options are:
  "MaxQuant", if a proteinGroups.txt file is used, or "standard" if a
  standard input file is used. Default is "MaxQuant."

- data_type:

  Type of sample protein intensity data columns to use from the
  proteinGroups.txt file. Some available options are "LFQ", "iBAQ",
  "Intensity". Default is "LFQ." User-defined prefixes in the
  proteinGroups.txt file are also allowed. The `data_type` argument is
  case-sensitive, and only applies when `input_type = "MaxQuant"`.

- filter_na:

  Logical. If `TRUE`(default), filters out empty rows and columns from
  the data frame.

- filter_prot:

  Logical. If `TRUE` (default), filters out reverse proteins, proteins
  only identified by site, potential contaminants, and proteins
  identified with less than the minimum number of unique peptides
  indicated by `uniq_pep`. Only applies when `input_type = "MaxQuant"`.

- uniq_pep:

  Numerical. Proteins that are identified by this number or fewer number
  of unique peptides are filtered out (default is 2).Only applies when
  `input_type = "MaxQuant"`.

- tech_reps:

  Logical. Indicate as `TRUE` if technical replicates are present in the
  data. Default is `FALSE`.

- zero_na:

  Logical. If `TRUE` (default), zeros are considered missing values and
  replaced with NAs.

- log_tr:

  Logical. If `TRUE` (default), intensity values are log transformed to
  the base indicated by `base`.

- base:

  Numerical. Logarithm base. Default is 2.

## Value

A `raw_df` object which is a data frame containing protein intensities.
Proteins or protein groups are indicated by rows and samples by columns.

## Details

- It then reads in the expDesign.txt file provided as `exp_design` and
  extracts relevant information from it to add to the data frame. an
  example of the expDesign.txt is provided here:
  <https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt>.

- First, empty rows and columns are removed from the data frame.

- Next, if a proteinGroups.txt file is used, it filters out reverse
  proteins, proteins that were only identified by site, and potential
  contaminants.Then it removes proteins identified with less than the
  number of unique peptides indicated by `uniq_pep` from the data frame.

- Next, it extracts the intensity columns indicated by `data type` and
  the selected protein rows from the data frame.

- Converts missing values (zeros) to NAs.

- Finally, the function log transforms the intensity values.

## Author

Chathurani Ranathunge

## Examples

``` r
# \donttest{

### Using a proteinGroups.txt file produced by MaxQuant as input.
## Generate a raw_df object with default settings. No technical replicates.
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
  input_type = "MaxQuant"
)
#> 0 empty row(s) removed.
#> 0 empty column(s) removed.
#> 80 protein(s) (rows) only identified by site removed.
#> 65 reverse protein(s) (rows) removed.
#> 42 protein potential contaminant(s) (rows) removed.
#> 1923 protein(s) identified by 2 or fewer unique peptides removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.

## Data containing technical replicates
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg2.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed2.txt",
  input_type = "MaxQuant",
  tech_reps = TRUE
)
#> 0 empty row(s) removed.
#> 1 empty column(s) removed.
#> 12 reverse protein(s) (rows) removed.
#> 29 protein contaminant(s) (rows) removed.
#> 188 protein(s) identified by 2 or fewer unique peptides removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.

## Alter the number of unique peptides needed to retain a protein
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
  input_type = "MaxQuant",
  uniq_pep = 1
)
#> 0 empty row(s) removed.
#> 0 empty column(s) removed.
#> 80 protein(s) (rows) only identified by site removed.
#> 65 reverse protein(s) (rows) removed.
#> 42 protein potential contaminant(s) (rows) removed.
#> 961 protein(s) identified by 1 or fewer unique peptides removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.

## Use "iBAQ" values instead of "LFQ" values
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
  input_type = "MaxQuant",
  data_type = "iBAQ"
)
#> 0 empty row(s) removed.
#> 0 empty column(s) removed.
#> 80 protein(s) (rows) only identified by site removed.
#> 65 reverse protein(s) (rows) removed.
#> 42 protein potential contaminant(s) (rows) removed.
#> 1923 protein(s) identified by 2 or fewer unique peptides removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.

### Using a universal standard input file instead of MaxQuant output.
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/st.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
  input_type = "standard"
)
#> 0 empty row(s) removed.
#> 0 empty column(s) removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.
# }
```
