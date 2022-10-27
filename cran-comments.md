# promor 0.1.1

### Re-submission
In this re-submission:

* Replaced "https://www.biorxiv.org/content/early/2022/08/18/2022.08.17.503867" with "https://www.biorxiv.org/content/10.1101/2022.08.17.503867v1" in the inst/CITATION file
as suggested by the CRAN team.

### Submission

### R CMD check results - 10-26-2022
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

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

### Re-submission

### R CMD check results - 07-19-2022
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

In this re-submission:

* Added references to the DESCRIPTION file.
* Limited the use of `\dontrun` to an example that would otherwise write to the user's working directory* 
* Added `\donttest` to examples that are not executable in < 5 sec and unwrapped all the remaining examples* 
* Removed unnecessary information messages that were previously written to the console* 
* Provided seed as an argument to avoid setting seed to a specific number within functions* 
* Made changes in functions to write files to `tempdir()` when `file_path` is not specified.


### First submission

### R CMD check results - 07-08-2022
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
