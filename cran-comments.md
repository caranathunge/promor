# promor 0.1.0

## Re-submission
## R CMD check results - 07-19-2022
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

In this re-submission:

• Added references to the DESCRIPTION file.

• Limited the use of `\dontrun` to an example that would otherwise write to the user's working directory.

• Added `\donttest` to examples that are not executable in < 5 sec and unwrapped all the remaining examples.

• Removed unnecessary information messages that were previously written to the console.

• Provided seed as an argument to avoid setting seed to a specific number within functions.

• Made changes in functions to write files to `tempdir()` when `file_path` is not specified.


## First submission
## R CMD check results - 07-08-2022
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
