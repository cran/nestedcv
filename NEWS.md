News
=====

# nestedcv 0.2.1
###### 15/06/2022

* Parallelisation on windows added
* hsstan model has been added (Athina Spiliopoulou)
* outer_folds can be specified for consistent model comparisons
* Checks on x, y added
* NA handling
* summary and print methods
* Implemented LOOCV
* Collinearity filter
* Implement lm and glm as models in outercv()
* Runnable examples have been added throughout

# nestedcv 0.0.9100
###### 02/03/2022

* Major update to include nestedcv.train function which adds nested CV to the 
`train` function of `caret`
* Note passing of extra arguments to filter functions specified by `filterFUN`
is no longer done through `...` but with a list of arguments passed through a
new argument `filter_options`.

# nestedcv 0.0.9003
###### 02/03/2022

* Initial build of nestedcv
* Added outercv.rf function for measuring performance of rf
* Added cv.rf for tuning mtry parameter
* Added plot_caret for plotting caret objects with error bars on the tuning 
metric
