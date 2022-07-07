## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(nestedcv)
library(pROC)

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("myles-lewis/nestedcv", auth_token = "API token...")
#  library(nestedcv)

## -----------------------------------------------------------------------------
## Example binary classification problem with P >> n
x <- matrix(rnorm(150 * 2e+04), 150, 2e+04)  # predictors
y <- factor(rbinom(150, 1, 0.5))  # binary response

## Partition data into 2/3 training set, 1/3 test set
trainSet <- caret::createDataPartition(y, p = 0.66, list = FALSE)

## t-test filter using whole test set
filt <- ttest_filter(y, x, nfilter = 100)
filx <- x[, filt]

## Train glmnet on training set only using filtered predictor matrix
library(glmnet)
fit <- cv.glmnet(filx[trainSet, ], y[trainSet], family = "binomial")

## Predict response on test set
predy <- predict(fit, newx = filx[-trainSet, ], s = "lambda.min", type = "class")
predy <- as.vector(predy)
predyp <- predict(fit, newx = filx[-trainSet, ], s = "lambda.min", type = "response")
predyp <- as.vector(predyp)
output <- data.frame(testy = y[-trainSet], predy = predy, predyp = predyp)

## Results on test set
## shows bias since univariate filtering was applied to whole dataset
predSummary(output)

## Nested CV
fit2 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 7:10 / 10,
                      filterFUN = ttest_filter,
                      filter_options = list(nfilter = 100))
fit2

testroc <- pROC::roc(output$testy, output$predyp, direction = "<", quiet = TRUE)
inroc <- innercv_roc(fit2)
plot(fit2$roc)
lines(inroc, col = 'blue')
lines(testroc, col = 'red')
legend('bottomright', legend = c("Nested CV", "Left-out inner CV folds", 
                                 "Test partition, non-nested filtering"), 
       col = c("black", "blue", "red"), lty = 1, lwd = 2, bty = "n")

## ---- out.width='75%', fig.align="center", echo=FALSE-------------------------
knitr::include_graphics("fig1.svg")

## ---- out.width='75%', fig.align="center", echo=FALSE-------------------------
knitr::include_graphics("fig2.svg")

## ----eval = FALSE-------------------------------------------------------------
#  # Raw RNA-Seq data for this example is located at:
#  # https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11611/
#  
#  # set up data
#  load("/../R4RA_270821.RData")
#  
#  index <- r4ra.meta$Outliers_Detected_On_PCA != "outlier" & r4ra.meta$Visit == 3
#            & !is.na(r4ra.meta$Visit)
#  metadata <- r4ra.meta[index, ]
#  dim(metadata)  # 133 individuals
#  
#  medians <- Rfast::rowMedians(as.matrix(r4ra.vst))
#  data <- t(as.matrix(r4ra.vst))
#  # remove low expressed genes
#  data <- data[index, medians > 6]
#  dim(data)  # 16254 genes
#  
#  # Rituximab cohort only
#  yrtx <- metadata$CDAI.response.status.V7[metadata$Randomised.medication == "Rituximab"]
#  yrtx <- factor(yrtx)
#  data.rtx <- data[metadata$Randomised.medication == "Rituximab", ]
#  
#  # no filter
#  res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx,
#                           family = "binomial", cv.cores = 8,
#                           alphaSet = seq(0.7, 1, 0.05))
#  res.rtx

## ----eval = FALSE-------------------------------------------------------------
#  # t-test filter
#  res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, filterFUN = ttest_filter,
#                           filter_options = list(nfilter = 300, p_cutoff = NULL),
#                           family = "binomial", cv.cores = 8,
#                           alphaSet = seq(0.7, 1, 0.05))
#  summary(res.rtx)

## ----eval = FALSE-------------------------------------------------------------
#  plot_alphas(res.rtx)
#  plot_lambdas(res.rtx)

## ---- out.width='100%', fig.align="center", echo=FALSE------------------------
knitr::include_graphics("plot_alpha_lam.png")

## ----eval = FALSE-------------------------------------------------------------
#  # Fold 1 line plot
#  plot(res.rtx$outer_result[[1]]$cvafit)
#  
#  # Scatter plot
#  plot(res.rtx$outer_result[[1]]$cvafit, type = 'p')
#  
#  # Number of non-zero coefficients
#  plot(res.rtx$outer_result[[1]]$cvafit, xaxis = 'nvar')

## ---- out.width='100%', fig.align="center", echo=FALSE------------------------
knitr::include_graphics("plot_cva.png")

## ----eval = FALSE-------------------------------------------------------------
#  # Outer CV ROC
#  plot(res.rtx$roc, main = "Outer fold ROC", font.main = 1, col = 'blue')
#  legend("bottomright", legend = paste0("AUC = ", signif(pROC::auc(res.rtx$roc), 3)), bty = 'n')
#  
#  # Inner CV ROC
#  rtx.inroc <- innercv_roc(res.rtx)
#  plot(rtx.inroc, main = "Inner fold ROC", font.main = 1, col = 'red')
#  legend("bottomright", legend = paste0("AUC = ", signif(pROC::auc(rtx.inroc), 3)), bty = 'n')

## ---- out.width='100%', fig.align="center", echo=FALSE------------------------
knitr::include_graphics("roc.png")

## ----eval = FALSE-------------------------------------------------------------
#  boxplot_model(res.rtx, data.rtx, ylab = "VST")

## ---- out.width='70%', fig.align="center", echo=FALSE-------------------------
knitr::include_graphics("boxplot.png")

## ----eval = FALSE-------------------------------------------------------------
#  # Outer LOOCV
#  res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = ttest_filter,
#                           filter_options = list(nfilter = 300, p_cutoff = NULL),
#                           outer_method = "LOOCV",
#                           family = "binomial", cv.cores = 8,
#                           alphaSet = seq(0.7, 1, 0.05))
#  summary(res.rtx)

## ----eval = FALSE-------------------------------------------------------------
#  # Random forest filter
#  res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0.5, filterFUN = rf_filter,
#                           filter_options = list(nfilter = 300),
#                           family = "binomial", cv.cores = 8,
#                           alphaSet = seq(0.7, 1, 0.05))
#  summary(res.rtx)
#  
#  # ReliefF algorithm filter
#  res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = relieff_filter,
#                           filter_options = list(nfilter = 300),
#                           family = "binomial", cv.cores = 8,
#                           alphaSet = seq(0.7, 1, 0.05))
#  summary(res.rtx)

## ----eval = FALSE-------------------------------------------------------------
#  filter <- function(y, x, ...) {}

## ----eval = FALSE-------------------------------------------------------------
#  # nested CV using caret
#  tg <- expand.grid(lambda = exp(seq(log(2e-3), log(1e0), length.out = 100)),
#                    alpha = seq(0.8, 1, 0.1))
#  ncv <- nestcv.train(y = yrtx, x = data.rtx,
#                 method = "glmnet",
#                 savePredictions = "final",
#                 filterFUN = ttest_filter, filter_options = list(nfilter = 300),
#                 tuneGrid = tg, cv.cores = 8)
#  ncv$summary
#  
#  # Plot ROC on outer folds
#  plot(ncv$roc)
#  
#  # Plot ROC on inner LO folds
#  inroc <- innercv_roc(ncv)
#  plot(inroc)
#  pROC::auc(inroc)
#  
#  # Extract coefficients of final fitted model
#  glmnet_coefs(ncv$final_fit$finalModel, s = ncv$finalTune$lambda)

## ----eval = FALSE-------------------------------------------------------------
#  # Example tuning plot for outer fold 1
#  plot(ncv$outer_result[[1]]$fit, xTrans = log)
#  
#  # ggplot2 version
#  ggplot(ncv$outer_result[[1]]$fit) +
#    scale_x_log10()

## ----eval = FALSE-------------------------------------------------------------
#  # for nestcv.glmnet object
#  preds <- predict(res.rtx, newdata = data.rtx, type = 'response')
#  
#  # for nestcv.train object
#  preds <- predict(ncv, newdata = data.rtx)

## -----------------------------------------------------------------------------

# specify options for running the code
nfolds <- 3

# specify number of cores for parallelising computation
# the product of cv.cores and mc.cores will be used in total
# number of cores for parallelising over CV folds
cv.cores <- 1
# number of cores for parallelising hsstan sampling
# to pass CRAN package checks we need to create object oldopt
# in the end we reset options to the old configuration
oldopt <- options(mc.cores = 2)

# load iris dataset and simulate a continuous outcome
data(iris)
dt <- iris[, 1:4]
colnames(dt) <- c("marker1", "marker2", "marker3", "marker4")
dt <- as.data.frame(apply(dt, 2, scale))
dt$outcome.cont <- -3 + 0.5 * dt$marker1 + 2 * dt$marker2 + rnorm(nrow(dt), 0, 2)

# unpenalised covariates: always retain in the prediction model
uvars <- "marker1"
# penalised covariates: coefficients are drawn from hierarchical shrinkage prior
pvars <- c("marker2", "marker3", "marker4") # penalised covariates
# run cross-validation with univariate filter and hsstan
res.cv.hsstan <- outercv(y = dt$outcome.cont, x = dt[, c(uvars, pvars)],
                         model = model.hsstan,
                         filterFUN = lm_filter,
                         filter_options = list(force_vars = uvars,
                                               nfilter = 2,
                                               p_cutoff = NULL,
                                               rsq_cutoff = 0.9),
                         n_outer_folds = nfolds, cv.cores = cv.cores,
                         unpenalized = uvars, warmup = 1000, iter = 2000)
# view prediction performance based on testing folds
summary(res.cv.hsstan)

# load hsstan package to examine the Bayesian model
library(hsstan)
sampler.stats(res.cv.hsstan$final_fit)
print(projsel(res.cv.hsstan$final_fit), digits = 4)  # adding marker2
options(oldopt)  # reset options

## -----------------------------------------------------------------------------
# sigmoid function
sigmoid <- function(x) {1 / (1 + exp(-x))}

# specify options for running the code
nfolds <- 3

# specify number of cores for parallelising computation
# the product of cv.cores and mc.cores will be used in total
# number of cores for parallelising over CV folds
cv.cores <- 1
# number of cores for parallelising hsstan sampling
# to pass CRAN package checks we need to create object oldopt
# in the end we reset options to the old configuration
oldopt <- options(mc.cores = 2)

# load iris dataset and create a binary outcome
set.seed(267)
data(iris)
dt <- iris[, 1:4]
colnames(dt) <- c("marker1", "marker2", "marker3", "marker4")
dt <- as.data.frame(apply(dt, 2, scale))
rownames(dt) <- paste0("sample", c(1:nrow(dt)))
dt$outcome.bin <- sigmoid(0.5 * dt$marker1 + 2 * dt$marker2) > runif(nrow(dt))
dt$outcome.bin <- factor(dt$outcome.bin)


# unpenalised covariates: always retain in the prediction model
uvars <- "marker1"
# penalised covariates: coefficients are drawn from hierarchical shrinkage prior
pvars <- c("marker2", "marker3", "marker4") # penalised covariates
# run cross-validation with univariate filter and hsstan
res.cv.hsstan <- outercv(y = dt$outcome.bin,
                         x = as.matrix(dt[, c(uvars, pvars)]),
                         model = model.hsstan,
                         filterFUN = ttest_filter,
                         filter_options = list(force_vars = uvars,
                                               nfilter = 2,
                                               p_cutoff = NULL,
                                               rsq_cutoff = 0.9),
                         n_outer_folds = nfolds, cv.cores = cv.cores,
                         unpenalized = uvars, warmup = 1000, iter = 2000)


# view prediction performance based on testing folds
summary(res.cv.hsstan)

# examine the Bayesian model
sampler.stats(res.cv.hsstan$final_fit)
print(projsel(res.cv.hsstan$final_fit), digits = 4)  # adding marker2
options(oldopt)  # reset options

