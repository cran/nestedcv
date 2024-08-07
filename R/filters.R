# Filters to reduce number of predictors in nested CV

#' Univariate filters
#'
#' A selection of simple univariate filters using t-test, Wilcoxon test, one-way
#' ANOVA or correlation (Pearson or Spearman) for ranking variables. These
#' filters are designed for speed. `ttest_filter` uses the `Rfast` package,
#' `wilcoxon_filter` (Mann-Whitney) test uses
#' [matrixTests::row_wilcoxon_twosample], `anova_filter` uses
#' [matrixTests::col_oneway_welch] (Welch's F-test) from the `matrixTests`
#' package. Can be applied to all or a subset of predictors. For mixed datasets
#' (combined continuous & categorical) see [stat_filter()]
#'
#' @param y Response vector
#' @param x Matrix or dataframe of predictors
#' @param force_vars Vector of column names within `x` which are always retained
#'   in the model (i.e. not filtered). Default `NULL` means all predictors will
#'   be passed to `filterFUN`.
#' @param nfilter Number of predictors to return. If `NULL` all predictors with
#'   p-values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param rsq_cutoff r^2 cutoff for removing predictors due to collinearity.
#'   Default `NULL` means no collinearity filtering. Predictors are ranked based
#'   on t-test. If 2 or more predictors are collinear, the first ranked
#'   predictor by t-test is retained, while the other collinear predictors are
#'   removed. See [collinear()].
#' @param type Type of vector returned. Default "index" returns indices, "names"
#'   returns predictor names, "full" returns a matrix of p values.
#' @param keep_factors Logical affecting factors with 3 or more levels.
#'   Dataframes are coerced to a matrix using [data.matrix]. Binary
#'   factors are converted to numeric values 0/1 and analysed as such. If
#'   `keep_factors` is `TRUE` (the default), factors with 3 or more levels are
#'   not filtered and are retained. If `keep_factors` is `FALSE`, they are
#'   removed.
#' @param exact Logical whether exact or approximate p-value is calculated.
#'   Default is `FALSE` for speed.
#' @param method Type of correlation, either "pearson" or "spearman".
#' @param ... optional arguments, including `rsq_method` passed to [collinear()]
#'   or arguments passed to [matrixTests::row_wilcoxon_twosample] in
#'   [wilcoxon_filter()].
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters in order
#'   of t-test p-value. If `type` is `"full"` full output from
#'   [Rfast::ttests] is returned.
#' @seealso [lm_filter()] [stat_filter()]
#' @examples
#' ## sigmoid function
#' sigmoid <- function(x) {1 / (1 + exp(-x))}
#' 
#' ## load iris dataset and simulate a binary outcome
#' data(iris)
#' dt <- iris[, 1:4]
#' colnames(dt) <- c("marker1", "marker2", "marker3", "marker4")
#' dt <- as.data.frame(apply(dt, 2, scale))
#' y2 <- sigmoid(0.5 * dt$marker1 + 2 * dt$marker2) > runif(nrow(dt))
#' y2 <- factor(y2, labels = c("C1", "C2"))
#' 
#' ttest_filter(y2, dt)  # returns index of filtered predictors
#' ttest_filter(y2, dt, type = "name")  # shows names of predictors
#' ttest_filter(y2, dt, type = "full")  # full results table
#'
#' data(iris)
#' dt <- iris[, 1:4]
#' y3 <- iris[, 5]
#' anova_filter(y3, dt)  # returns index of filtered predictors
#' anova_filter(y3, dt, type = "full")  # shows names of predictors
#' anova_filter(y3, dt, type = "name")  # full results table
#'
#' @importFrom Rfast ttests
#' @export
#'
ttest_filter <- function(y,
                         x,
                         force_vars = NULL,
                         nfilter = NULL,
                         p_cutoff = 0.05,
                         rsq_cutoff = NULL,
                         type = c("index", "names", "full"),
                         keep_factors = TRUE,
                         ...) {
  type <- match.arg(type)
  y <- factor(y)
  indx1 <- as.numeric(y) == 1
  indx2 <- as.numeric(y) == 2
  factor_ind <- which_factor(x)
  if (is.data.frame(x)) x <- data.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- seq_len(ncol(x))
  res <- Rfast::ttests(x[indx1, ], x[indx2, ])
  rownames(res) <- colnames(x)
  if (type == "full") {
    if (length(factor_ind) == 0) {
      return(res)
    } else return(res[-factor_ind, ])
  }
  filter_end(res[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff, type, keep_factors, 
             factor_ind, ...)
}


filter_end <- function(pval, x, force_vars, nfilter, p_cutoff, rsq_cutoff,
                          type, keep_factors, factor_ind, ...) {
  if (keep_factors) {
    force_vars <- unique(c(force_vars, colnames(x)[factor_ind]))
  }
  check_vars <- which(!colnames(x) %in% force_vars)
  outp <- pval[check_vars]
  outorder <- order(outp)
  out <- check_vars[outorder]
  outp <- outp[outorder]
  if (!is.null(p_cutoff)) out <- out[outp < p_cutoff]
  if (!is.null(rsq_cutoff)) {
    co <- collinear(x[, out], rsq_cutoff = rsq_cutoff, ...)
    if (length(co) > 0) out <- out[-co]
  }
  if (!keep_factors && !is.null(factor_ind)) out <- out[!out %in% factor_ind]
  if (!is.null(nfilter) && length(out) > nfilter) out <- out[1:nfilter]
  if (length(out) == 0) stop("No predictors left after filtering")
  out <- c(out, which(colnames(x) %in% force_vars))
  out <- out[!is.na(out)]
  switch(type,
         index = out, names = colnames(x)[out])
}


#' @rdname ttest_filter
#' @importFrom matrixTests col_oneway_welch
#' @export
#' 
anova_filter <- function(y,
                         x,
                         force_vars = NULL,
                         nfilter = NULL,
                         p_cutoff = 0.05,
                         rsq_cutoff = NULL,
                         type = c("index", "names", "full"),
                         keep_factors = TRUE,
                         ...) {
  type <- match.arg(type)
  y <- factor(y)
  factor_ind <- which_factor(x)
  if (is.data.frame(x)) x <- data.matrix(x)
  res <- matrixTests::col_oneway_welch(x, y)
  rownames(res) <- colnames(x)
  if (type == "full") {
    if (length(factor_ind) == 0) {
      return(res)
    } else return(res[-factor_ind, ])
  }
  filter_end(res[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff, type, keep_factors, 
             factor_ind, ...)
}


#' @rdname ttest_filter
#' @importFrom matrixTests row_wilcoxon_twosample
#' @export
#' 
wilcoxon_filter <- function(y,
                            x,
                            force_vars = NULL,
                            nfilter = NULL,
                            p_cutoff = 0.05,
                            rsq_cutoff = NULL,
                            type = c("index", "names", "full"),
                            exact = FALSE,
                            keep_factors = TRUE,
                            ...) {
  type <- match.arg(type)
  y <- factor(y)
  indx1 <- as.numeric(y) == 1
  indx2 <- as.numeric(y) == 2
  factor_ind <- which_factor(x)
  if (is.data.frame(x)) x <- data.matrix(x)
  res <- suppressWarnings(
    matrixTests::row_wilcoxon_twosample(t(x[indx1, ]), t(x[indx2, ]),
                                        exact = exact, ...)
  )
  if (type == "full") {
    if (length(factor_ind) == 0) {
      return(res)
    } else return(res[-factor_ind, ])
  }
  filter_end(res[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff, type,
             keep_factors, factor_ind)
}


#' Correlation between a vector and a matrix
#'
#' Fast Pearson/Spearman correlation where `y` is vector, `x` is matrix, adapted 
#' from [stats::cor.test].
#' 
#' @param y Numerical vector
#' @param x Matrix
#' @param method Type of correlation, either "pearson" or "spearman".
#' @param use Optional character string giving a method for computing 
#' covariances in the presence of missing values. See [cor]
#' @details For speed, p-values for Spearman's test are computed by 
#' asymptotic t approximation, equivalent to [cor.test] with `exact = FALSE`.
#' @return Matrix with columns containing the correlation statistic, either
#'   Pearson r or Spearman rho, and p-values for each column of `x` correlated
#'   against vector `y`
#' @importFrom stats complete.cases cor pt
#' @export
#' 
correls2 <- function(y, x,
                     method = "pearson",
                     use = "complete.obs") {
  ok <- complete.cases(x, y)
  x <- x[ok, ]
  y <- y[ok]
  n <- length(y)
  df <- n - 2
  if (method == "spearman") {
    r <- as.vector( cor(Rfast::Rank(y), Rfast::colRanks(x)) )
    q <- (n^3 - n) * (1 - r) / 6  # AS 89 algorithm
    den <- (n*(n^2-1))/6
    r <- 1 - q/den
    pval <- 2 * pt(abs(r) / sqrt((1 - r^2)/df), df, lower.tail=FALSE)
    rname <- 'rho'
  } else {
    r <- suppressWarnings(as.vector(cor(y, x, use = use)))
    STATISTIC <- sqrt(df) * r / sqrt(1 - r^2)
    pval <- 2 * pt(abs(STATISTIC), df, lower.tail = FALSE)  # two-tailed
    rname <- 'r'
  }
  res <- cbind(r, pval)
  colnames(res) <- c(rname, 'pvalue')
  rownames(res) <- colnames(x)
  res
}


#' @rdname ttest_filter
#' @export
#' 
correl_filter <- function(y,
                          x,
                          method = "pearson",
                          force_vars = NULL,
                          nfilter = NULL,
                          p_cutoff = 0.05,
                          rsq_cutoff = NULL,
                          type = c("index", "names", "full"),
                          keep_factors = TRUE,
                          ...) {
  type <- match.arg(type)
  factor_ind <- which_factor(x)
  if (is.data.frame(x)) x <- data.matrix(x)
  res <- correls2(y, x, method = method, ...)
  if (type == "full") {
    if (length(factor_ind) == 0) {
      return(res)
    } else return(res[-factor_ind, ])
  }
  filter_end(res[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff, type,
             keep_factors, factor_ind)
}


#' Random forest filter
#' 
#' Fits a random forest model and ranks variables by variable importance. 
#' 
#' @param y Response vector
#' @param x Matrix or dataframe of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param ntree Number of trees to grow. See [randomForest::randomForest].
#' @param mtry Number of predictors randomly sampled as candidates at each 
#' split. See [randomForest::randomForest].
#' @param ... Optional arguments passed to [randomForest::randomForest].
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` a named vector of variable importance is returned.
#' @details
#' This filter uses the `randomForest()` function from the `randomForest`
#' package. Variable importance is calculated using the
#' [randomForest::importance] function, specifying type 1 = mean decrease in
#' accuracy. See [randomForest::importance].
#' @export
#' 
rf_filter <- function(y, x, nfilter = NULL,
                      type = c("index", "names", "full"),
                      ntree = 1000,
                      mtry = ncol(x) * 0.2,
                      ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' must be installed to use this filter",
         call. = FALSE)
  }
  type <- match.arg(type)
  fit <- randomForest::randomForest(x, y, importance = TRUE,
                                    ntree = ntree, mtry = mtry, ...)
  vi <- as.vector(randomForest::importance(fit, type = 1))
  names(vi) <- if (type == "index") 1:ncol(x) else colnames(x)
  if (type == "full") return(vi)
  vi <- sort(vi, decreasing = TRUE)
  vi <- vi[vi != 0]
  if (!is.null(nfilter)) vi <- vi[1:min(nfilter, length(vi))]
  out <- names(vi)
  if (type == "index") out <- as.integer(out)
  out
}


#' Random forest ranger filter
#' 
#' Fits a random forest model via the `ranger` package and ranks variables by
#' variable importance.
#' 
#' @param y Response vector
#' @param x Matrix or dataframe of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param num.trees Number of trees to grow. See [ranger::ranger].
#' @param mtry Number of predictors randomly sampled as candidates at each 
#' split. See [ranger::ranger].
#' @param ... Optional arguments passed to [ranger::ranger].
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` a named vector of variable importance is returned.
#' @details
#' This filter uses the `ranger()` function from the `ranger` package. Variable
#' importance is calculated using mean decrease in gini impurity.
#' @export
#' 
ranger_filter <- function(y, x, nfilter = NULL,
                          type = c("index", "names", "full"),
                          num.trees = 1000,
                          mtry = ncol(x) * 0.2,
                          ...) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' must be installed to use this filter",
         call. = FALSE)
  }
  type <- match.arg(type)
  fit <- ranger::ranger(x = x, y = y, importance = "impurity",
                        num.trees = num.trees, mtry = mtry,
                        write.forest = FALSE,
                        verbose = FALSE,
                        num.threads = 1, ...)
  vi <- fit$variable.importance
  names(vi) <- if (type == "index") 1:ncol(x) else colnames(x)
  if (type == "full") return(vi)
  vi <- sort(vi, decreasing = TRUE)
  vi <- vi[vi != 0]
  if (!is.null(nfilter)) vi <- vi[1:min(nfilter, length(vi))]
  out <- names(vi)
  if (type == "index") out <- as.integer(out)
  out
}


#' ReliefF filter
#' 
#' Uses ReliefF algorithm from the CORElearn package to rank predictors in order 
#' of importance.
#' 
#' @param y Response vector
#' @param x Matrix or dataframe of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param estimator Type of algorithm used, see [CORElearn::attrEval]
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param ... Other arguments passed to [CORElearn::attrEval]
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` a named vector of variable importance is returned.
#' @seealso [CORElearn::attrEval()]
#' @export
#'
relieff_filter <- function(y, x, nfilter = NULL, 
                           estimator = "ReliefFequalK",
                           type = c("index", "names", "full"), ...) {
  if (!requireNamespace("CORElearn", quietly = TRUE)) {
    stop("Package 'CORElearn' must be installed to use this filter",
         call. = FALSE)
  }
  type <- match.arg(type)
  df <- as.data.frame(x)
  df$y <- y
  ref <- CORElearn::attrEval('y', df, estimator = estimator, ...)
  names(ref) <- if (type == "index") 1:ncol(x) else colnames(x)
  if (type == "full") return(ref)
  ref <- sort(ref, decreasing = TRUE)
  if (!is.null(nfilter)) ref <- ref[1:min(nfilter, length(ref))]
  out <- names(ref)
  if (type == "index") out <- as.integer(out)
  out
}


#' Combo filter
#' 
#' Filter combining univariate (t-test or anova) filtering and reliefF filtering 
#' in equal measure.
#' 
#' @param y Response vector
#' @param x Matrix or dataframe of predictors
#' @param nfilter Number of predictors to return, using 1/2 from `ttest_filter` 
#' or `anova_filter` and 1/2 from `relieff_filter`. Since `unique` is applied, 
#' the final number returned may be less than `nfilter`.
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns full output.
#' @param ... Optional arguments passed via [relieff_filter] to 
#' [CORElearn::attrEval]
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If `type` 
#' is `"full"` a list containing full outputs from either [ttest_filter] or 
#' [anova_filter] and [relieff_filter] is returned.
#' @export
#' 
combo_filter <- function(y, x,
                         nfilter,
                         type = c("index", "names", "full"), ...) {
  uni_set <- if (nlevels(y) == 2) {ttest_filter(y, x, nfilter, type = type)
  } else anova_filter(y, x, nfilter, type = type)
  relf_set <- relieff_filter(y, x, nfilter, type = type, ...)
  if (type == "full") {
    return(list(unifilt = uni_set, relieff_filter = relf_set))
  }
  n <- round(nfilter / 2)
  unique(c(uni_set[1:min(n, length(uni_set))],
           relf_set[1:min(n, length(relf_set))]))
}


#' glmnet filter
#'
#' Filter using sparsity of elastic net regression using glmnet to calculate
#' variable importance.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param family Either a character string representing one of the built-in
#'   families, or else a `glm()` family object. See [glmnet::glmnet()]. If not
#'   specified, the function tries to set this automatically to one of either
#'   "gaussian", "binomial" or "multinomial".
#' @param force_vars Vector of column names `x` which have no shrinkage and are
#'   always included in the model.
#' @param nfilter Number of predictors to return
#' @param method String indicating method of determining variable importance.
#'   "mean" (the default) uses the mean absolute coefficients across the range
#'   of lambdas; "nonzero" counts the number of times variables are retained in
#'   the model across all values of lambda.
#' @param type Type of vector returned. Default "index" returns indices, "names"
#'   returns predictor names, "full" returns full output.
#' @param ... Other arguments passed to [glmnet::glmnet]
#' @details The glmnet elastic net mixing parameter alpha can be varied to
#'   include a larger number of predictors. Default alpha = 1 is pure LASSO,
#'   resulting in greatest sparsity, while alpha = 0 is pure ridge regression,
#'   retaining all predictors in the regression model. Note, the `family`
#'   argument is commonly needed, see [glmnet::glmnet].
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters. If
#'   `type` is `"full"` a named vector of variable importance is returned.
#' @seealso [glmnet::glmnet]
#' @importFrom glmnet glmnet
#' @export
#' 
glmnet_filter <- function(y,
                          x,
                          family = NULL,
                          force_vars = NULL,
                          nfilter = NULL,
                          method = c("mean", "nonzero"),
                          type = c("index", "names", "full"),
                          ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  if (is.null(family)) {
    family <- if (is.factor(y) | is.character(y)) {
      if (nlevels(factor(y)) == 2) "binomial" else "multinomial"
    } else "gaussian"
  }
  penalty.factor <- rep(1, ncol(x))
  if (!is.null(force_vars)) {
    keep <- which(colnames(x) %in% force_vars)
    penalty.factor[keep] <- 0
  } else keep <- NULL
  fit <- glmnet(x, y, family = family, penalty.factor = penalty.factor, ...)
  cf <- as.matrix(coef(fit))
  if (method == "mean") {
    cf <- abs(cf)
    out <- rowMeans(cf)  # mean abs coefs
  } else {
    cf <- cf != 0
    out <- rowSums(cf)  # number of nonzero coefs
  }
  out <- out[-1]  # remove intercept
  if (type == "full") return(out)
  if (type == "index") names(out) <- 1:ncol(x)
  if (!is.null(force_vars)) {
    out <- out[c(keep, order(out[-keep], decreasing = TRUE))]
  } else out <- sort(out, decreasing = TRUE)
  out <- out[out != 0]
  if (!is.null(nfilter)) out <- out[1:min(nfilter, length(out))]
  out <- names(out)
  if (length(out) == 0) stop("No predictors selected")
  if (type == "index") out <- as.integer(out)
  out
}


# Code modified from caret::findCorrelation
# https://github.com/topepo/caret/blob/master/pkg/caret/R/findCorrelation.R

#' Filter to reduce collinearity in predictors
#'
#' This function identifies predictors with r^2 above a given cut-off and
#' produces an index of predictors to be removed. The function takes a matrix or
#' data.frame of predictors, and the columns need to be ordered in terms of
#' importance - first column of any pair that are correlated is retained and
#' subsequent columns which correlate above the cut-off are flagged for removal.
#'
#' @param x A matrix or data.frame of values. The order of columns is used to
#'   determine which columns to retain, so the columns in `x` should be sorted
#'   with the most important columns first.
#' @param rsq_cutoff Value of cut-off for r-squared
#' @param rsq_method character string indicating which correlation coefficient
#'   is to be computed. One of "pearson" (default), "kendall", or "spearman".
#'   See [cor()].
#' @param verbose Boolean whether to print details
#' @return Integer vector of the indices of columns in `x` to remove due to
#'   collinearity
#' @export
#' 
collinear <- function(x, rsq_cutoff = 0.9, rsq_method = "pearson",
                      verbose = FALSE) {
  rsq <- cor(x, method = rsq_method)^2
  rsq[lower.tri(rsq, diag = TRUE)] <- NA
  combsAboveCutoff <- which(rsq > rsq_cutoff)
  colsToCheck <- ceiling(combsAboveCutoff / nrow(rsq))
  rowsToCheck <- combsAboveCutoff %% nrow(rsq)
  colsToDiscard <- colsToCheck > rowsToCheck
  rowsToDiscard <- !colsToDiscard
  if (any(rowsToDiscard)) warning("Unexpected rows to discard")
  if (verbose) {
    df <- data.frame(keep = rowsToCheck, remove = colsToCheck)
    df <- df[order(df$keep), ]
    remd <- NULL
    for (i in unique(df$keep)) {
      rem <- df$remove[df$keep %in% i]
      rem <- rem[!rem %in% remd]
      if (length(rem) > 0) {
        if (!i %in% df$remove) cat("Keep ") else cat ("Removed ")
        cat(paste0(colnames(x)[i], ", remove "))
        cat(paste(colnames(x)[rem], collapse = ", "))
        remd <- c(remd, rem)
        cat("\n")
      }
    }
  }
  deletecol <- colsToCheck
  deletecol <- unique(deletecol)
  deletecol
}

#' Linear model filter
#'
#' Linear models are fitted on each predictor, with inclusion of variable names
#' listed in `force_vars` in the model. Predictors are ranked by Akaike
#' information criteria (AIC) value, or can be filtered by the p-value on the
#' estimate of the coefficient for that predictor in its model.
#'
#' @param y Numeric or integer response vector
#' @param x Matrix of predictors. If `x` is a data.frame it will be turned into
#'   a matrix. But note that factors will be reduced to numeric values, but a
#'   full design matrix is not generated, so if factors have 3 or more levels,
#'   it is recommended to convert `x` into a design (model) matrix first.
#' @param force_vars Vector of column names `x` which are incorporated into the
#'   linear model.
#' @param nfilter Number of predictors to return. If `NULL` all predictors with
#'   p-values < `p_cutoff` are returned.
#' @param p_cutoff p-value cut-off. P-values are calculated by t-statistic on
#'   the estimated coefficient for the predictor being tested.
#' @param rsq_cutoff r^2 cutoff for removing predictors due to collinearity.
#'   Default `NULL` means no collinearity filtering. Predictors are ranked based
#'   on AIC from a linear model. If 2 or more predictors are collinear, the
#'   first ranked predictor by AIC is retained, while the other collinear
#'   predictors are removed. See [collinear()].
#' @param rsq_method character string indicating which correlation coefficient
#'   is to be computed. One of "pearson" (default), "kendall", or "spearman".
#'   See [collinear()].
#' @param type Type of vector returned. Default "index" returns indices, "names"
#'   returns predictor names, "full" returns a matrix of p values.
#' @param keep_factors Logical affecting factors with 3 or more levels.
#'   Dataframes are coerced to a matrix using [data.matrix]. Binary
#'   factors are converted to numeric values 0/1 and analysed as such. If
#'   `keep_factors` is `TRUE` (the default), factors with 3 or more levels are
#'   not filtered and are retained. If `keep_factors` is `FALSE`, they are
#'   removed.
#' @param method Integer determining linear model method. See
#'   [RcppEigen::fastLmPure()]
#' @param mc.cores Number of cores for parallelisation using
#'   [parallel::mclapply()].
#' @return Integer vector of indices of filtered parameters (`type = "index"`)
#'   or character vector of names (`type = "names"`) of filtered parameters in
#'   order of linear model AIC. Any variables in `force_vars` which are
#'   incorporated into all models are listed first. If `type = "full"` a matrix
#'   of AIC value, sigma (residual standard error, see [summary.lm]),
#'   coefficient, t-statistic and p-value for each tested predictor is returned.
#' @details
#' This filter is based on the model `y ~ xvar + force_vars` where `y` is the
#' response vector, `xvar` are variables in columns taken sequentially from `x`
#' and `force_vars` are optional covariates extracted from `x`. It uses
#' [RcppEigen::fastLmPure()] with `method = 0` as default since it is
#' rank-revealing. `method = 3` is significantly faster but can give errors in
#' estimation of p-value with variables of zero variance. The algorithm attempts
#' to detect these and set their stats to `NA`. `NA` in `x` are not tolerated.
#' 
#' Parallelisation is available via [mclapply()]. This is provided mainly for
#' the use case of the filter being used as standalone. Nesting parallelisation
#' inside of parallelised [nestcv.glmnet()] or [nestcv.train()] loops is not
#' recommended.
#' 
#' @export
#'
lm_filter <- function(y, x,
                      force_vars = NULL,
                      nfilter = NULL,
                      p_cutoff = 0.05,
                      rsq_cutoff = NULL,
                      rsq_method = "pearson",
                      type = c("index", "names", "full"),
                      keep_factors = TRUE,
                      method = 0L,
                      mc.cores = 1) {
  if (!requireNamespace("RcppEigen", quietly = TRUE)) {
    stop("Package 'RcppEigen' must be installed to use this filter",
         call. = FALSE)
  }
  type <- match.arg(type)
  ok <- !is.na(y)
  y <- y[ok]
  x <- x[ok, ]
  factor_ind <- which_factor(x)
  if (is.data.frame(x)) x <- data.matrix(x)
  check_vars <- colnames(x)[!colnames(x) %in% force_vars]
  colvars <- matrixStats::colVars(x[, check_vars])
  var0 <- colvars == 0  # fastLmPure incorrect if var(x) = 0
  startx <- matrix(rep(1, nrow(x) *2), ncol = 2)
  colnames(startx) <- c("(Intercept)", ".test")
  if (length(force_vars) > 0) {
    xset0 <- x[, force_vars, drop = FALSE]
    xset <- cbind(startx, xset0)
  } else xset <- startx
  res <- mclapply(check_vars, function(i) {
    xset[, 2] <- x[, i]
    fit <- RcppEigen::fastLmPure(xset, y, method = method)
    if (any(is.na(fit$coefficients))) return(rep(NA, 3))  # check for rank deficient
    rss <- sum(fit$residuals^2)
    c(rss, fit$coefficients[2], fit$se[2])
  }, mc.cores = mc.cores)
  res <- simplify2array(res)
  colnames(res) <- check_vars
  rss <- res[1,]
  cf <- res[2,]
  tval <- res[2,] / res[3,]
  n <- length(y)
  P <- ncol(xset)
  ## from stats::logLik.lm
  loglik <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(rss)))
  aic <- -2 * loglik + 2 * (P + 1)  ## from stats::AIC
  rdf <- n - P
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  pval <- 2*pt(abs(tval), rdf, lower.tail = FALSE)  ## from stats::summary.lm
  result <- cbind(aic, sigma, coef = cf, tval, pval)
  result[var0, ] <- NA
  if (type == "full") {
    if (length(factor_ind) == 0) {
      return(result)
    } else return(result[-factor_ind, ])
  }
  filter_end(pval,
             x, force_vars, nfilter, p_cutoff, rsq_cutoff, type, keep_factors, 
             factor_ind)
}
