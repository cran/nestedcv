
#' Boruta filter
#'
#' Filter using Boruta algorithm.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param select Which type of features to retain. Options include "Confirmed"
#'   and/or "Tentative".
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param ... Other arguments passed to [Boruta::Boruta()]
#' @details
#' Boruta works differently from other filters in that it does not rank
#' variables by variable importance, but tries to determine relevant features
#' and divides features into Rejected, Tentative or Confirmed.
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters. If
#'   `type` is `"full"` full output from `Boruta` is returned.
#' 
#' @export

boruta_filter <- function(y, x, select = c('Confirmed', 'Tentative'),
                           type = c("index", "names", "full"), ...) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    stop("Package 'Boruta' must be installed", call. = FALSE)
  }
  type <- match.arg(type)
  ref <- Boruta::Boruta(x, y, ...)$finalDecision
  out <- which(ref %in% select)
  if (length(out) == 0) stop("No predictors left after filtering")
  switch(type,
         index = out,
         names = colnames(x)[out],
         full = ref)
}



#' Partial Least Squares filter
#'
#' Filter using coefficients from partial least squares (PLS) regression to
#' select optimal predictors.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param force_vars Vector of column names within `x` which are always retained
#'   in the model (i.e. not filtered). Default `NULL` means all predictors will
#'   be filtered.
#' @param nfilter Either a single value for the total number of predictors to
#'   return. Or a vector of length `ncomp` to manually return predictors from
#'   each PLS component.
#' @param ncomp the number of components to include in the PLS model.
#' @param scale_x Logical whether to scale predictors before fitting the PLS
#'   model. This is recommended.
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param ... Other arguments passed to [pls::plsr()]
#' @details
#' The best predictors may overlap between components, so if `nfilter` is
#' specified as a vector, the total number of unique predictors returned may be
#' variable.
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters. If
#'   `type` is `"full"` full output of coefficients from `plsr` is returned as a
#'   list for each model component ordered by highest absolute coefficient.
#' @export

pls_filter <- function(y, x,
                       force_vars = NULL,
                       nfilter,
                       ncomp = 5,
                       scale_x = TRUE,
                       type = c("index", "names", "full"), ...) {
  if (!requireNamespace("pls", quietly = TRUE)) {
    stop("Package 'pls' must be installed", call. = FALSE)
  }
  type <- match.arg(type)
  if (is.factor(y) && nlevels(y) > 2) stop("Classes > 2 not supported")
  y <- as.numeric(y)
  x0 <- x
  if (scale_x) {
    x <- scale(x)
    sd0 <- which(attr(x, "scaled:scale") == 0)
    if (length(sd0) > 0) x <- x[, -sd0]
  }
  fit <- pls::plsr(y ~ x, ncomp = ncomp, ...)
  cf <- fit$coefficients
  cf <- lapply(seq_len(ncomp), function(i) {
    cfi <- cf[,, i]
    cfi[order(abs(cfi), decreasing = TRUE)]
  })
  if (type == "full") return(cf)
  
  nfilter <- pmin(nfilter, ncol(x))
  if (length(nfilter) == 1) {
    # find sufficient vars from each comp
    topvars <- ""
    n <- ceiling(nfilter / ncomp)
    while (length(topvars) < nfilter) {
      topvars <- unique(unlist(lapply(seq_len(ncomp), function(i) {
        names(cf[[i]][seq_len(n)])
      })))
      n <- n +1
    }
    topvars <- topvars[seq_len(nfilter)]
  } else {
    # nfilter as vector
    if (length(nfilter) != ncomp) stop("nfilter is not the same length as ncomp")
    topvars <- unique(unlist(lapply(seq_len(ncomp), function(i) {
      names(cf[[i]][1:nfilter[i]])
    })))
  }
  
  topvars <- unique(c(topvars, force_vars))
  if (type == "names") return(topvars)
  which(colnames(x0) %in% topvars)
}

