% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.hsstan.R
\name{predict.hsstan}
\alias{predict.hsstan}
\title{Predict from hsstan model fitted within cross-validation}
\usage{
\method{predict}{hsstan}(object, newdata = NULL, type = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{hsstan}.}

\item{newdata}{Optional data frame containing the variables to use to
predict. If \code{NULL} (default), the model matrix is used. If specified, its
continuous variables should be standardized, since the model coefficients
are learnt on standardized data.}

\item{type}{Option for binary outcomes only. Default \code{NULL} will return a
class with the highest probability for each sample. If set to \code{probs}, it
will return the probabilities for outcome = 0 and for outcome = 1 for each
sample.}

\item{...}{Optional arguments passed to \code{hsstan::posterior_predict}}
}
\value{
For a binary outcome and type = \code{NULL}, a character vector with the
name of the class that has the highest probability for each sample.
For a binary outcome and type = \code{prob}, a 2-dimensional matrix with the
probability of class 0 and of class 1 for each sample.
For a continuous outcome a numeric vector with the predicted value for
each sample.
}
\description{
Draws from the posterior predictive distribution of the outcome.
}
\author{
Athina Spiliopoulou
}
