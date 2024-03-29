% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precision_recall.R
\name{prc}
\alias{prc}
\alias{prc.default}
\alias{prc.data.frame}
\alias{prc.nestcv.glmnet}
\alias{prc.nestcv.train}
\alias{prc.nestcv.SuperLearner}
\alias{prc.outercv}
\alias{prc.repeatcv}
\title{Build precision-recall curve}
\usage{
prc(...)

\method{prc}{default}(response, predictor, positive = 2, ...)

\method{prc}{data.frame}(output, ...)

\method{prc}{nestcv.glmnet}(object, ...)

\method{prc}{nestcv.train}(object, ...)

\method{prc}{nestcv.SuperLearner}(object, ...)

\method{prc}{outercv}(object, ...)

\method{prc}{repeatcv}(object, ...)
}
\arguments{
\item{...}{other arguments}

\item{response}{binary factor vector of response of default order controls,
cases.}

\item{predictor}{numeric vector of probabilities}

\item{positive}{Either an integer 1 or 2 for the level of response factor
considered to be 'positive' or 'relevant', or a character value for that
factor.}

\item{output}{data.frame with columns \code{testy} containing observed response
from test folds, and \code{predyp} predicted probabilities for classification}

\item{object}{a 'nestcv.glmnet', 'nestcv.train', 'nestcv.SuperLearn',
'outercv' or 'repeatcv' S3 class results object.}
}
\value{
An object of S3 class 'prc' containing the following fields:
\item{recall}{vector of recall values}
\item{precision}{vector of precision values}
\item{auc}{area under precision-recall curve value using trapezoid method}
\item{baseline}{baseline precision value}
}
\description{
Builds a precision-recall curve for a 'nestedcv' model using \code{prediction()}
and \code{performance()} functions from the ROCR package and returns an object of
class 'prc' for plotting.
}
\examples{
\donttest{
library(mlbench)
data(Sonar)
y <- Sonar$Class
x <- Sonar[, -61]

fit1 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1, cv.cores = 2)

fit1$prc <- prc(fit1)  # calculate precision-recall curve
fit1$prc$auc  # precision-recall AUC value

fit2 <- nestcv.train(y, x, method = "gbm", cv.cores = 2)
fit2$prc <- prc(fit2)
fit2$prc$auc

plot(fit1$prc, ylim = c(0, 1))
lines(fit2$prc, col = "red")

res <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1) |>
  repeatcv(n = 4, rep.cores = 2)

res$prc <- prc(res)  # precision-recall curve on repeated predictions
plot(res$prc)
}
}
