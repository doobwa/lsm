\name{predict.lsm}
\alias{predict.lsm}
\title{Make predictions for a fitted lsm model}
\usage{
  predict.lsm(fit, iter = length(fit$ws), subs)
}
\arguments{
  \item{fit}{an lsm fit object}

  \item{iter}{which iteration from Gibbs sampling to use
  for predictions.  Cannot be larger than the number used
  when the model was fit.  Defaults to the last iteration.}

  \item{subs}{subscripts of the orginal matrix where
  predictions are deisred}
}
\value{
  a vector of predictions (y.hat), one for each row of subs
}
\description{
  Uses the point estimates from the given iteration of
  Gibbs sampling to compute predictions at the provided
  elements.  The equation is: y.hat[i,j] = 1 - exp(sum_k
  w[i,k] * z[j,k] * thetas[j,k] + epsilon)
}
\author{
  Chris DuBois
}

