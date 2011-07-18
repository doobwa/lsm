\name{lsm}
\alias{lsm}
\title{Gibbs sampling for the latent set model from DuBois, et al 2011.}
\usage{
  lsm(y, K, priors, niter = 20, zcode = 0, verbose = FALSE)
}
\arguments{
  \item{y}{binary Matrix of observed two-mode data (e.g.
  events x individuals)}

  \item{K}{number of latent sets to learn}

  \item{priors}{named list: \describe{ \item{z}{a vector
  (length two) with the Beta priors on probabilities of
  nonzero entries in each column of W and Z}
  \item{epsilon}{a vector (length two) with the Beta priors
  on noise parameter epsilon} \item{theta}{a vector (length
  two) with the parameters of the lognormal prior}}}

  \item{niter}{number of iterations to perform}

  \item{verbose}{if true, diagnostics will be shown after
  each iteration}

  \item{zcode}{determines how 0's are encoded in the
  provided matrix. (Default's to 0.)}

  \item{dims}{names list with T the number of events, N the
  number of}
}
\value{
  lsm object with parameter estimates from the final
  iteration, parameter values from each iteration,
  loglikelihood at each iteration, and the observed data.
}
\description{
  Perform Gibbs sampling for the latent set model from
  DuBois, et al 2011.
}
\author{
  Chris DuBois
}

