\name{post.predictive}
\alias{post.predictive}
\usage{
  post.predictive(fits, subs, keep)
}
\arguments{
  \item{fits}{a list of lsm fit objects}

  \item{subs}{subscripts of the matrix where predictions
  are desired}

  \item{keep}{which iterations from each Gibbs sampling
  chain to use for computing the posterior loglikelihood}
}
\value{
  vector of predictions averaged over the provided
  iterations from possibly multiple chains
}
\author{
  Chris DuBois
}

