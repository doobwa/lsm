\name{post.predictive}
\alias{post.predictive}
\title{Compute the posterior predictive using samples from multiple chains (ie. using multiple lsm fit objects).}
\usage{post.predictive(fits, subs, keep)
}
\value{vector of predictions averaged over the provided iterations from possibly multiple chains}
\author{Chris DuBois}
\arguments{\item{fits}{a list of lsm fit objects}
\item{subs}{subscripts of the matrix where predictions are desired}
\item{keep}{which iterations from each Gibbs sampling chain to use for computing the posterior loglikelihood}
}

