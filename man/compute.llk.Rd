\name{compute.llk}
\alias{compute.llk}
\title{Compute loglikelihood of data using the given parameters}
\usage{compute.llk(y, params, wsum="all", zcode=0, v)
}
\description{One of the main advantages of the link function used for this
latent set model is that we only need to compute predictions at
the elements (i,j) where person j is a member of at least one of
the active sets (ie.   {(i,j): sum_k w{ik}z_{jk} > 0}).  For all
other elements we know the prediction is 1-exp(epsilon).  When
computing the loglikelihood, then, the sum over all these other
elements is quickly available just by using the row (or column)
sums of the observed data, y, and the dimnesions of y.}
\details{I have yet to find a way for a Matrix object to default to NA
rather than 0.  In the case where the majority of the matrix is
unobserved, you may want to represent NA's using 0's.  This can be
done by setting zcode to -1.

##'}
\value{loglikelihood of the data given the provided model parameters}
\author{Chris DuBois}
\arguments{\item{y}{Matrix of observed data}
\item{params}{list of model parameter values}
\item{wsum}{all if the likelihood of all entries is desired, row: loglikelihoods of each row of y, col: loglikelihoods of each column of y}
\item{zcode}{determines what value represents 0.  0: 0's.}
\item{v}{Matrix of W x (Z * Theta).  Provide this if it's already known to save on computation time.}
}

