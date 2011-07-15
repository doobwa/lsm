gibbs <- function(y,dims,priors,covars=NULL,niter=20,
                  verbose=TRUE,z.init="prior",zcode=0) {

  zs <- ws <- thetas <- list()
  llks <- rep(0,niter)

  # Initialize z and theta and sample w first.
  params <- list()
  phis <- rbeta(dims$K,priors$z[1],priors$z[2])
  z <- Matrix(0,dims$N,dims$K)
  for (k in 1:dims$K) {
    z[,k] <- 1*(runif(dims$N) < phis[k])
  }
  params$z <- z
  params$w <- Matrix(0, dims$T, dims$K)
  params$psi <- 1  
  params$theta <- Matrix(- initialize.theta(dims$N, dims$K))  # should be positive
  params$epsilon <- log(1 - priors$epsilon[1]/priors$epsilon[2])
  for (iter in 1:niter) {
    params <- sampleW(y,dims,params,priors,verbose,zcode)
    params <- sampleZ(y,dims,params,priors,verbose,zcode)
    params <- sampleTheta(y,dims,params,priors,verbose,zcode)
    params <- sampleEpsilon(y,dims,params,priors,verbose,zcode)

    # Save samples
    zs[[iter]] <- params$z
    ws[[iter]] <- params$w
    thetas[[iter]] <- params$theta

    fit <- list(y=y,params=params,llks=llks,zs=zs,ws=ws,thetas=thetas,dims=dims)
    class(fit) <- "lsm"

    # Display progress
    cat("\niter:",iter," llk:")
    if (verbose) {
      plotProgress(fit)
      llks[iter] <- compute.llk(y,params,wsum="all",zcode=zcode)
      cat(llks[iter],"|")
    }
  }
  return(fit)
}
predict.lsm <- function (fit,subs) {
  params <- fit$params
    if (nrow(subs) > 1e+06) 
        error("Too many predictions to make.")
    ws <- params$w[subs[, 1], ]
    zs <- params$z[subs[, 2], ]
    thetas <- params$theta[subs[, 2], ]
    v <- - ws * zs * thetas
    yhat <- 1 - exp(rowSums(v) + params$epsilon)
    return(yhat)
  }

post.predictive.list <- function (fits, subs, keep) 
{
  y <- fits[[1]]$y
  F <- length(fits)
  K <- ncol(fits[[1]]$params$z)
  yhats <- matrix(0, nrow(subs), length(keep) * F)
  for (f in 1:F) {
    fit <- fits[[f]]
    for (i in 1:length(keep)) {
      k <- keep[i]
      yhats[, (f - 1) * length(keep) + i] <- predict(fit,subs)
    }
  }
  yhat <- matrix(yhats,nrow(y),ncol(y))
  if (!is.null(rownames(y))) rownames(yhat) <- rownames(y)
  if (!is.null(colnames(y))) colnames(yhat) <- colnames(y)  
  return(yhat)
}
post.predictive <- function(object,...) UseMethod("post.predictive")

llk <- function(fit) {
  
}

compute.llk <- function(y,params,wsum="all",zcode=0,v=NULL) {
  # zero.code: the value in y that indicates a 0.  If zero.code = -1, then 0's are
  #            interpreted as NA's.
  w <- params$w
  z <- params$z
  epsilon <- params$epsilon
  theta <- params$theta
  psi <- params$psi
  
  if (zcode == 0) {

    # Get likelihood contribution for each actor
    sum2 <- function(x) sum(x,na.rm=TRUE)
    rowSums2 <- function(x) rowSums(x,na.rm=TRUE)
    colSums2 <- function(x) colSums(x,na.rm=TRUE)
    # Rowwise, columnwise counts of observed values
    rowcounts <- ncol(y) - rowSums(is.na(y))
    colcounts <- nrow(y) - colSums(is.na(y))

    # Compute exponent
    if (is.null(v)) {
      v <- (w * psi) %*% t(z * theta)
      v@x <- -v@x
    }

    # Indicates which elements are not zero: I(vij<0), I(yij=1,vij<0)
    v.nz <- (v < 0)
    v.nz[which(is.na(y))] <- NA  # This line fixed a bug!
    yv.nz <- (y * v.nz)

    # Only work on elements of matrix
    a <- b <- v
    a@x <- a@x + epsilon
    b@x <- log(exp(-a@x) - 1)
    lij <- yv.nz * b + v.nz * a

    # Respective sums
    llk <- switch(wsum,
                  "col" = colSums2(lij) +
                          log(exp(-epsilon) - 1) * (colSums2(y) - colSums2(yv.nz)) +
                          epsilon * (colcounts - colSums2(v.nz)),
                  "row" = rowSums2(lij) +
                          log(exp(-epsilon) - 1) * (rowSums2(y) - rowSums2(yv.nz)) +
                          epsilon * (rowcounts - rowSums2(v.nz)),
                  "all" = sum2(lij) +
                          sum2(log(exp(-epsilon) - 1) * (colSums2(y) - colSums2(yv.nz))) +
                          sum2(epsilon * (colcounts - colSums2(v.nz))))
  } else if (zcode==-2) {
    # Means 0's are coded -1, but there are no 0's
    # Get likelihood contribution for each actor
    y1 <- drop0(y == 1)

    # Compute exponent
    if (is.null(v)) {
      v <- -(w * psi) %*% t(z * theta)
    }

    v1 <- (v < 0)  # indicates sum_k w_ik z_jk theta_jk is nonzero
    y1v1 <- (y1 * v1)
    y1v0 <- (y1 - y1v1)

    # Only work on elements of matrix
    a <- b <- v
    a@x <- a@x + epsilon
    b@x <- log(1 - exp(a@x))
    lij <- y1v1 * b + y1v0 * log(1 - exp(epsilon))

    # Respective sums
    llk <- switch(wsum,
                  "col" = colSums(lij),
                  "row" = rowSums(lij),
                  "all" = sum(lij))

  } else {
    # Get likelihood contribution for each actor
    y1 <- drop0(y == 1)
    y0 <- drop0(y == -1)

    # Compute exponent
    if (is.null(v)) {
      v <- -(w * psi) %*% t(z * theta)
    }

    v1 <- (v < 0)  # indicates sum_k w_ik z_jk theta_jk is nonzero
    y1v1 <- (y1 * v1)
    y0v1 <- (y0 * v1)
    y1v0 <- drop0(y1 - y1v1)
    y0v0 <- drop0(y0 - y0v1)

    # Only work on elements of matrix
    a <- b <- v
    a@x <- a@x + epsilon
    b@x <- log(1 - exp(a@x))
    lij <- y1v1 * b + y0v1 * a + y1v0 * log(1 - exp(epsilon)) + y0v0 * epsilon

    # Respective sums
    llk <- switch(wsum,
                  "col" = colSums(lij),
                  "row" = rowSums(lij),
                  "all" = sum(lij))
  }
  return(llk)
}

sampleTheta <- function(y,dims,params,priors,verbose=FALSE,zcode=0) {
  if (verbose) cat("\nSampling theta")  
  w <- params$w
  z <- params$z
  ix <- which(z == 0)
  params$theta[ix] <- -initialize.theta(dims$N,dims$K)[ix]
  psi <- params$psi
  epsilon <- params$epsilon
  # For each latent set
  for (k in 1:dims$K) {
    if (verbose) cat(".")
    js <- which(z[,k] == 1)
    if (length(js) > 0) {
      theta.current <- theta.proposed <- params$theta

      # Get proposed thetas where z_jk=1
      dtheta <- rnorm(length(js))  # TODO: Better proposal distribution
      theta.proposed[js,k] <- theta.proposed[js,k] - dtheta #using minus to keep symmetry with other version of sampleTheta
      theta.proposed[which(theta.proposed < 0)] <- 0

      # Get vector of prob. accept
      params$theta <- theta.current
      lj.current <- compute.llk(y,params,"col",zcode)
      params$theta <- theta.proposed
      lj.propose <- compute.llk(y,params,"col",zcode)
      params$theta <- theta.current
      currentPrior <- log(dnorm(theta.current[js,k],priors$theta[1],priors$theta[2]))
      proposalPrior <- log(dnorm(theta.proposed[js,k],priors$theta[1],priors$theta[2]))
      r <- exp(lj.propose[js] - lj.current[js] + proposalPrior - currentPrior)

      # Update lj's to have lj.proposal where accepted
      for (i in 1:length(js)) {
        j <- js[i]
        if (runif(1) < r[i]) {
          params$theta[j,k] <- theta.proposed[j,k]
        } else {
          params$theta[j,k] <- theta.current[j,k]
        }
      }
    }
  }
  return(params)
}


sampleW <- function(y,dims,params,priors,verbose=FALSE,zcode=0) {
  # TODO: Put a different prior on the w's (ie. not the same prior as z's).
  
  if (verbose) cat("\nSampling w")

  # For each k, sample
  for (k in 1:dims$K) {
    if (verbose) cat(".")
    pz <- matrix(0,dims$T,2)
    for (a in 0:1) {
      params$w[,k] <- a
      q <- rowSums(params$w)
      r <- compute.llk(y,params,"row",zcode)
      pz[,a+1] <- r + lbeta(q + priors$z[1], dims$K - q + priors$z[2])
    }
    params$w[,k] <- samplelpmat(pz)
  }
  return(params)
}

sampleZ <- function(y,dims,params,priors,verbose=FALSE,zcode=0) {
  if (verbose) cat("\nSampling z")
  # Draw thetas where z=0 from prior
  theta <- - initialize.theta(dims$N,dims$K)
  iz <- which(params$z == 1)
  theta[iz] <- params$theta[iz]  # use current thetas where z==1

  # For each k, sample
  for (k in 1:dims$K) {
    if (verbose) cat(".")
    pz <- matrix(0,dims$N,2)
    for (a in 0:1) {
      params$z[,k] <- a
      q <- rowSums(params$z)
      r <- compute.llk(y,params,"col",zcode)
      pz[,a+1] <- r + lbeta(q + priors$z[1], dims$K - q + priors$z[2])
    }
    params$z[,k] <- samplelpmat(pz)
  }
  return(params)
}

sampleEpsilon <- function(y,dims,params,priors,verbose=FALSE,zcode=0) {
  if (verbose) cat("\nSampling epsilon")

  # Draw proposed epsilon from prior
  e1 <- params$epsilon
#  e2 <- log(1 - rbeta(1,priors$epsilon[1],priors$epsilon[2]))
  e2 <- 1-exp(-e1)
  e2 <- params$epsilon + rnorm(1,e1,.02)
  e2 <- ifelse(e2>0,e2,.0001)
  e2 <- log(1-e2)
  # Get acceptance probability
  params$epsilon <- e1
  llk.current <- compute.llk(y,params,"all",zcode)
  params$epsilon <- e2
  llk.propose <- compute.llk(y,params,"all",zcode)
  currentPrior <- log(pbeta(1 - exp(e1), priors$epsilon[1], priors$epsilon[2]))
  proposalPrior <- log(pbeta(1 - exp(e2), priors$epsilon[1], priors$epsilon[2]))
  r <- exp(llk.propose - llk.current + proposalPrior - currentPrior)

  # Accept/reject
  if (runif(1) < r) {
    params$epsilon <- e2
  } else {
    params$epsilon <- e1
  }
  return(params)
}
